from atmos import atmosphere, R_e
from drag import KdCurve
from math import pi, sin, cos, atan2, acos

from num import RKF78, gss, bisect


def dec_to_dms(deg):
    sign = deg >= 0

    deg = abs(deg)
    h = int(deg)
    m = int((deg - h) * 60)
    s = int((deg - h - m / 60) * 3600)

    # return sign, h, m, s

    return "{:}{:03}°{:02}'{:02}\"".format("+" if sign else "-", h, m, s)


class Bullet:
    def __init__(self, name, mass, diam, form, Kd_curve):
        self.name = name
        self.W = mass
        self.D = diam
        self.i = form
        self.f_Kd = Kd_curve

        self.C = self.W / (self.i * self.D**2)  # ballistic coefficient

    def _ode_t(self, t, x, y, vx, vy):
        """
        t: time of flight
        x, y: cooridnate in 2-d geocentric coordinate system:
        vx, vy: component of velocity
        """

        v = (vx**2 + vy**2) ** 0.5
        r = (x**2 + y**2) ** 0.5
        h = r - R_e
        dx = vx
        dy = vy

        _, _, c, rho, g = self.f_env(h)

        M = v / c  # mach
        Kd = self.f_Kd.get(M)  # drag coefficient
        a = Kd * rho * v**2 / self.C
        dvx = -a * vx / v - g * x / r
        dvy = -a * vy / v - g * y / r

        return dx, dy, dvx, dvy

    def record_to_data(self, record, prettyprint=True):
        """convert a forward trajectory record into more sensible values:"""

        data = []
        for line in record:
            t, (x, y, vx, vy) = line
            r = (x**2 + y**2) ** 0.5
            h = r - R_e

            theta = -(atan2(y, x) - 0.5 * pi)

            """
            Geocentric angle measuring from shot start to shot spalsh
            shot
            start    shot
                +--_splash
                |  /
            R_e |θ/ R_e
                |/
                +
            """
            gr = theta * R_e
            v = (vx**2 + vy**2) ** 0.5

            """
            calculate the vector angle between r_vec and v_vec, which
            is the complementary angle of the shot-to-horizon
            """
            phi = 90 - acos((x * vx + y * vy) / (r * v)) * 180 / pi
            data.append((t, h, gr, v, phi))

        if prettyprint:
            from tabulate import tabulate

            print(
                tabulate(
                    data,
                    headers=("ToF", "H", "G.R.", "V.", "Angle"),
                )
            )

        return data

    def rangeTable(
        self,
        tol,
        vel,
        minR,
        maxR,
        deltaR,
        gunH=0,
        tgtH=0,
        env=atmosphere,
        t_max=1000,
        prettyprint=True,
    ):
        R = minR
        Rs = []
        while R <= maxR:
            Rs.append(R)
            R += deltaR

        lTrajs, hTrajs = self.inverse(tol, vel, Rs, gunH, tgtH, env, t_max)

        lTable, hTable = [], []

        for R, lTraj in zip(Rs, lTrajs):
            if lTraj is not None:
                elev, t, (x, y, vx, vy) = lTraj
                r = (x**2 + y**2) ** 0.5
                v = (vx**2 + vy**2) ** 0.5
                phi = 90 - acos((x * vx + y * vy) / (r * v)) * 180 / pi
                lTable.append((R, dec_to_dms(elev), v, dec_to_dms(phi), t))

        for R, hTraj in zip(Rs, hTrajs):
            if hTraj is not None:
                elev, t, (x, y, vx, vy) = hTraj
                r = (x**2 + y**2) ** 0.5
                v = (vx**2 + vy**2) ** 0.5
                phi = 90 - acos((x * vx + y * vy) / (r * v)) * 180 / pi
                hTable.append((R, dec_to_dms(elev), v, dec_to_dms(phi), t))

        if prettyprint:
            from tabulate import tabulate

            headers = (
                "Ground\nRange m",
                "Launch\nAngle",
                "Velocity\nm/s",
                "Impact\nAngle",
                "Time of\nFlight s",
            )

            print("Low")
            print(tabulate(lTable, headers=headers))

            print("High")
            print(tabulate(hTable, headers=headers))

    def inverse(
        self,
        tol,
        vel,
        tgtR,
        gunH=0,
        tgtH=0,
        env=atmosphere,
        t_max=1000,
        elev_min=-90,
        elev_max=90,
    ):
        """
        Inverse calculation: given shot splash range, calculate in inverse the
        angle necessary to achieve said range.
        """
        elev_min = max(elev_min, -90)
        elev_max = min(elev_max, 90)

        def f_r(elev, r=0):
            try:
                record = self.forward(tol, vel, elev, gunH, tgtH, env, t_max)
                t, (x, y, vx, vy) = record[-1]
            except ValueError:
                return -r, None

            theta = -(atan2(y, x) - 0.5 * pi)
            gr = theta * R_e
            # print(elev, r, gr - r)
            return gr - r, record[-1]

        """we assume 0 deg 0 min 1 sec is the maximum practical precision in
        terms of elevation angle that we care about.
        """
        elev_max = 0.5 * sum(
            gss(lambda ang: f_r(ang)[0], 0, 90, x_tol=3600**-1, findMin=False)
        )
        # print(elev_max)

        r_max = f_r(elev_max)[0]
        r_min = f_r(elev_min)[0]
        if isinstance(tgtR, int) or isinstance(tgtR, float):
            tgtR = [tgtR]

        lTrajs = []
        hTrajs = []
        for R in tgtR:
            if r_min < R < r_max:
                elev_i, elev_j = bisect(
                    lambda ang: f_r(ang, R)[0],
                    elev_min,
                    elev_max,
                    x_tol=3600**-1,
                )

                l_elev = 0.5 * (elev_i + elev_j)
                rec = f_r(l_elev)[1]
                if rec is not None:
                    lTrajs.append((l_elev, *rec))
                else:
                    lTrajs.append(None)

                elev_i, elev_j = bisect(
                    lambda ang: f_r(ang, R)[0], elev_max, 90, x_tol=3600**-1
                )
                h_elev = 0.5 * (elev_i + elev_j)
                print(R, h_elev)
                rec = f_r(h_elev)[1]
                if rec is not None:
                    hTrajs.append((h_elev, *rec))
                else:
                    hTrajs.append(None)

            else:
                lTrajs.append(None)
                hTrajs.append(None)

        return lTrajs, hTrajs

    def forward(
        self, tol, vel, elev, gunH=0, tgtH=0, env=atmosphere, t_max=1000
    ):
        """
        Forward calculation: given ballistic parameters, determine the
        trajector flown by the shot.
        """
        self.f_env = env

        """
         dv
        ---- = Kd(Mach) * rho  * V^2 / C
         dt

         C = M / (i D^2)
        """
        t_0 = 0
        if gunH == tgtH:
            gunH += tol

        x_0, y_0 = 0, R_e + gunH
        phi = elev * pi / 180

        def abortTgt(x, ys, o_x, o_ys):
            x, y, vx, vy = ys
            h = (x**2 + y**2) ** 0.5 - (R_e + tgtH)

            # only abort the calculation on downward crossing of target plane
            # height
            return h < 0 and ((-x * vx + y * vy) < 0)

        vx_0 = vel * cos(phi)
        vy_0 = vel * sin(phi)

        record = [[0, [x_0, y_0, vx_0, vy_0]]]

        t_2, vec_2, _ = RKF78(
            self._ode_t,
            (x_0, y_0, vx_0, vy_0),
            0,
            t_max,
            relTol=tol,
            abortFunc=abortTgt,
            adaptTo=(False, False, True, True),
            record=record,
        )  # Coarse integration to maximum time to find approximate ToF
        if t_2 == t_max:
            raise ValueError("Projectile Maximum Time-of-Flight t_max Exceeded")

        if len(record) > 1:
            t_1, (x_1, y_1, vx_1, vy_1) = record[-1]

        else:
            t_1, (x_1, y_1, vx_1, vy_1) = record[0]

        def f_tgt(t):
            _, (x, y, _, _), _ = RKF78(
                self._ode_t,
                (x_1, y_1, vx_1, vy_1),
                t_1,
                t,
                relTol=tol,
                adaptTo=(False, False, True, True),
            )  # fine integration from last point before impact
            return (x**2 + y**2) ** 0.5 - (R_e + tgtH)

        # print(t_1, x_1, y_1, vx_1, vy_1)
        # print(t_2, *vec_2)

        if (x_1**2 + y_1**2) ** 0.5 - (R_e + tgtH) < 0:
            """
            then we are probably cresting below the target plane.
            In this case we need to try raise the point t_1 to
            above the target plane *if possible*.

            This is a very edge case scenario but nevertheless
            in the name of accuracy and rigouroness it needs
            be done.
            """

            t_1_prime = 0.5 * sum(
                gss(f_tgt, t_1, t_2, x_tol=max(t_2, 1) * tol, findMin=False)
            )
            h_prime = f_tgt(t_1_prime)
            if h_prime > 0:
                # the new peak barely crest the target plane.
                t_t = 0.5 * sum(
                    bisect(f_tgt, t_1_prime, t_2, x_tol=max(t_2, 1) * tol)
                )

                # print(t_t, h_prime)

            else:
                # even the new peak point found cannot crest the target plane.
                raise ValueError(
                    "Projectile Cresting Below Target at {:.3f} m".format(
                        h_prime
                    )
                )

        else:
            t_t = 0.5 * sum(bisect(f_tgt, t_1, t_2, x_tol=max(t_2, 1) * tol))

        _, _, _ = RKF78(
            self._ode_t,
            (x_1, y_1, vx_1, vy_1),
            t_1,
            t_t,
            relTol=tol,
            adaptTo=(False, False, True, True),
            record=record,
        )

        return record


if __name__ == "__main__":
    test = Bullet(
        "test", mass=9.0990629, diam=88e-3, Kd_curve=KdCurve["G8"], form=0.925
    )
    # test.record_to_data(test.forward(tol=1e-6, vel=819.912, elev=45, tgtH=1))
    """
    print(
        *test.inverse(tol=1e-3, vel=819.92, tgtR=[1000, 2000, 8000], tgtH=10),
        sep="\n"
    )"""

    test.rangeTable(tol=1e-3, vel=819.2, minR=0, maxR=3000, deltaR=100, tgtH=10)

    """
    test = Bullet(
        "M2 ball",
        mass=0.045942,
        form=0.86,
        diam=12.7e-3,
        Kd_curve=KdCurve["G5"],
    )

    test.record_to_data(test.forward(tol=1e-9, vel=856 * 0.99, elev=0.725 * 2))
    """
