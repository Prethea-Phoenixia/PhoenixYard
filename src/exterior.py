from atmos import atmosphere, R_e
from drag import KdCurve
from math import pi, sin, cos, atan2, acos

from num import RKF78, gss, bisect


def dec_to_dms(deg):
    if deg >= 0:
        sign = 1
    else:
        sign = -1

    deg = abs(deg)

    h = int(deg)
    m = int((deg - h) * 60)
    s = int((deg - h - m / 60) * 3600)

    return h * sign, m, s


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

        print(Rs)

        lTrajs, hTrajs = self.inverse(tol, vel, Rs, gunH, tgtH, env, t_max)

        lTable, hTable = [], []
        for R, lTraj, hTraj in zip(Rs, lTrajs, hTrajs):
            if (lTraj is not None) and (hTraj is not None):
                elev, t, (x, y, vx, vy) = lTraj
                lTable.append(
                    (
                        R,
                        "{:}°{:}'{:}\"".format(*dec_to_dms(elev)),
                        t,
                    )
                )
                elev, t, (x, y, vx, vy) = hTraj
                hTable.append(
                    (
                        R,
                        "{:}°{:}'{:}\"".format(*dec_to_dms(elev)),
                        t,
                    )
                )

        if prettyprint:
            from tabulate import tabulate

            print("Low")
            print(tabulate(lTable))

            print("High")
            print(tabulate(hTable))

    def inverse(
        self, tol, vel, tgtR, gunH=0, tgtH=0, env=atmosphere, t_max=1000
    ):
        """
        Inverse calculation: given shot splash range, calculate in inverse the
        angle necessary to achieve said range.
        """

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
        print(elev_max)

        r_max = f_r(elev_max)[0]
        r_min = f_r(0)[0]

        if isinstance(tgtR, int) or isinstance(tgtR, float):
            tgtR = [tgtR]

        lTrajs = []
        hTrajs = []
        for R in tgtR:
            if r_min < R < r_max:
                l_elev = 0.5 * sum(
                    bisect(
                        lambda ang: f_r(ang, R)[0],
                        0,
                        elev_max,
                        x_tol=3600**-1,
                    )
                )
                h_elev = 0.5 * sum(
                    bisect(
                        lambda ang: f_r(ang, R)[0],
                        elev_max,
                        90,
                        x_tol=3600**-1,
                    )
                )

                lTrajs.append((l_elev, *f_r(l_elev)[1]))
                hTrajs.append((h_elev, *f_r(h_elev)[1]))

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

            return h < 0

        vx_0 = vel * cos(phi)
        vy_0 = vel * sin(phi)

        record = [[0, [x_0, y_0, vx_0, vy_0]]]

        t_2, _, _ = RKF78(
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

    test.record_to_data(
        test.forward(tol=1e-3, vel=819.912, elev=4.256, tgtH=100)
    )

    print(
        *test.inverse(tol=1e-3, vel=819.92, tgtR=[1000, 2000, 8000], tgtH=100),
        sep="\n"
    )

    """
    test.rangeTable(
        tol=1e-3, vel=819.2, minR=0, maxR=5000, deltaR=1000
    )
    """
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
