from atmos import atmosphere, R_e
from drag import KdCurve
from math import pi, sin, cos, atan2, acos

from num import RKF78, secant, GSS, gss


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
        # convert the record into more sensible values:

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
                    headers=(
                        "ToF",
                        "H",
                        "G.R.",
                        "V.",
                        "Angle",
                    ),
                )
            )

        return data

    def inverse(
        self, tol, vel, tgtR, gunH=0, tgtH=0, env=atmosphere, t_max=1000
    ):
        """
        Inverse calculation: given shot splash range, calculate in inverse the
        angle necessary to achieve said range.
        """

        def f_r(elev):
            record = self.forward(tol, vel, elev, gunH, tgtH, env, t_max)
            t, (x, y, vx, vy) = record[-1]

            theta = -(atan2(y, x) - 0.5 * pi)
            gr = theta * R_e

            return gr

        # print(GSS(f_r, 0, 90, yRelTol=tol, findMin=False))
        print(gss(f_r, 0, 90, tol * 90, findMin=False))

    def forward(
        self, tol, vel, elev, gunH=0, tgtH=0, env=atmosphere, t_max=1000
    ):
        """
        Forward calculation: given ballistic parameters, determine the
        trajector flown by the shot.
        Approximate Error in relation to tol:

         tol  |   ε
        1e-6  |  6.5 m
        1e-7  |  6.5 hm
        1e-8  |  6.5 dm
        1e-9  |  6.5 cm
        1e-10 |  6.5 mm

        """
        self.f_env = env

        """
         dv
        ---- = Kd(Mach) * rho  * V^2 / C
         dt

         C = M / (i D^2)
        """
        x_0, y_0 = 0, R_e + gunH
        phi = elev * pi / 180

        def abortTgt(x, ys, o_x, o_ys):
            x, y, vx, vy = ys
            h = (x**2 + y**2) ** 0.5 - (R_e + tgtH)
            # v = (vx**2 + vy**2) ** 0.5
            return h < 0

        vx_0 = vel * cos(phi)
        vy_0 = vel * sin(phi)

        record = []
        t_1, _, _ = RKF78(
            self._ode_t,
            (x_0, y_0, vx_0, vy_0),
            0,
            t_max,
            relTol=tol,
            minTol=1e-14,
            abortFunc=abortTgt,
            record=record,
        )  # Coarse integration to maximum time to find approximate ToF

        t_0, (x_0, y_0, vx_0, vy_0) = record[-2]

        if t_1 == t_max:
            raise ValueError("Projectile Maximum Time-of-Flight t_max Exceeded")

        def f_tgt(t):
            _, (x, y, _, _), _ = RKF78(
                self._ode_t,
                (x_0, y_0, vx_0, vy_0),
                t_0,
                t,
                relTol=tol,
                minTol=1e-14,
            )  # fine integration from last point before impact
            return (x**2 + y**2) ** 0.5 - (R_e + tgtH)

        t_t, _ = secant(
            f_tgt, 0, t_0, t_1, x_tol=t_0 * tol
        )  # time to target is determined via secant search

        _, (x_1, y_1, vx_1, vy_1), _ = RKF78(
            self._ode_t,
            (x_0, y_0, vx_0, vy_0),
            t_0,
            t_t,
            relTol=tol,
            minTol=1e-14,
            record=record,
        )

        return record


if __name__ == "__main__":
    test = Bullet(
        "test", mass=9.0990629, diam=88e-3, Kd_curve=KdCurve["G8"], form=0.925
    )

    test.record_to_data(test.forward(tol=1e-9, vel=819.912, elev=25))

    # test.inverse(tol=1e-9, vel=819.92, tgtR=8000)
    """
    test = Bullet(
        "M2 ball",
        mass=0.045942,
        form=0.86,
        diam=12.7e-3,
        Kd_curve=KdCurve["G5"],
    )

    test.record_to_data(test.forward(tol=1e-9, vel=856, elev=0.725))
    """
