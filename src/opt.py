from num import *
from prop import *


class Constrained:
    def __init__(
        self,
        caliber,
        shotMass,
        propellant,
        startPressure,
        dragCoe,
        designPressure,
        designVelocity,
    ):
        # constants for constrained designs
        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.propellant = propellant
        self.p_0 = startPressure
        self.phi_1 = 1 + dragCoe

        # design limits
        self.p_d = designPressure
        self.v_d = designVelocity

    def __getattr__(self, attrName):
        try:
            return getattr(self.propellant, attrName)
        except:
            raise AttributeError("object has no '%s'" % attrName)

    def solve(
        self, loadFraction, chargeMassRatio, tol, minWeb=1e-6, webTol=1e-6
    ):
        """
        minWeb  : represents minimum possible grain size
        """
        omega = self.m * chargeMassRatio
        V_0 = omega / (self.rho_p * self.maxLF * loadFraction)
        Delta = omega / V_0
        p_bar_0 = self.p_0 / (Delta * self.f)
        l_0 = V_0 / self.S
        phi = self.phi_1 + omega / (3 * self.m)
        """
        it is impossible to account for the chamberage effect given unspecified
        barrel length, in our formulation
        """
        v_j = (2 * self.f * omega / (self.theta * phi * self.m)) ** 0.5

        if v_j < 0.75 * self.v_d:
            raise ValueError(
                "Propellant load too low to achieve " "design velocity."
            )

        psi_0 = (1 / Delta - 1 / self.rho_p) / (
            self.f / self.p_0 + self.alpha - 1 / self.rho_p
        )

        Zs = cubic(
            a=self.chi * self.mu,
            b=self.chi * self.labda,
            c=self.chi,
            d=-psi_0,
        )
        # pick a valid solution between 0 and 1
        Zs = tuple(
            Z
            for Z in Zs
            if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float

        Z_0 = Zs[0]

        def _fp_bar(Z, l_bar, v_bar):
            psi = self.f_psi_Z(Z)
            l_psi_bar = (
                1
                - Delta / self.rho_p
                - Delta * (self.alpha - 1 / self.rho_p) * psi
            )

            p_bar = (
                self.f * omega * psi
                - 0.5 * self.theta * phi * self.m * (v_bar * v_j) ** 2
            ) / (self.S * l_0 * (l_bar + l_psi_bar) * self.f * Delta)

            return p_bar

        """
        step 1, find grain size that satisifies design pressure
        """
        p_bar_d = self.p_d / (self.f * Delta)  # convert to unitless

        def _fp_e_1(e_1, tol):
            """
            calculate either the peak pressure, given the arc thickness,
            or until the system develops 2x design pressure. The latter
            is designed to provide a flat
            """

            # e_1 = round(e_1 / webTol) * webTol

            B = (
                self.S**2
                * e_1**2
                / (self.f * phi * omega * self.m * self.u_1**2)
                * (self.f * Delta) ** (2 * (1 - self.n))
            )

            # integrate this to end of burn

            def _pf_Z(x, *ys):
                Z, t_bar, l_bar, v_bar, p_bar = x, *ys
                p_bar_prime = _fp_bar(Z, l_bar, v_bar)
                return (*ys[:-1], p_bar_prime)

            def _ode_Z(Z, t_bar, l_bar, v_bar, p_bar):
                """burnout domain ode of internal ballistics"""
                psi = self.f_psi_Z(Z)
                dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
                l_psi_bar = (
                    1
                    - Delta / self.rho_p
                    - Delta * (self.alpha - 1 / self.rho_p) * psi
                )
                # dp_bar/dt_bar

                if Z <= self.Z_b:
                    dt_bar = (
                        2 * B / self.theta
                    ) ** 0.5 * p_bar**-self.n  # dt_bar/dZ

                    dl_bar = (
                        v_bar * (2 * B / self.theta) ** 0.5 * p_bar**-self.n
                    )  # dl_bar/dZ

                    dv_bar = (B * self.theta * 0.5) ** 0.5 * p_bar ** (
                        1 - self.n
                    )

                    dp_bar = (
                        (
                            (1 + p_bar * Delta * (self.alpha - 1 / self.rho_p))
                            * dpsi
                            / dt_bar
                            - p_bar * v_bar * (1 + self.theta)
                        )
                        * dt_bar
                        / (l_bar + l_psi_bar)
                    )

                else:
                    dt_bar = 0
                    dl_bar = 0
                    dv_bar = 0
                    dp_bar = 0

                return (dt_bar, dl_bar, dv_bar, dp_bar)

            def abort(x, ys, dys):
                Z = x
                t_bar, l_bar, v_bar, p_bar = ys
                dt_bar, dl_bar, dv_abr, dp_bar = dys

                return dp_bar < 0 or p_bar > 2 * p_bar_d

            record = []

            Z_j, (t_bar_j, l_bar_j, v_bar_j, p_bar_j), (_, _, _, _) = RKF78(
                dFunc=_ode_Z,
                iniVal=(0, 0, 0, p_bar_0),
                x_0=Z_0,
                x_1=self.Z_b,
                relTol=tol,
                absTol=tol,
                abortFunc=abort,
                parFunc=_pf_Z,
                record=record,
            )

            Z_i, (t_bar_i, l_bar_i, v_bar_i, p_bar_i) = record[-2]

            def _fp_Z(Z):
                _, (t_bar, l_bar, v_bar, p_bar), (_, _, _, _) = RKF78(
                    dFunc=_ode_Z,
                    iniVal=(t_bar_i, l_bar_i, v_bar_i, p_bar_i),
                    x_0=Z_i,
                    x_1=Z,
                    relTol=tol,
                    absTol=tol,
                )
                return p_bar

            Z_1, Z_2 = GSS(
                _fp_Z,
                Z_i,
                Z_j,
                relTol=tol,
                findMin=False,
                absTol=1e-16,  # limitation of double.
            )  # find the peak pressure point.

            Z_p = 0.5 * (Z_1 + Z_2)

            return (
                _fp_Z(Z_p) - p_bar_d,
                Z_j,
                (t_bar_j, l_bar_j, v_bar_j, p_bar_j),
            )

        """
        The two initial guesses are good enough for the majority of cases,
        guess one: 0.1mm, guess two: 1mm
        """
        """
        if _fp_e_1(minWeb, tol)[0] < 0:
            raise ValueError(
                "Design pressure cannot be achieved within minimum web"
            )
        """
        dp_bar_probe = _fp_e_1(minWeb, tol)[0]
        probeWeb = minWeb

        if dp_bar_probe < 0:
            raise ValueError(
                "Design pressure cannot be achieved by varying web down to minimum"
            )

        while dp_bar_probe > 0:
            probeWeb *= 2
            dp_bar_probe = _fp_e_1(probeWeb, tol)[0]

        # print("x_min=", 0.1 * probeWeb)

        e_1, _ = secant(
            lambda x: _fp_e_1(x, tol)[0],
            probeWeb,  # >0
            0.5 * probeWeb,  # ?0
            tol=p_bar_d * tol,
            x_min=0.5 * probeWeb,  # <=0
        )  # this is the e_1 that satisifies the pressure specification.

        """
        e_1, e_1_prime = bisect(
            lambda x: _fp_e_1(x, tol)[0],
            probeWeb,  # >0
            0.1 * probeWeb,  # <=0
            tol=webTol,
        )
        """

        p_bar_dev, Z_i, (t_bar_i, l_bar_i, v_bar_i, p_bar_i) = _fp_e_1(e_1, tol)

        if abs(p_bar_dev / p_bar_d) > tol:
            raise ValueError("Design pressure is not met")

        """
        step 2, find the requisite muzzle length to achieve design velocity
        """
        v_bar_d = self.v_d / v_j

        if v_bar_i > v_bar_d:
            return ValueError("Design velocity exceeded before peak pressure")
        else:
            pass

        B = (
            self.S**2
            * e_1**2
            / (self.f * phi * omega * self.m * self.u_1**2)
            * (self.f * Delta) ** (2 * (1 - self.n))
        )

        def _pf_v(x, *ys):
            v_bar, t_bar, Z, l_bar, p_bar = x, *ys
            p_bar_prime = _fp_bar(Z, l_bar, v_bar)
            return (*ys[:-1], p_bar_prime)

        def _ode_v(v_bar, t_bar, Z, l_bar, p_bar):
            # p_bar = _fp_bar(Z, l_bar, v_bar)
            psi = self.f_psi_Z(Z)
            dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
            l_psi_bar = (
                1
                - Delta / self.rho_p
                - Delta * (self.alpha - 1 / self.rho_p) * psi
            )
            # dp_bar/dt_bar

            if Z <= self.Z_b:
                dZ = (2 / (B * self.theta)) ** 0.5 * p_bar ** (self.n - 1)
            else:
                dZ = 0
            dl_bar = 2 * v_bar / (self.theta * p_bar)
            dt_bar = 2 / (self.theta * p_bar)
            dp_bar = (
                (
                    (1 + p_bar * Delta * (self.alpha - 1 / self.rho_p))
                    * dpsi
                    * dZ
                    / dt_bar
                    - p_bar * v_bar * (1 + self.theta)
                )
                * dt_bar
                / (l_bar + l_psi_bar)
            )

            return (dt_bar, dZ, dl_bar, dp_bar)

        v_bar_g, (t_bar_g, Z_g, l_bar_g, p_bar_g), (_, _, _, _) = RKF78(
            dFunc=_ode_v,
            iniVal=(t_bar_i, Z_i, l_bar_i, p_bar_i),
            x_0=v_bar_i,
            x_1=v_bar_d,
            relTol=tol,
            absTol=None,
            parFunc=_pf_v,
        )

        if abs(v_bar_g - v_bar_d) / (v_bar_d) > tol:
            raise ValueError("Velocity specification is not met")

        return e_1, l_bar_g * l_0


if __name__ == "__main__":
    compositions = GrainComp.readFile("data/propellants.csv")
    JA2 = compositions["JA2"]
    JA2SHC = Propellant(JA2, MultPerfGeometry.SEVEN_PERF_CYLINDER, 2, 2.5)

    test = Constrained(
        caliber=50e-3,
        shotMass=1,
        propellant=JA2SHC,
        startPressure=10e6,
        dragCoe=5e-2,
        designPressure=300e6,
        designVelocity=1200,
    )
    """
    print(
        test.solve(
            loadFraction=0.3,
            chargeMassRatio=0.8,
            tol=1e-3,
        )
    )
    """

    for i in range(9):
        loadFraction = (i + 1) / 10
        for j in range(10):
            chargeMassRatio = 0.2 + j / 5
            print("lf=", loadFraction, "c.m.=", chargeMassRatio)
            try:
                print(
                    test.solve(
                        loadFraction=loadFraction,
                        chargeMassRatio=chargeMassRatio,
                        tol=1e-3,
                    )
                )
            except ValueError as e:
                print(e)
