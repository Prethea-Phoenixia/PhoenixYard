from num import *
from gun import *


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

    def solve(self, loadFraction, chargeMassRatio, tol):
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
        """
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

        """
        step 1, find grain size that satisifies design pressure
        """
        p_bar_d = self.p_d / (self.f * Delta)  # convert to unitless

        def pressurePeaked(Z, t_bar, l_bar, v_bar):
            pass

        def _fp_e_1(e_1, tol):
            B = (
                self.S**2
                * e_1**2
                / (self.f * phi * omega * self.m * self.u_1**2)
                * (self.f * Delta) ** (2 * (1 - self.n))
            )

            # integrate this to end of burn

            def _ode_Z(Z, t_bar, l_bar, v_bar, p_bar):
                """burnout domain ode of internal ballistics"""
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

            def pressurePeaked(x, ys):
                Z = x
                t_bar, l_bar, v_bar, p_bar = ys

                _, _, _, dp = _ode_Z(Z, t_bar, l_bar, v_bar, p_bar)
                return dp < 0

            Z, (t_bar, l_bar, v_bar, p_bar), (_, _, _, _) = RKF78(
                dFunc=_ode_Z,
                iniVal=(0, 0, 0, p_bar_0),
                x_0=Z_0,
                x_1=self.Z_b,
                relTol=tol,
                absTol=tol,
                abortFunc=pressurePeaked,
            )

            """
            Expectantly, integrating backward may reduce the nbr.
            of steps required.
            """

            def _fp_Z(Z):
                _, (t_bar, l_bar, v_bar, p_bar), (_, _, _, _) = RKF78(
                    dFunc=_ode_Z,
                    iniVal=(0, 0, 0, p_bar_0),
                    x_0=0,
                    x_1=Z,
                    relTol=tol,
                    absTol=tol,
                )
                return p_bar

            Z_b_1, Z_b_2 = GSS(_fp_Z, Z_0, self.Z_b, relTol=tol, findMin=False)
            Z_b = 0.5 * (Z_b_1 + Z_b_2)
            return _fp_Z(Z_b) - p_bar_d, (t_bar, l_bar, v_bar, p_bar)

        """
        The two initial guesses are good enough for the majority of cases,
        guess one: 0.1mm, guess two: 1mm
        """

        print(_fp_e_1(0.1e-3, tol))

        print(_fp_e_1(1e-3, tol))

        e_1, _ = secant(
            lambda x: _fp_e_1(x, tol)[0],
            0.1e-3,
            1e-3,
            tol=p_bar_d * tol,
            x_min=1e-14,
        )  # this is the e_1 that satisifies the pressure specification.

        """
        step 2, find the requisite muzzle length to achieve design velocity
        """

        B = (
            self.S**2
            * e_1**2
            / (self.f * phi * omega * self.m * self.u_1**2)
            * (self.f * Delta) ** (2 * (1 - self.n))
        )

        def _ode_v(self, v_bar, t_bar, Z, l_bar):
            p_bar = _fp_bar(Z, l_bar, v_bar)
            if Z <= self.Z_b:
                dZ = (2 / (B * self.theta)) ** 0.5 * p_bar ** (self.n - 1)
            else:
                dZ = 0
            dl_bar = 2 * v_bar / (self.theta * p_bar)
            dt_bar = 2 / (self.theta * p_bar)

            return (dt_bar, dZ, dl_bar)


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

    print(test.solve(loadFraction=0.5, chargeMassRatio=0.5, tol=1e-3))
