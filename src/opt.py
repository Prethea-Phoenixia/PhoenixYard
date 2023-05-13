class Constrained:
    def __init__(
        self,
        caliber,
        shotMass,
        startPressure,
        chamberExpansion,
        dragCoe,
        designPressure,
        designVelocity,
    ):
        # constants for constrained designs
        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.propellant = propellant
        self.p_0 = startPressure
        self.chi_k = chamberExpansion
        self.phi_1 = 1 + dragCoe

        # design limits
        self.p_d = designPressure
        self.v_d = designVelocity

    def __getattr__(self, attrName):
        try:
            return getattr(self.propellant, attrName)
        except:
            raise AttributeError("object has no '%s'" % attrName)

    def f(loadFraction, chargeMassRatio):
        self.omega = self.m * chargeMassRatio
        self.Delta = self.omega / self.V_0
        self.phi = self.phi_1 + self.omega / (3 * self.m)
        """
        it is impossible to account for the chamberage effect given unspecified
        barrel length, in our formulation
        """
        self.v_j = (
            2 * self.f * self.omega / (self.theta * self.phi * self.m)
        ) ** 0.5

        self.psi_0 = (1 / self.Delta - 1 / self.rho_p) / (
            self.f / self.p_0 + self.alpha - 1 / self.rho_p
        )

        self.B = (
            self.S**2
            * e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_1**2)
            * (self.f * self.Delta) ** (2 * (1 - self.n))
        )

        Zs = cubic(
            self.chi * self.mu, self.chi * self.labda, self.chi, -self.psi_0
        )
        # pick a valid solution between 0 and 1
        Zs = tuple(
            Z
            for Z in Zs
            if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float

        self.Z_0 = Zs[0]

        def _fp_bar(self, Z, l_bar, v_bar):
            psi = self.f_psi_Z(Z)
            l_psi_bar = (
                1
                - self.Delta / self.rho_p
                - self.Delta * (self.alpha - 1 / self.rho_p) * psi
            )

            p_bar = (
                self.f * self.omega * psi
                - 0.5 * self.theta * self.phi * self.m * (v_bar * self.v_j) ** 2
            ) / (self.S * self.l_0 * (l_bar + l_psi_bar) * self.f * self.Delta)

            return p_bar

        def _ode_Z(self, Z, t_bar, l_bar, v_bar):
            """burnout domain ode of internal ballistics"""
            p_bar = self._fp_bar(Z, l_bar, v_bar)
            if Z <= self.Z_b:
                dt_bar = (2 * self.B / self.theta) ** 0.5 * p_bar**-self.n
                dl_bar = (
                    v_bar * (2 * self.B / self.theta) ** 0.5 * p_bar**-self.n
                )
                dv_bar = (self.B * self.theta * 0.5) ** 0.5 * p_bar ** (
                    1 - self.n
                )
            else:
                # technically speaking it is undefined in this area
                dt_bar = 0
                dl_bar = 0
                dv_bar = 0

            return (dt_bar, dl_bar, dv_bar)
