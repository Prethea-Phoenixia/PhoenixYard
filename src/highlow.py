from math import pi
from num import gss, RKF78, cubic, bisect
from prop import GrainComp, Propellant

from gun import DOMAIN_TIME, DOMAIN_LENG
from gun import (
    POINT_START,
    POINT_PEAK,
    POINT_PEAK_SHOT,
    POINT_PEAK_BREECH,
    POINT_FRACTURE,
    POINT_BURNOUT,
    POINT_EXIT,
)


class Highlow:
    def __init__(
        self,
        caliber,
        shotMass,
        propellant,
        grainSize,
        chargeMass,
        chamberVolume,  # high pressure chamber
        expansionVolume,  # low pressure chamber
        expansionStartPressure,  # low pressure chamber starting pressure
        startPressure,  # shot start pressure
        portAreaRatio,  # total area of passable "mesh" relative to bore.
        chambrage,
        lengthGun,
        nozzleExpansion,
        dragCoefficient=0,
        **_,
    ):
        if any(
            (
                caliber <= 0,
                shotMass <= 0,
                chargeMass <= 0,
                grainSize <= 0,
                chamberVolume <= 0,
                expansionVolume <= 0,
                lengthGun <= 0,
                nozzleExpansion < 1,
                dragCoefficient < 0,
                dragCoefficient >= 1,
                expansionStartPressure < 0,
                startPressure < 0,
            )
        ):
            raise ValueError("Invalid gun parameters")

        if chargeMass > (propellant.maxLF * propellant.rho_p * chamberVolume):
            raise ValueError(
                "Specified Load Fraction Violates Geometrical Constraint"
            )

        self.propellant = propellant

        e_1 = 0.5 * grainSize
        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.omega = chargeMass

        self.V_0 = chamberVolume
        self.V_1 = expansionVolume

        self.p_0_e = expansionStartPressure
        self.p_0_s = startPressure

        self.l_g = lengthGun
        self.A_bar = nozzleExpansion
        self.chi_k = chambrage

        self.Delta = self.omega / self.V_0
        self.l_0 = self.V_0 / self.S

        self.phi_1 = 1 / (1 - dragCoefficient)  # drag work coefficient
        self.phi = self.phi_1 + self.omega / (3 * self.m)

        self.v_j = (
            2 * self.f * self.omega / (self.theta * self.phi * self.m)
        ) ** 0.5

        if self.p_0_e == 0:
            raise NotImplementedError(
                "Current implementation use exponential burn rate and does not"
                + " allow for solving the case with 0 shot start pressure."
            )
        else:
            self.psi_0 = (1 / self.Delta - 1 / self.rho_p) / (
                self.f / self.p_0_e + self.alpha - 1 / self.rho_p
            )
            if self.psi_0 <= 0:
                raise ValueError(
                    "Initial burnup fraction is solved to be negative."
                    + " In practice this implies a detonation of the gun breech"
                    + " will likely occur."
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
        Zs = sorted(
            Z
            for Z in Zs
            if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to overcome"
                + " start pressure, or has burnt to post fracture."
            )

        self.Z_0 = Zs[0]  # pick the smallest solution
        self.S_j_bar = portAreaRatio
        self.S_j = self.S * portAreaRatio

        gamma = self.theta + 1
        phi_2 = 0.15  # for small ports 1.5mm to 2mm in size
        self.C_A = (
            phi_2
            * (0.5 * self.theta * self.phi * self.m / self.omega) ** 0.5
            * gamma**0.5
            * (2 / (gamma + 1)) ** (0.5 * (gamma + 1) / self.theta)
        )
        self.C_B = (
            phi_2
            * (0.5 * self.theta * self.phi * self.m / self.omega) ** 0.5
            * ((2 * gamma) / self.theta) ** 0.5
        )

        self.cfpr = (2 / (gamma + 1)) ** (
            gamma / self.theta
        )  # pressure ratio threshold for critical flow
        self.V_bar = self.V_1 / self.V_0

    def __getattr__(self, attrName):
        if "propellant" in vars(self) and not (
            attrName.startswith("__") and attrName.endswith("__")
        ):
            try:
                return getattr(self.propellant, attrName)
            except AttributeError:
                AttributeError(
                    "%r object has no attribute %r"
                    % (self.__class__.__name__, attrName)
                )
        else:
            raise AttributeError

    def _f_p_1_bar(self, Z, eta, tau_1, psi=None):
        psi = psi if psi else self.f_psi_Z(Z)
        l_psi_bar = 1 - self.Delta * (
            (1 - psi) / self.rho_p + self.alpha * (psi - eta)
        )

        return tau_1 / l_psi_bar * (psi - eta)

    def _f_p_2_bar(self, l_bar, eta, tau_2):
        return (
            tau_2 * eta / (l_bar + self.V_bar - self.alpha * self.Delta * eta)
        )

    def _ode_t(self, t_bar, Z, l_bar, v_bar, eta, tau_1, tau_2):
        # VALUES FOR HIGH PRESSURE CHAMBER
        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
        p_1_bar = self._f_p_1_bar(Z, eta, tau_1, psi)

        # VALUES FOR LOW PRESSURE CHAMBER
        p_2_bar = self._f_p_2_bar(l_bar, eta, tau_2)

        # COMBUSTION OF POWDER
        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_1_bar**self.n
        else:
            dZ = 0  # dZ/dt_bar

        # FLOW INTO LOW PRESSURE CHAMBER
        pr = p_2_bar / p_1_bar
        if pr < self.cfpr:
            deta = self.C_A * self.S_j_bar * p_1_bar  # deta / dt_bar
        else:
            gamma = self.theta + 1
            deta = (
                self.C_B
                * (pr ** (2 / gamma) - pr ** ((gamma + 1) / gamma)) ** 0.5
                * self.S_j_bar
                * p_1_bar
            )

        # CHANGE IN HIGH PRESSURE CHAMBER
        dtau_1 = ((1 - tau_1) * (dpsi * dZ) - self.theta * tau_1 * deta) / (
            psi - eta
        )  # dtau/dt_bar

        # ATMOSPHERIC DRAG
        if self.c_a_bar != 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                +0.25 * k * (k + 1) * v_r**2
                + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar
        else:
            p_d_bar = 0

        # PROJECTILE MOTION
        dv_bar = self.theta * 0.5 * (p_2_bar - p_d_bar)
        dl_bar = v_bar

        # CHANGE IN LOW PRESSURE CHAMBER

        if eta == 0:
            dtau_2 = 0
        else:
            dtau_2 = (self.theta * tau_1 * deta - 2 * v_bar * dv_bar) / eta

        return (dZ, dl_bar, dv_bar, deta, dtau_1, dtau_2)

    def _ode_ts(self, t_bar, Z, eta, tau_1, tau_2):
        # VALUES FOR HIGH PRESSURE CHAMBER
        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
        p_1_bar = self._f_p_1_bar(Z, eta, tau_1, psi)
        # VALUES FOR LOW PRESSURE CHAMBER
        p_2_bar = self._f_p_2_bar(0, eta, tau_2)

        # COMBUSTION OF POWDER
        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_1_bar**self.n
        else:
            dZ = 0  # dZ/dt_bar

        # FLOW INTO LOW PRESSURE CHAMBER
        pr = p_2_bar / p_1_bar

        if pr < self.cfpr:
            # FLOW IS CRITICAL
            deta = self.C_A * self.S_j_bar * p_1_bar  # deta / dt_bar

        else:
            # FLOW IS SUB CRITICAL
            gamma = self.theta + 1
            deta = (
                self.C_B
                * (pr ** (2 / gamma) - pr ** ((gamma + 1) / gamma)) ** 0.5
                * self.S_j_bar
                * p_1_bar
            )

        # CHANGE IN HIGH PRESSURE CHAMBER
        dtau_1 = ((1 - tau_1) * (dpsi * dZ) - self.theta * tau_1 * deta) / (
            psi - eta
        )  # dtau/dt_bar

        # CHANGE IN LOW PRESSURE CHAMBER

        if eta == 0:
            dtau_2 = 0
        else:
            dtau_2 = (self.theta * tau_1 * deta) / eta

        return (dZ, deta, dtau_1, dtau_2)

    def _ode_l(self, l_bar, t_bar, Z, v_bar, eta, tau_1, tau_2):
        # VALUES FOR HIGH PRESSURE CHAMBER
        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
        p_1_bar = self._f_p_1_bar(Z, eta, tau_1, psi)

        # VALUES FOR LOW PRESSURE CHAMBER
        p_2_bar = self._f_p_2_bar(l_bar, eta, tau_2)

        # COMBUSTION OF POWDER
        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_1_bar**self.n / v_bar
        else:
            dZ = 0  # dZ /dl_bar

        # FLOW INTO LOW PRESSURE CHAMBER
        pr = p_2_bar / p_1_bar
        if pr < self.cfpr:
            # FLOW IS CRITICAL
            deta = self.C_A * self.S_j_bar * p_1_bar / v_bar  # deta / dl_bar
        else:
            # FLOW IS SUBCRITICAL
            gamma = self.theta + 1
            deta = (
                self.C_B
                * (pr ** (2 / gamma) - pr ** ((gamma + 1) / gamma)) ** 0.5
                * self.S_j_bar
                * p_1_bar
            ) / v_bar

        # CHANGE IN HIGH PRESSURE CHAMBER
        dtau_1 = ((1 - tau_1) * (dpsi * dZ) - self.theta * tau_1 * (deta)) / (
            psi - eta
        )  # dtau/dl_bar

        # ATMOSPHERIC DRAG
        if self.c_a_bar != 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                +0.25 * k * (k + 1) * v_r**2
                + k * v_r * (1 + (0.25 * (k + 1) * v_r) ** 2) ** 0.5
            ) * self.p_a_bar
        else:
            p_d_bar = 0

        # PROJECTILE MOTION
        dt_bar = 1 / v_bar  # dt_bar / dl_bar
        dv_bar = (
            self.theta * 0.5 * (p_2_bar - p_d_bar) * dt_bar
        )  # dv_bar/dl_bar

        # CHANGE IN LOW PRESSURE CHAMBER
        if eta == 0:
            dtau_2 = 0
        else:
            dtau_2 = (
                self.theta * tau_1 * deta - 2 * v_bar * dv_bar
            ) / eta  # dtau_2/dl_bar

        return (dt_bar, dZ, dv_bar, deta, dtau_1, dtau_2)

    def _ode_Z(self, Z, t_bar, l_bar, v_bar, eta, tau_1, tau_2):
        # VALUES FOR HIGH PRESSURE CHAMBER
        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
        p_1_bar = self._f_p_1_bar(Z, eta, tau_1, psi)
        # VALUES FOR LOW PRESSURE CHAMBER
        p_2_bar = self._f_p_2_bar(l_bar, eta, tau_2)

        # PASSAGE OF TIME
        dt_bar = (2 * self.B / self.theta) ** 0.5 * p_1_bar**-self.n

        # FLOW INTO LOW PRESSURE CHAMBER
        pr = p_2_bar / p_1_bar
        if pr < self.cfpr:
            # FLOW IS CRITICAL
            deta = self.C_A * self.S_j_bar * p_1_bar * dt_bar  # deta / dZ
        else:
            # FLOW IS SUB CRITICAL
            gamma = self.theta + 1
            deta = (
                self.C_B
                * (pr ** (2 / gamma) - pr ** ((gamma + 1) / gamma)) ** 0.5
                * self.S_j_bar
                * p_1_bar
            ) * dt_bar

        # CHANGE IN HIGH PRESSURE CHAMBER

        dtau_1 = ((1 - tau_1) * (dpsi) - self.theta * tau_1 * (deta)) / (
            psi - eta
        )  # dtau_1 / dZ

        # ATMOSPHERIC DRAG
        if self.c_a_bar != 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                +0.25 * k * (k + 1) * v_r**2
                + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar
        else:
            p_d_bar = 0

        # PROJECTILE MOTION

        dl_bar = v_bar * dt_bar
        dv_bar = 0.5 * self.theta * (p_2_bar - p_d_bar) * dt_bar  # dv_bar / dZ

        # CHANGE IN LOW PRESSURE CHAMBER

        if eta == 0:
            dtau_2 = 0
        else:
            dtau_2 = (
                self.theta * tau_1 * deta - 2 * v_bar * dv_bar
            ) / eta  # dtau_2/dZ

        return (dt_bar, dl_bar, dv_bar, deta, dtau_1, dtau_2)

    def _ode_Zs(self, Z, t_bar, eta, tau_1, tau_2):
        # VALUES FOR HIGH PRESSURE CHAMBER
        # print("in", Z, t_bar, eta, tau_1, tau_2)
        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
        p_1_bar = self._f_p_1_bar(Z, eta, tau_1, psi)
        # VALUES FOR LOW PRESSURE CHAMBER
        p_2_bar = self._f_p_2_bar(0, eta, tau_2)

        # PASSAGE OF TIME
        dt_bar = (2 * self.B / self.theta) ** 0.5 * p_1_bar**-self.n

        # FLOW INTO LOW PRESSURE CHAMBER
        pr = p_2_bar / p_1_bar
        if pr < self.cfpr:
            # FLOW IS CRITICAL
            deta = self.C_A * self.S_j_bar * p_1_bar * dt_bar  # deta / dZ
        else:
            # FLOW IS SUB CRITICAL
            gamma = self.theta + 1
            deta = (
                self.C_B
                * (pr ** (2 / gamma) - pr ** ((gamma + 1) / gamma)) ** 0.5
                * self.S_j_bar
                * p_1_bar
            ) * dt_bar

        dtau_1 = ((1 - tau_1) * (dpsi) - self.theta * tau_1 * (deta)) / (
            psi - eta
        )  # dtau_1 / dZ

        # CHANGE IN LOW PRESSURE CHAMBER
        if eta == 0:
            dtau_2 = 0
        else:
            dtau_2 = (self.theta * tau_1 * deta) / eta  # dtau_2/dZ

        return (dt_bar, deta, dtau_1, dtau_2)

    def integrate(
        self,
        step=10,
        tol=1e-5,
        dom=DOMAIN_TIME,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        **_,
    ):
        """
        Runs a full numerical solution for the gun in the specified domain sampled
        evenly at specified number of steps, using a scaled numerical tolerance as
        specified.

        tolerance is meant to be interpreted as the maximum relative deviation each
        component is allowed to have, at each step of integration point.

        Through significant trials and errors, it was determined that for this particular
        system, the error due to compounding does not appear to be significant,
        usually on the order of 1e-16 - 1e-14 as compared to much larger for component
        errors.
        """
        record = []
        minTol = 1e-16  # based on experience

        if any((step < 0, tol < 0)):
            raise ValueError("Invalid integration specification")

        if any((ambientP < 0, ambientRho < 0, ambientGamma < 1)):
            raise ValueError("Invalid ambient condition")

        tScale = self.l_0 / self.v_j
        pScale = self.f * self.Delta

        self.p_0_e_bar = self.p_0_e / pScale
        self.p_0_s_bar = self.p_0_s / pScale

        # ambient conditions
        self.p_a_bar = ambientP / pScale

        if ambientRho != 0:
            self.c_a_bar = (
                ambientGamma * ambientP / ambientRho
            ) ** 0.5 / self.v_j
        else:
            self.c_a_bar = 0

        self.k_1 = ambientGamma

        l_g_bar = self.l_g / self.l_0
        p_bar_0 = self.p_0_e / pScale
        Z_b = self.Z_b
        Z_0 = self.Z_0

        bar_data = []
        bar_err = []

        # fmt: off
        def updBarData(
            tag, t_bar, l_bar, Z, v_bar, eta, tau_1, tau_2, t_bar_err, l_bar_err,
            Z_err, v_bar_err, eta_err, tau_1_err, tau_2_err,
        ):
            p_1_bar, p_2_bar = (
                self._f_p_1_bar(Z, eta, tau_1),
                self._f_p_2_bar(l_bar, eta, tau_2),
            )
            p_bar_err = "N/A"
            bar_data.append(
                (
                    tag, t_bar, l_bar, Z, v_bar, p_1_bar, p_2_bar, eta, tau_1,
                    tau_2,
                )
            )
            bar_err.append(
                (
                    "L", t_bar_err, l_bar_err, Z_err, v_bar_err, p_bar_err,
                    p_bar_err, eta_err, tau_1_err, tau_2_err
                )
            )
        # fmt: on
        updBarData(
            tag=None,
            t_bar=0,
            l_bar=0,
            Z=Z_0,
            v_bar=0,
            eta=0,
            tau_1=1,
            tau_2=0,
            t_bar_err=0,
            l_bar_err=0,
            Z_err=0,
            v_bar_err=0,
            eta_err=0,
            tau_1_err=0,
            tau_2_err=0,
        )

        def abort_s(x, ys, o_x, o_ys):
            Z, t_bar, eta, tau_1, tau_2 = x, *ys
            p_2_bar = self._f_p_2_bar(0, eta, tau_2)

            return p_2_bar > 2 * self.p_0_s_bar

        def f(Z):
            t_bar, eta, tau_1, tau_2 = RKF78(
                self._ode_Zs,
                (0, tol, 1, 0),
                Z_0,
                Z,
                relTol=tol,
                absTol=tol**2,
                minTol=minTol,
                debug=False,
                abortFunc=abort_s,
            )[1]

            p_2_bar = self._f_p_2_bar(0, eta, tau_2)

            return p_2_bar - self.p_0_s_bar

        Z_1 = 0.5 * sum(bisect(f, Z_0, Z_b, y_abs_tol=self.p_0_s_bar * tol))

        t_bar_1, eta_1, tau_1_1, tau_2_1 = RKF78(
            self._ode_Zs,
            (0, tol, 1, 0),
            Z_0,
            Z_1,
            relTol=tol,
            absTol=tol**2,
            minTol=minTol,
            debug=False,
            abortFunc=abort_s,
        )[1]

        updBarData(
            tag=POINT_START,
            t_bar=t_bar_1,
            l_bar=0,
            Z=Z_1,
            v_bar=0,
            eta=eta_1,
            tau_1=tau_1_1,
            tau_2=tau_2_1,
            t_bar_err=0,
            l_bar_err=0,
            Z_err=0,
            v_bar_err=0,
            eta_err=0,
            tau_1_err=0,
            tau_2_err=0,
        )

        Z_i = Z_1
        Z_j = Z_b
        N = 1
        Delta_Z = Z_b - Z_0

        # fmt: off
        t_bar_i, l_bar_i, v_bar_i, p_bar_i, eta_i, tau_1_i, tau_2_i = (
            t_bar_1, 0, 0, self.p_0_s_bar, eta_1, tau_1_1, tau_2_1
        )
        # fmt: on

        isBurnOutContained = True

        """
        Instead of letting the integrator handle the heavy lifting, we
        partition Z and integrate upwards until either barrel exit or
        burnout has been achieved. This seek-and-validate inspired
        from bisection puts the group of value denoted by subscript
        i within or on the muzzle, with the propellant either still
        burning or right on the burnout point..
        """
        ztlvett_record = [(Z_1, (t_bar_1, 0, 0, eta_1, tau_1_1, tau_2_1))]
        p_max = 1e9  # 1GPa
        p_bar_max = p_max / pScale

        def abort(x, ys, o_x, o_ys):
            Z, t_bar, l_bar, v_bar, eta, tau_1, tau_2 = x, *ys

            p_1_bar = self._f_p_1_bar(Z, eta, tau_1)

            return any(
                (l_bar > l_g_bar, p_1_bar > p_bar_max, v_bar <= 0, p_1_bar < 0)
            )

        while Z_i < Z_b:  # terminates if burnout is achieved
            ztlvett_record_i = []
            if Z_j == Z_i:
                raise ValueError(
                    "Numerical accuracy exhausted in search of exit/burnout point."
                )
            try:
                if Z_j > Z_b:
                    Z_j = Z_b

                (
                    Z,
                    (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_1_j, tau_2_j),
                    _,
                ) = RKF78(
                    self._ode_Z,
                    (t_bar_i, l_bar_i, v_bar_i, eta_i, tau_1_i, tau_2_i),
                    Z_i,
                    Z_j,
                    relTol=tol,
                    absTol=tol**2,
                    minTol=minTol,
                    abortFunc=abort,
                    record=ztlvett_record_i,
                )

                p_bar_1_j = self._f_p_1_bar(Z_j, eta_j, tau_1_j)

            except ValueError as e:
                raise e

            if l_bar_j >= l_g_bar:
                if abs(l_bar_i - l_g_bar) / (l_g_bar) > tol or l_bar_i == 0:
                    N *= 2
                    Z_j = Z_i + Delta_Z / N
                else:
                    isBurnOutContained = False
                    break  # l_bar_i is solved to within a tol of l_bar_g

            else:
                ztlvett_record.extend(ztlvett_record_i)
                if p_bar_1_j > p_bar_max:
                    raise ValueError(
                        "Nobel-Abel EoS is generally accurate enough below 600MPa. However,"
                        + " Unreasonably high pressure (>{:.0f} MPa) was encountered.".format(
                            p_max / 1e6
                        )  # in practice most of the pressure-realted spikes are captured here.
                    )

                if v_bar_j <= 0:
                    Z, t_bar, l_bar, v_bar, eta, tau_1, tau_2 = (
                        ztlvett_record[-1][0],
                        *ztlvett_record[-1][1],
                    )

                    raise ValueError(
                        "Squib load condition detected: Shot stopped in bore.\n"
                        + "Shot is last calculated at {:.0f} mm at {:.0f} mm/s after {:.0f} ms".format(
                            l_bar * self.l_0 * 1e3,
                            v_bar * self.v_j * 1e3,
                            t_bar * tScale * 1e3,
                        )
                    )

                if any(v < 0 for v in (t_bar_j, l_bar_j, p_bar_1_j)):
                    raise ValueError(
                        "Numerical Integration diverged: negative"
                        + " values encountered in results.\n"
                        + "{:.0f} ms, {:.0f} mm, {:.0f} m/s, {:.0f} MPa".format(
                            t_bar_j * tScale * 1e3,
                            l_bar_j * self.l_0 * 1e3,
                            v_bar_j * self.v_j,
                            p_bar_1_j * pScale / 1e6,
                        )
                    )  # this will catch any case where t, l, p are negative

                # fmt: off
                t_bar_i, l_bar_i, v_bar_i, eta_i, tau_1_i, tau_2_i = (
                    t_bar_j, l_bar_j, v_bar_j, eta_j, tau_1_j, tau_2_j
                )  # fmt: on
                Z_i = Z_j
                """
                this way the group of values denoted by _i is always updated
                as a group.
                """
                Z_j += Delta_Z / N

        if t_bar_i == 0:
            raise ValueError("burnout point found to be at the origin.")

        # input()
        """
        Cludge code to force the SoE past the discontinuity at Z = Z_b, since
        we wrote the SoE to be be piecewise continous from (0, Z_b] and (Z_b, +inf)
        it is necessary to do this to prevent the RKF integrator coming up with
        irreducible error estimates and driving the step size to 0 around Z = Z_b
        """
        if isBurnOutContained:
            Z_i = Z_b + tol

        # fmt: off
        record.extend(
            (t_bar, (l_bar, self.f_psi_Z(Z), v_bar,
                     self._f_p_1_bar(Z, eta, tau_1),
                     self._f_p_2_bar(l_bar, eta, tau_2), eta, tau_1, tau_2))
            for (Z, (t_bar, l_bar, v_bar, eta, tau_1, tau_2)) in ztlvett_record
        )
        # fmt: on
        """
        Subscript e indicate exit condition.
        At this point, since its guaranteed that point i will be further
        towards the chamber of the firearm than point e, we do not have
        to worry about the dependency of the correct behaviour of this ODE
        on the positive direction of integration.
        """

        ltzvett_record = []
        (
            _,
            (t_bar_e, Z_e, v_bar_e, eta_e, tau_1_e, tau_2_e),
            (t_bar_err, Z_err, v_bar_err, eta_err, tau_1_err, tau_2_err),
        ) = RKF78(
            self._ode_l,
            (t_bar_i, Z_i, v_bar_i, eta_i, tau_1_i, tau_2_i),
            l_bar_i,
            l_g_bar,
            relTol=tol,
            absTol=tol**2,
            minTol=minTol,
            record=ltzvett_record,
        )
        # fmt: off
        record.extend(
            (t_bar, (l_bar, self.f_psi_Z(Z), v_bar,
                     self._f_p_1_bar(Z, eta, tau_1),
                     self._f_p_2_bar(l_bar, eta, tau_2), eta, tau_1, tau_2))
            for (l_bar, (t_bar, Z, v_bar, eta, tau_1, tau_2)) in ltzvett_record
        )
        # fmt: on
        updBarData(
            tag=POINT_EXIT,
            t_bar=t_bar_e,
            l_bar=l_g_bar,
            Z=Z_e,
            v_bar=v_bar_e,
            eta=eta_e,
            tau_1=tau_1_e,
            tau_2=tau_2_e,
            t_bar_err=t_bar_err,
            l_bar_err=0,
            Z_err=Z_err,
            v_bar_err=v_bar_err,
            eta_err=eta_err,
            tau_1_err=tau_1_err,
            tau_2_err=tau_2_err,
        )

        t_bar_f = None
        if Z_b > 1.0 and Z_e >= 1.0:  # fracture point exist and is contained
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile
            movement to charge fracture
            """
            # fmt: off
            (
                _,
                (t_bar_f, l_bar_f, v_bar_f, eta_f, tau_1_f, tau_2_f),
                (t_bar_err, l_bar_err, v_bar_err, eta_err_f, tau_1_err, tau_2_err),
            ) = RKF78(
                self._ode_Z,
                (t_bar_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                Z_1,
                1,
                relTol=tol,
                absTol=tol**2,
                minTol=minTol,
            )
            # fmt: on
            updBarData(
                tag=POINT_FRACTURE,
                t_bar=t_bar_f,
                l_bar=l_bar_f,
                Z=1,
                v_bar=v_bar_f,
                eta=eta_f,
                tau_1=tau_1_f,
                tau_2=tau_2_f,
                t_bar_err=t_bar_err,
                l_bar_err=l_bar_err,
                Z_err=0,
                v_bar_err=v_bar_err,
                eta_err=eta_err,
                tau_1_err=tau_1_err,
                tau_2_err=tau_2_err,
            )

        t_bar_b = None
        if isBurnOutContained:
            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to Z_b, from onset of projectile
            movement to charge burnout.
            """
            # fmt:off

            (
                _,
                (t_bar_b, l_bar_b, v_bar_b, eta_b, tau_1_b, tau_2_b),
                (t_bar_err, l_bar_err, v_bar_err, eta_err, tau_1_err, tau_2_err),
            ) = RKF78(
                self._ode_Z,
                (t_bar_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                Z_1,
                Z_b,
                relTol=tol,
                absTol=tol**2,
                minTol=minTol,
            )

            updBarData(
                tag=POINT_BURNOUT,
                t_bar=t_bar_b,
                l_bar=l_bar_b,
                Z=Z_b,
                v_bar=v_bar_b,
                eta=eta_b,
                tau_1=tau_1_b,
                tau_2=tau_2_b,
                t_bar_err=t_bar_err,
                l_bar_err=l_bar_err,
                Z_err=0,
                v_bar_err=v_bar_err,
                eta_err=eta_err,
                tau_1_err=tau_1_err,
                tau_2_err=tau_2_err
            )
            # fmt: on

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.

        we hereby simply solve p golden section searching it from origin
        to point e, i.e. inside the barrel.
        """

        def f(t, m="a"):
            Z, l_bar, v_bar, eta, tau_1, tau_2 = RKF78(
                self._ode_t,
                (Z_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                t_bar_1,
                t,
                relTol=tol,
                absTol=tol**2,
                minTol=minTol,
                # debug=True,
            )[1]
            p_bar = self._f_p_1_bar(Z, eta, tau_1)

            if m == "a":
                return p_bar

        def findPeak(f, tag):
            """
            tolerance is specified a bit differently for gold section search
            GSS tol is the length between the upper bound and lower bound
            of the maxima/minima, thus by our definition of tolerance (one sided),
            we take the median value.
            """

            t_bar_tol = tol * min(
                t for t in (t_bar_e, t_bar_b, t_bar_f) if t is not None
            )

            t_bar = 0.5 * sum(
                gss(
                    f,
                    t_bar_1,
                    t_bar_e if t_bar_b is None else t_bar_b,
                    x_tol=t_bar_tol,
                    findMin=False,
                )
            )

            # fmt: off
            (
                _,
                (Z, l_bar, v_bar, eta, tau_1, tau_2),
                (Z_err, l_bar_err, v_bar_err, eta_err, tau_1_err, tau_2_err),
            ) = RKF78(
                self._ode_t,
                (Z_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                t_bar_1,
                t_bar,
                relTol=tol,
                absTol=tol**2,
                minTol=minTol,
            )
            # fmt: on
            t_bar_err = 0.5 * t_bar_tol

            updBarData(
                tag=tag,
                t_bar=t_bar,
                l_bar=l_bar,
                Z=Z,
                v_bar=v_bar,
                eta=eta,
                tau_1=tau_1,
                tau_2=tau_2,
                t_bar_err=t_bar_err,
                l_bar_err=l_bar_err,
                Z_err=Z_err,
                v_bar_err=v_bar_err,
                eta_err=eta_err,
                tau_1_err=tau_1_err,
                tau_2_err=tau_2_err,
            )

        findPeak(lambda x: f(x, "a"), POINT_PEAK)

        """
        populate data for output purposes
        """
        try:
            if dom == DOMAIN_TIME:
                # fmt: off
                for j in range(step):
                    t_bar_j = t_bar_e / (step + 1) * (j + 1)
                    if t_bar_j < t_bar_1:
                        (
                            _,
                            (Z_j, eta_j, tau_1_j, tau_2_j),
                            (Z_err, eta_err, tau_1_err, tau_2_err),
                        ) = RKF78(
                            self._ode_ts,
                            (Z_0, 0, 1, 0),
                            0,
                            t_bar_j,
                            relTol=tol,
                            absTol=tol**2,
                            minTol=minTol,
                        )
                        # fmt: on
                        #print(self._f_p_2_bar(l_bar_j, eta_j, tau_2_j) * pScale)
                        l_bar_j, v_bar_j = 0, 0
                        l_bar_err, v_bar_err = 0, 0

                    else:
                        (
                            _,
                            (Z_j, l_bar_j, v_bar_j, eta_j, tau_1_j, tau_2_j),
                            (Z_err, l_bar_err, v_bar_err, eta_err, tau_1_err, tau_2_err),
                        ) = RKF78(
                            self._ode_t,
                            (Z_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                            t_bar_1,
                            t_bar_j,
                            relTol=tol,
                            absTol=tol**2,
                            minTol=minTol,
                        )
                        # fmt: on

                    updBarData(
                        tag="",
                        t_bar=t_bar_j,
                        l_bar=l_bar_j,
                        Z=Z_j,
                        v_bar=v_bar_j,
                        eta=eta_j,
                        tau_1=tau_1_j,
                        tau_2=tau_2_j,
                        t_bar_err=0,
                        l_bar_err=l_bar_err,
                        Z_err=Z_err,
                        v_bar_err=v_bar_err,
                        eta_err=eta_err,
                        tau_1_err=tau_1_err,
                        tau_2_err=tau_2_err,
                    )

            elif dom == DOMAIN_LENG:
                """
                Due to two issues, i.e. 1.the length domain ODE
                cannot be integrated from the origin point, and 2.the
                correct behaviour can only be expected when starting from
                a point with active burning else dZ flat lines.
                we do another Z domain integration to seed the initial values
                to a point where ongoing burning is guaranteed.
                (in the case of gun barrel length >= burn length, the group
                 of value by subscipt i will not guarantee burning is still
                 ongoing).
                """
                t_bar_j = 0.5 * t_bar_i
                Z_j, l_bar_j, v_bar_j, eta_j, tau_1_j, tau_2_j = RKF78(
                    self._ode_t,
                    (Z_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                    t_bar_1,
                    t_bar_j,
                    relTol=tol,
                    absTol=tol**2,
                    minTol=minTol,
                )[1]

                for j in range(step):
                    l_bar_k = l_g_bar / (step + 1) * (j + 1)

                    # fmt: off
                    (
                        _,
                        (t_bar_j, Z_j, v_bar_j, eta_j, tau_1_j, tau_2_j),
                        (t_bar_err, Z_err, v_bar_err, eta_err, tau_1_err, tau_2_err),
                    ) = RKF78(
                        self._ode_l,
                        (t_bar_j, Z_j, v_bar_j, eta_j, tau_1_j, tau_2_j),
                        l_bar_j,
                        l_bar_k,
                        relTol=tol,
                        absTol=tol**2,
                        minTol=minTol,
                    )  # fmt: on

                    l_bar_j = l_bar_k

                    updBarData(
                        tag="",
                        t_bar=t_bar_j,
                        l_bar=l_bar_j,
                        Z=Z_j,
                        v_bar=v_bar_j,
                        eta=eta_j,
                        tau_1=tau_1_j,
                        tau_2=tau_2_j,
                        t_bar_err=t_bar_err,
                        l_bar_err=0,
                        Z_err=Z_err,
                        v_bar_err=v_bar_err,
                        eta_err=eta_err,
                        tau_1_err=tau_1_err,
                        tau_2_err=tau_2_err,
                    )
            else:
                raise ValueError("Unknown domain")

        except ValueError as e:
            raise e
        finally:
            pass

        """
        sort the data points
        """

        data, error = [], []

        for bar_dataLine, bar_errorLine in zip(bar_data, bar_err):
            # fmt: off
            dtag, t_bar, l_bar, Z, v_bar, p_1_bar, p_2_bar, eta, tau_1, tau_2 = bar_dataLine
            (etag, t_bar_err, l_bar_err, Z_err, v_bar_err, p_1_bar_err, p_2_bar_err, eta_err,
             tau_1_err, tau_2_err) = bar_errorLine
            # fmt: on
            t, t_err = (t_bar * tScale, t_bar_err * tScale)
            l, l_err = (l_bar * self.l_0, l_bar_err * self.l_0)
            psi, psi_err = (self.f_psi_Z(Z), abs(self.f_sigma_Z(Z) * Z_err))
            v, v_err = v_bar * self.v_j, v_bar_err * self.v_j

            p_1, p_1_err = p_1_bar * pScale, (
                p_1_bar_err
                if isinstance(p_1_bar_err, str)
                else p_1_bar_err * pScale
            )
            p_2, p_2_err = p_2_bar * pScale, (
                p_2_bar_err
                if isinstance(p_2_bar_err, str)
                else p_2_bar_err * pScale
            )
            T_1 = tau_1 * self.T_v
            T_2 = tau_2 * self.T_v
            T_1_err = tau_1_err * self.T_v
            T_2_err = tau_2_err * self.T_v
            # fmt: off
            data.append((dtag, t, l, psi, v, p_1, p_2, eta, T_1, T_2))
            error.append(
                (etag, t_err, l_err, psi_err, v_err, p_1_err, p_2_err, eta_err,
                 T_1_err, T_2_err)
            )
            # fmt: on

        """
        scale the records too
        """
        # fmt: off
        for t_bar, (l_bar, psi, v_bar, p_1_bar, p_2_bar, eta, tau_1, tau_2) in record:
            data.append((
                "*", t_bar * tScale, l_bar * self.l_0, psi, v_bar * self.v_j,
                p_1_bar * pScale, p_2_bar * pScale,
                eta, tau_1 * self.T_v, tau_2 * self.T_v,
            ))

            errLine = ("L", *("--" for _ in range(9)))
            error.append(errLine)
        # fmt: ongoing
        data, error = zip(
            *(
                (a, b)
                for a, b in sorted(zip(data, error), key=lambda x: x[0][1])
            )
        )

        return data, error

    def getEff(self, vg):
        """
        te: thermal efficiency
        be: ballistic efficiency
        """
        te = (vg / self.v_j) ** 2
        be = te / self.phi
        return te, be


if __name__ == "__main__":
    """standard 7 port cylinder has d_0=e_1, port dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """
    from tabulate import tabulate

    compositions = GrainComp.readFile("data/propellants.csv")

    M17 = compositions["M17"]
    M1 = compositions["M1"]
    from prop import SimpleGeometry

    M17C = Propellant(M17, SimpleGeometry.CYLINDER, None, 2.5)
    M1C = Propellant(M1, SimpleGeometry.CYLINDER, None, 10)
    lf = 0.5
    print("DELTA/rho:", lf)
    test = Highlow(
        caliber=0.082,
        shotMass=2,
        propellant=M1C,
        grainSize=0.75e-3,
        chargeMass=0.3,
        chamberVolume=0.3 / M1C.rho_p / lf,
        expansionVolume=0.3 / M1C.rho_p / lf,
        startPressure=20e6,
        expansionStartPressure=10e6,
        lengthGun=3.5,
        nozzleExpansion=2.0,
        portAreaRatio=1,
        chambrage=1,
    )
    record = []

    print("\nnumerical: time")
    print(
        tabulate(
            test.integrate(100, 1e-4, dom=DOMAIN_TIME)[0],
            headers=(
                "tag",
                "t",
                "l",
                "psi",
                "v",
                "p1",
                "p2",
                "eta",
                "T1",
                "T2",
            ),
        )
    )
    print("\nnumerical: length")
    print(
        tabulate(
            test.integrate(10, 1e-4, dom=DOMAIN_LENG)[0],
            headers=(
                "tag",
                "t",
                "l",
                "psi",
                "v",
                "p1",
                "p2",
                "eta",
                "T1",
                "T2",
            ),
        )
    )
    # print(test.getEff(942))
