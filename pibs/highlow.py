from math import pi, log, inf
from num import gss, RKF78, cubic, bisect, dekker

from gun import Gun
from gun import DOMAIN_TIME, DOMAIN_LENG
from gun import (
    POINT_START,
    POINT_PEAK_AVG,
    POINT_PEAK_SHOT,
    POINT_FRACTURE,
    POINT_BURNOUT,
    POINT_EXIT,
)

POINT_PEAK_HIGH = "PEAK_HIGH_P"
POINT_PEAK_BLEED = "PEAK_BLEED_P"

HIGHLOW_PEAKS = [POINT_PEAK_HIGH, POINT_PEAK_BLEED, POINT_PEAK_AVG, POINT_PEAK_SHOT]


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
        burstPressure,  # low pressure chamber starting pressure
        startPressure,  # shot start pressure
        portArea,
        chambrage,
        lengthGun,
        dragCoefficient=0,
        structuralMaterial=None,
        structuralSafetyFactor=1.1,
        autofrettage=True,
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
                dragCoefficient < 0,
                dragCoefficient >= 1,
                burstPressure < 0,
                startPressure < 0,
            )
        ):
            raise ValueError("Invalid gun parameters")

        if chargeMass > (propellant.maxLF * propellant.rho_p * chamberVolume):
            raise ValueError("Specified Load Fraction Violates Geometrical Constraint")

        self.caliber = caliber
        self.propellant = propellant

        self.e_1 = 0.5 * grainSize
        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.omega = chargeMass

        self.V_0 = chamberVolume
        self.V_1 = expansionVolume

        self.p_0_e = burstPressure
        self.p_0_s = startPressure

        self.l_g = lengthGun
        self.chi_k = chambrage

        self.l_0 = self.V_0 / self.S
        self.l_1 = self.V_1 / self.S

        self.phi_1 = 1 / (1 - dragCoefficient)  # drag work coefficient
        self.phi = self.phi_1 + self.omega / (3 * self.m)

        # pressure ratio threshold for critical flow
        self.v_j = (2 * self.f * self.omega / (self.theta * self.phi * self.m)) ** 0.5

        self.material = structuralMaterial
        self.ssf = structuralSafetyFactor
        self.is_af = autofrettage

        if self.p_0_e == 0:
            raise NotImplementedError(
                "Current implementation use exponential burn rate and does not"
                + " allow for solving the case with 0 shot start pressure."
            )
        else:
            self.psi_0 = (
                self.p_0_e
                * (self.V_0 - self.omega / self.rho_p)
                / (
                    self.f * self.omega
                    - self.p_0_e * (self.omega / self.rho_p - self.alpha * self.omega)
                )
            )
            if self.psi_0 <= 0:
                raise ValueError(
                    "Initial burnup fraction is solved to be negative."
                    + " In practice this implies a detonation of the gun breech"
                    + " will likely occur."
                )

        Zs = cubic(self.chi * self.mu, self.chi * self.labda, self.chi, -self.psi_0)
        # pick a valid solution between 0 and 1
        Zs = sorted(
            Z for Z in Zs if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to overcome"
                + " port open pressure, or has burnt to post fracture."
            )

        self.Z_0 = Zs[0]  # pick the smallest solution

        self.S_j = portArea
        self.S_j_bar = portArea / self.S

        gamma = self.theta + 1
        self.phi_2 = 0.15  # for small ports 1.5mm to 2mm in size

        self.cfpr = (2 / (gamma + 1)) ** (gamma / self.theta)
        self.K_0 = gamma**0.5 * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * self.theta))

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

    def _f_p_1(self, Z, eta, tau_1, psi=None):
        psi = psi if psi else self.f_psi_Z(Z)
        V_psi = (
            self.V_0
            - self.omega / self.rho_p * (1 - psi)
            - self.alpha * self.omega * (psi - eta)
        )

        return self.f * self.omega * tau_1 / V_psi * (psi - eta)

    def _f_p_2(self, l, eta, tau_2):
        l_star = (self.V_1 - self.alpha * self.omega * eta) / self.S
        return self.f * self.omega * tau_2 * eta / (self.S * (l_star + l))

    def _ode_t(self, t, Z, l, v, eta, tau_1, tau_2, _):
        psi = self.f_psi_Z(Z)
        p_1 = self._f_p_1(Z, eta, tau_1, psi)

        if Z <= self.Z_b:
            dZ = self.u_1 / self.e_1 * p_1**self.n
        else:
            dZ = 0

        p_2 = self._f_p_2(l, eta, tau_2)

        pr = p_2 / p_1
        if pr <= self.cfpr:
            deta = (self.phi_2 * self.K_0 * p_1 * self.S_j) / (
                (self.f * tau_1) ** 0.5 * self.omega
            )
        else:
            gamma = self.theta + 1
            deta = (
                (self.phi_2 * p_1 * self.S_j)
                / ((self.f * tau_1) ** 0.5 * self.omega)
                * (
                    (2 * gamma / self.theta)
                    * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                )
                ** 0.5
            )

        dpsi = self.f_sigma_Z(Z) * dZ
        dtau_1 = ((1 - tau_1) * dpsi - self.theta * tau_1 * deta) / (psi - eta)

        if self.c_a != 0 and v > 0:
            k = self.k_1  # gamma
            v_r = v / self.c_a
            p_d = (
                0.25 * k * (k + 1) * v_r**2
                + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a
        else:
            p_d = 0

        dv = self.S / (self.phi * self.m) * (p_2 - p_d)
        dl = v

        dtau_2 = (
            ((1 + self.theta) * tau_1 - tau_2) * deta
            - (self.theta * self.phi * self.m) / (self.f * self.omega) * v * dv
        ) / eta

        return dZ, dl, dv, deta, dtau_1, dtau_2

    def _ode_ts(self, t, Z, eta, tau_1, tau_2, _):
        psi = self.f_psi_Z(Z)
        p_1 = self._f_p_1(Z, eta, tau_1, psi)

        if Z <= self.Z_b:
            dZ = self.u_1 / self.e_1 * p_1**self.n
        else:
            dZ = 0

        p_2 = self._f_p_2(0, eta, tau_2)

        pr = p_2 / p_1
        if pr <= self.cfpr:
            deta = (self.phi_2 * self.K_0 * p_1 * self.S_j) / (
                (self.f * tau_1) ** 0.5 * self.omega
            )
        else:
            gamma = self.theta + 1
            deta = (
                (self.phi_2 * p_1 * self.S_j)
                / ((self.f * tau_1) ** 0.5 * self.omega)
                * (
                    (2 * gamma / self.theta)
                    * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                )
                ** 0.5
            )

        dpsi = self.f_sigma_Z(Z) * dZ
        dtau_1 = ((1 - tau_1) * dpsi - self.theta * tau_1 * deta) / (psi - eta)

        dtau_2 = (((1 + self.theta) * tau_1 - tau_2) * deta) / eta
        return dZ, deta, dtau_1, dtau_2

    def _ode_l(self, l, t, Z, v, eta, tau_1, tau_2, _):
        dt = 1 / v  # dt / dl

        psi = self.f_psi_Z(Z)
        p_1 = self._f_p_1(Z, eta, tau_1, psi)

        if Z <= self.Z_b:
            dZ = self.u_1 / self.e_1 * p_1**self.n * dt  # dZ / dl
        else:
            dZ = 0

        p_2 = self._f_p_2(l, eta, tau_2)

        pr = p_2 / p_1
        if pr <= self.cfpr:
            deta = (
                (self.phi_2 * self.K_0 * p_1 * self.S_j)
                / ((self.f * tau_1) ** 0.5 * self.omega)
            ) * dt  # deta / dl
        else:
            gamma = self.theta + 1
            deta = (
                (self.phi_2 * p_1 * self.S_j)
                / ((self.f * tau_1) ** 0.5 * self.omega)
                * (
                    (2 * gamma / self.theta)
                    * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                )
                ** 0.5
            ) * dt  # deta / dl

        dpsi = self.f_sigma_Z(Z) * dZ  # dpsi / dl
        dtau_1 = ((1 - tau_1) * dpsi - self.theta * tau_1 * deta) / (psi - eta)

        if self.c_a != 0 and v > 0:
            k = self.k_1  # gamma
            v_r = v / self.c_a
            p_d = (
                0.25 * k * (k + 1) * v_r**2
                + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a
        else:
            p_d = 0

        dv = self.S / (self.phi * self.m) * (p_2 - p_d) * dt

        dtau_2 = (
            ((1 + self.theta) * tau_1 - tau_2) * deta
            - (self.theta * self.phi * self.m) / (self.f * self.omega) * v * dv
        ) / eta

        return dt, dZ, dv, deta, dtau_1, dtau_2

    def _ode_Z(self, Z, t, l, v, eta, tau_1, tau_2, _):
        psi = self.f_psi_Z(Z)
        p_1 = self._f_p_1(Z, eta, tau_1, psi)

        dt = 1 / (self.u_1 / self.e_1 * p_1**self.n)  # dt / dZ

        p_2 = self._f_p_2(l, eta, tau_2)

        pr = p_2 / p_1
        if pr <= self.cfpr:
            deta = (
                (self.phi_2 * self.K_0 * p_1 * self.S_j)
                / ((self.f * tau_1) ** 0.5 * self.omega)
            ) * dt
        else:
            gamma = self.theta + 1
            deta = (
                (self.phi_2 * p_1 * self.S_j)
                / ((self.f * tau_1) ** 0.5 * self.omega)
                * (
                    (2 * gamma / self.theta)
                    * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                )
                ** 0.5
            ) * dt

        dpsi = self.f_sigma_Z(Z)  # dpsi / dZ
        dtau_1 = ((1 - tau_1) * dpsi - self.theta * tau_1 * deta) / (psi - eta)

        if self.c_a != 0 and v > 0:
            k = self.k_1  # gamma
            v_r = v / self.c_a
            p_d = (
                0.25 * k * (k + 1) * v_r**2
                + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a
        else:
            p_d = 0

        dv = self.S / (self.phi * self.m) * (p_2 - p_d) * dt
        dl = v * dt

        dtau_2 = (
            ((1 + self.theta) * tau_1 - tau_2) * deta
            - (self.theta * self.phi * self.m) / (self.f * self.omega) * v * dv
        ) / eta

        return dt, dl, dv, deta, dtau_1, dtau_2

    def _ode_Zs(self, Z, t, eta, tau_1, tau_2, _):
        psi = self.f_psi_Z(Z)
        p_1 = self._f_p_1(Z, eta, tau_1, psi)

        dt = 1 / (self.u_1 / self.e_1 * p_1**self.n)  # dt / dZ

        p_2 = self._f_p_2(0, eta, tau_2)

        pr = p_2 / p_1
        if pr <= self.cfpr:
            deta = (
                self.phi_2
                * self.K_0
                * p_1
                * self.S_j
                / ((self.f * tau_1) ** 0.5 * self.omega)
            ) * dt
        else:
            gamma = self.theta + 1
            deta = (
                (self.phi_2 * p_1 * self.S_j)
                / ((self.f * tau_1) ** 0.5 * self.omega)
                * (
                    (2 * gamma / self.theta)
                    * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                )
                ** 0.5
            ) * dt

        dpsi = self.f_sigma_Z(Z)  # dpsi / dZ
        dtau_1 = ((1 - tau_1) * dpsi - self.theta * tau_1 * deta) / (psi - eta)

        if eta == 0:
            dtau_2 = None
        else:
            dtau_2 = (((1 + self.theta) * tau_1 - tau_2) * deta) / eta

        return dt, deta, dtau_1, dtau_2

    def integrate(
        self,
        step=10,
        tol=1e-5,
        dom=DOMAIN_TIME,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        peaks=HIGHLOW_PEAKS,
        progressQueue=None,
        **_,
    ):
        if progressQueue is not None:
            progressQueue.put(1)
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

        if any((step < 0, tol < 0)):
            raise ValueError("Invalid integration specification")

        if any((ambientP < 0, ambientRho < 0, ambientGamma < 1)):
            raise ValueError("Invalid ambient condition")

        p_max = 1e9  # 1GPa
        # ambient conditions
        self.p_a = ambientP

        if ambientRho != 0:
            self.c_a = (ambientGamma * ambientP / ambientRho) ** 0.5
        else:
            self.c_a = 0

        self.k_1 = ambientGamma

        Z_b = self.Z_b
        Z_0 = self.Z_0

        l_g = self.l_g

        bar_data = []
        bar_err = []

        # fmt: off
        def updBarData(
            tag="*", t=0, l=0, Z=0, v=0, eta=0, tau_1=0, tau_2=0,
            t_err=0, l_err=0, Z_err=0, v_err=0, eta_err=0, tau_1_err=0, tau_2_err=0,
        ):
            p_1, p_2 = (self._f_p_1(Z, eta, tau_1), self._f_p_2(l, eta, tau_2))
            p_err = "N/A"
            bar_data.append((tag, t, l, Z, v, p_1, p_2, eta, tau_1, tau_2))
            bar_err.append(
                (
                    "L", t_err, l_err, Z_err, v_err, p_err, p_err, eta_err,
                    tau_1_err, tau_2_err
                )
            )
        # fmt: on

        t_0, eta_0, tau_1_0, tau_2_0 = 0, 0, 1, 1 + self.theta

        updBarData(t=t_0, l=0, Z=Z_0, v=0, eta=eta_0, tau_1=tau_1_0, tau_2=tau_2_0)

        dt_0, deta_0, dtau_1_0, dtau_2_0 = self._ode_Zs(
            Z_0, t_0, eta_0, tau_1_0, tau_2_0, _
        )

        dZ = Z_0 * tol
        Z_0 += dZ
        t_0 += dt_0 * dZ
        eta_0 += deta_0 * dZ
        tau_1_0 += dtau_1_0 * dZ

        tett_record = [[Z_0, [t_0, eta_0, tau_1_0, tau_2_0]]]

        def abort(x, ys, record):
            Z, _, eta, tau_1, tau_2 = x, *ys

            p_1 = self._f_p_1(Z, eta, tau_1)
            p_2 = self._f_p_2(0, eta, tau_2)

            if len(record) < 1:
                return False

            ox, oys = record[-1]
            oZ, _, oeta, otau_1, _ = ox, *oys
            op_1 = self._f_p_1(oZ, oeta, otau_1)
            # op_2 = self._f_p_2(0, oeta, otau_2)

            return (
                p_2 > self.p_0_s
                or p_1 > p_max
                or (p_1 - p_2 < tol * p_1 and p_1 > op_1)
            )

        def f(Z):
            i = tett_record.index([v for v in tett_record if v[0] <= Z][-1])

            x = tett_record[i][0]
            ys = tett_record[i][1]

            r = []
            _, (t, eta, tau_1, tau_2), _ = RKF78(
                self._ode_Zs,
                ys,
                x,
                Z,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort,
                record=r,
            )

            xs = [v[0] for v in tett_record]
            tett_record.extend(v for v in r if v[0] not in xs)
            tett_record.sort()
            p_2 = self._f_p_2(0, eta, tau_2)

            return p_2

        Z, (t, eta, tau_1, tau_2), _ = RKF78(
            self._ode_Zs,
            tett_record[0][1],
            tett_record[0][0],
            Z_b,
            relTol=tol,
            absTol=tol**2,
            abortFunc=abort,
            record=tett_record,
        )
        p_1_sm = self._f_p_1(Z, eta, tau_1)
        p_2_sm = self._f_p_2(0, eta, tau_2)

        ox, oys = tett_record[-1]
        oZ, _, oeta, otau_1, _ = ox, *oys
        op_1_sm = self._f_p_1(oZ, oeta, otau_1)

        if (
            abs(p_1_sm - p_2_sm) < tol * p_1_sm
            and p_2_sm < self.p_0_s
            and p_1_sm > op_1_sm
        ):
            raise ValueError(
                "Pressure levels have coalesced between high and low chamber "
                + "before shot has started, "
                + f"at {p_1_sm * 1e-6:.6f} MPa (high), {p_2_sm * 1e-6:.6f} MPa (low)."
            )

        elif p_1_sm > p_max:
            raise ValueError(
                "Nobel-Abel EoS is generally accurate enough below 600MPa. However, "
                + f"Unreasonably high pressure ({p_1_sm * 1e-6:.0f}>{p_max * 1e-6:.0f}"
                + " MPa) was encountered.",  # in practice most of the pressure-realted spikes are captured here.
            )

        elif p_2_sm < self.p_0_s:
            raise ValueError(
                f"Maximum pressure developed in low-chamber ({p_2_sm * 1e-6:.6f} MPa) "
                + "not enough to start the shot."
            )

        Z_1, _ = dekker(lambda x: f(x) - self.p_0_s, Z_0, Z, y_abs_tol=self.p_0_s * tol)

        # fmt: off
        record.extend(
            (
                t,
                (
                    0, self.f_psi_Z(Z), 0, self._f_p_1(Z, eta, tau_1),
                    self._f_p_2(0, eta, tau_2), eta, tau_1, tau_2
                )
            )
            for (Z, (t, eta, tau_1, tau_2)) in tett_record
            if Z < Z_1
        )
        # fmt: on

        t_1, eta_1, tau_1_1, tau_2_1 = RKF78(
            self._ode_Zs,
            (t_0, eta_0, tau_1_0, tau_2_0),
            Z_0,
            Z_1,
            relTol=tol,
            absTol=tol**2,
        )[1]

        # fmt: off
        updBarData(
            tag=POINT_START, t=t_1, l=0, Z=Z_1, v=0, eta=eta_1, tau_1=tau_1_1, tau_2=tau_2_1
        )
        # fmt: on

        if progressQueue is not None:
            progressQueue.put(5)

        Z_i = Z_1
        Z_j = Z_b
        N = 1
        Delta_Z = Z_b - Z_0

        t_i, l_i, v_i, eta_i, tau_1_i, tau_2_i = t_1, 0, 0, eta_1, tau_1_1, tau_2_1

        isBurnOutContained = True

        """
        Instead of letting the integrator handle the heavy lifting, we
        partition Z and integrate upwards until either barrel exit or
        burnout has been achieved. This seek-and-validate inspired
        from bisection puts the group of value denoted by subscript
        i within or on the muzzle, with the propellant either still
        burning or right on the burnout point..
        """
        ztlvett_record = [(Z_1, (t_1, 0, 0, eta_1, tau_1_1, tau_2_1))]

        def abort(x, ys, record):
            Z, _, l, v, eta, tau_1, _ = x, *ys
            p_1 = self._f_p_1(Z, eta, tau_1)

            return any((l > l_g, p_1 > p_max, v <= 0, p_1 < 0))

        while Z_i < Z_b:  # terminates if burnout is achieved
            ztlvett_record_i = []
            if Z_j == Z_i:
                raise ValueError(
                    "Numerical accuracy exhausted in search of exit/burnout point."
                )
            try:
                if Z_j > Z_b:
                    Z_j = Z_b

                Z, (t_j, l_j, v_j, eta_j, tau_1_j, tau_2_j), _ = RKF78(
                    self._ode_Z,
                    (t_i, l_i, v_i, eta_i, tau_1_i, tau_2_i),
                    Z_i,
                    Z_j,
                    relTol=tol,
                    absTol=tol**2,
                    abortFunc=abort,
                    record=ztlvett_record_i,
                )

                p_1_j = self._f_p_1(Z_j, eta_j, tau_1_j)

            except ValueError as e:
                raise e

            if l_j >= l_g:
                if abs(l_i - l_g) > tol * l_g or l_i == 0:
                    N *= 2
                    Z_j = Z_i + Delta_Z / N
                else:
                    isBurnOutContained = False
                    break  # l_i is solved to within a tol of l_g

            else:
                ztlvett_record.extend(ztlvett_record_i)
                if p_1_j > p_max:
                    raise ValueError(
                        "Nobel-Abel EoS is generally accurate enough below 600MPa. However,"
                        + " Unreasonably high pressure (>{:.0f} MPa) was encountered.".format(
                            p_max * 1e-6
                        )  # in practice most of the pressure-realted spikes are captured here.
                    )

                if v_j <= 0:
                    Z, t, l, v, eta, tau_1, tau_2 = (
                        ztlvett_record[-1][0],
                        *ztlvett_record[-1][1],
                    )

                    raise ValueError(
                        "Squib load condition detected: Shot stopped in bore.\n"
                        + "Shot is last calculated at {:.0f} mm at {:.0f} mm/s after {:.0f} ms".format(
                            l * 1e3, v * 1e3, t * 1e3
                        )
                    )

                if any(v < 0 for v in (t_j, l_j, p_1_j)):
                    raise ValueError(
                        "Numerical Integration diverged: negative"
                        + " values encountered in results.\n"
                        + "{:.0f} ms, {:.0f} mm, {:.0f} mm/s, {:.0f} MPa".format(
                            t_j * 1e3, l_j * 1e3, v_j * 1e3, p_1_j * 1e-6
                        )
                    )  # this will catch any case where t, l, p are negative

                # fmt: off
                (
                    t_i, l_i, v_i, eta_i, tau_1_i, tau_2_i
                ) = (
                    t_j, l_j, v_j, eta_j, tau_1_j, tau_2_j
                )
                # fmt: on
                Z_i = Z_j
                """
                this way the group of values denoted by _i is always updated
                as a group.
                """
                Z_j += Delta_Z / N

        if t_i == 0:
            raise ValueError("burnout point found to be at the origin.")

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
            (
                t,
                (
                    l, self.f_psi_Z(Z), v, self._f_p_1(Z, eta, tau_1),
                    self._f_p_2(l, eta, tau_2), eta, tau_1, tau_2
                ),
            )
            for (Z, (t, l, v, eta, tau_1, tau_2)) in ztlvett_record
        )
        # fmt: on

        if progressQueue is not None:
            progressQueue.put(10)

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
            (t_e, Z_e, v_e, eta_e, tau_1_e, tau_2_e),
            (t_err, Z_err, v_err, eta_err, tau_1_err, tau_2_err),
        ) = RKF78(
            self._ode_l,
            (t_i, Z_i, v_i, eta_i, tau_1_i, tau_2_i),
            l_i,
            l_g,
            relTol=tol,
            absTol=tol**2,
            record=ltzvett_record,
        )
        # fmt: off
        record.extend(
            (
                t,
                (
                    l, self.f_psi_Z(Z), v, self._f_p_1(Z, eta, tau_1),
                    self._f_p_2(l, eta, tau_2), eta, tau_1, tau_2
                ),
            )
            for (l, (t, Z, v, eta, tau_1, tau_2)) in ltzvett_record
        )

        updBarData(
            tag=POINT_EXIT, t=t_e, l=l_g, Z=Z_e, v=v_e, eta=eta_e, tau_1=tau_1_e,
            tau_2=tau_2_e, t_err=t_err, l_err=0, Z_err=Z_err, v_err=v_err, eta_err=eta_err,
            tau_1_err=tau_1_err, tau_2_err=tau_2_err
        )
        # fmt: on

        if progressQueue is not None:
            progressQueue.put(20)

        t_f = None
        if Z_b > 1.0 and Z_e >= 1.0:  # fracture point exist and is contained
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile
            movement to charge fracture
            """

            (
                _,
                (t_f, l_f, v_f, eta_f, tau_1_f, tau_2_f),
                (t_err, l_err, v_err, eta_err, tau_1_err, tau_2_err),
            ) = RKF78(
                self._ode_Z,
                (t_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                Z_1,
                1,
                relTol=tol,
                absTol=tol**2,
            )
            # fmt: off
            updBarData(
                tag=POINT_FRACTURE, t=t_f, l=l_f, Z=1, v=v_f, eta=eta_f,
                tau_1=tau_1_f, tau_2=tau_2_f, t_err=t_err,
                l_err=l_err, Z_err=0, v_err=v_err, eta_err=eta_err,
                tau_1_err=tau_1_err, tau_2_err=tau_2_err,
            )
            # fmt: on

        t_b = None
        if isBurnOutContained:
            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to Z_b, from onset of projectile
            movement to charge burnout.
            """

            (
                _,
                (t_b, l_b, v_b, eta_b, tau_1_b, tau_2_b),
                (t_err, l_err, v_err, eta_err, tau_1_err, tau_2_err),
            ) = RKF78(
                self._ode_Z,
                (t_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                Z_1,
                Z_b,
                relTol=tol,
                absTol=tol**2,
            )
            # fmt:off
            updBarData(
                tag=POINT_BURNOUT, t=t_b, l=l_b, Z=Z_b, v=v_b, eta=eta_b,
                tau_1=tau_1_b, tau_2=tau_2_b, t_err=t_err, l_err=l_err, Z_err=0,
                v_err=v_err, eta_err=eta_err, tau_1_err=tau_1_err, tau_2_err=tau_2_err
            )
            # fmt: on

        if progressQueue is not None:
            progressQueue.put(30)

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.

        we hereby simply solve p golden section searching it from origin
        to point e, i.e. inside the barrel.
        """

        def g(t, m=POINT_PEAK_AVG):
            if t < t_1:
                Z, eta, tau_1, tau_2 = RKF78(
                    self._ode_ts,
                    (Z_0, eta_0, tau_1_0, tau_2_0),
                    t_0,
                    t,
                    relTol=tol,
                    absTol=tol**2,
                )[1]

                l, v = 0, 0
            else:
                Z, l, v, eta, tau_1, tau_2 = RKF78(
                    self._ode_t,
                    (Z_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                    t_1,
                    t,
                    relTol=tol,
                    absTol=tol**2,
                )[1]

            if m == POINT_PEAK_HIGH:
                p_high = self._f_p_1(Z, eta, tau_1)
                return p_high
            p_low = self._f_p_2(l, eta, tau_2)
            if m == POINT_PEAK_AVG:
                return p_low
            p_s, p_b = self.toPsPb(l, p_low, eta)
            if m == POINT_PEAK_SHOT:
                return p_s
            elif m == POINT_PEAK_BLEED:
                return p_b

        def findPeak(g, tag):
            """
            tolerance is specified a bit differently for gold section search
            GSS tol is the length between the upper bound and lower bound
            of the maxima/minima, thus by our definition of tolerance (one sided),
            we take the median value.
            """

            t_tol = tol * min(t for t in (t_e, t_b, t_f) if t is not None)

            t = 0.5 * sum(
                gss(g, t_0, t_e if t_b is None else t_b, x_tol=t_tol, findMin=False)
            )

            if t < t_1:
                l, v = 0, 0
                l_err, v_err = 0, 0
                (
                    _,
                    (Z, eta, tau_1, tau_2),
                    (Z_err, eta_err, tau_1_err, tau_2_err),
                ) = RKF78(
                    self._ode_ts,
                    (Z_0, eta_0, tau_1_0, tau_2_0),
                    t_0,
                    t,
                    relTol=tol,
                    absTol=tol**2,
                )

            else:
                # fmt: off
                (
                    _,
                    (Z, l, v, eta, tau_1, tau_2),
                    (Z_err, l_err, v_err, eta_err, tau_1_err, tau_2_err),
                ) = RKF78(
                    self._ode_t,
                    (Z_1, 0, 0, eta_1, tau_1_1, tau_2_1), t_1, t, relTol=tol, absTol=tol**2
                )
                # fmt: on

            t_err = 0.5 * t_tol

            # fmt: off
            updBarData(
                tag=tag, t=t, l=l, Z=Z, v=v, eta=eta, tau_1=tau_1, tau_2=tau_2,
                t_err=t_err, l_err=l_err, Z_err=Z_err, v_err=v_err, eta_err=eta_err,
                tau_1_err=tau_1_err, tau_2_err=tau_2_err
            )  # fmt:on

        for i, peak in enumerate(peaks):
            findPeak(lambda x: g(x, peak), peak)
            if progressQueue is not None:
                progressQueue.put(30 + i / (len(peaks) - 1) * 50)

        """
        populate data for output purposes
        """
        try:
            if dom == DOMAIN_TIME:
                # fmt: off
                (
                    t_j, Z_j, eta_j, tau_1_j, tau_2_j
                ) = (
                    t_0, Z_0, eta_0, tau_1_0, tau_2_0
                )
                # fmt: on
                l_j, v_j = 0, 0
                l_err, v_err = 0, 0

                hasStarted = False
                for j in range(step):
                    t_k = t_e / (step + 1) * (j + 1)

                    if t_k < t_1:
                        (
                            _,
                            (Z_j, eta_j, tau_1_j, tau_2_j),
                            (Z_err, eta_err, tau_1_err, tau_2_err),
                        ) = RKF78(
                            self._ode_ts,
                            (Z_j, eta_j, tau_1_j, tau_2_j),
                            t_j,
                            t_k,
                            relTol=tol,
                            absTol=tol**2,
                        )

                    else:
                        if not hasStarted:
                            # fmt: off
                            (
                                t_j, Z_j, eta_j, tau_1_j, tau_2_j
                            ) = (
                                t_1, Z_1, eta_1, tau_1_1, tau_2_1
                            )
                            hasStarted = True
                            # fmt: on

                        (
                            _,
                            (Z_j, l_j, v_j, eta_j, tau_1_j, tau_2_j),
                            (Z_err, l_err, v_err, eta_err, tau_1_err, tau_2_err),
                        ) = RKF78(
                            self._ode_t,
                            (Z_j, l_j, v_j, eta_j, tau_1_j, tau_2_j),
                            t_j,
                            t_k,
                            relTol=tol,
                            absTol=tol**2,
                        )

                    t_j = t_k

                    # fmt: off
                    updBarData(
                        tag="", t=t_j, l=l_j, Z=Z_j, v=v_j, eta=eta_j, tau_1=tau_1_j,
                        tau_2=tau_2_j, t_err=t_err, l_err=0, Z_err=Z_err, v_err=v_err,
                        eta_err=eta_err, tau_1_err=tau_1_err, tau_2_err=tau_2_err
                    )
                    # fmt: on

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

                t_j = 0.5 * (t_i - t_1) + t_1
                Z_j, l_j, v_j, eta_j, tau_1_j, tau_2_j = RKF78(
                    self._ode_t,
                    (Z_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                    t_1,
                    t_j,
                    relTol=tol,
                    absTol=tol**2,
                )[1]

                for j in range(step):
                    l_k = l_g / (step + 1) * (j + 1)

                    (
                        _,
                        (t_j, Z_j, v_j, eta_j, tau_1_j, tau_2_j),
                        (t_err, Z_err, v_err, eta_err, tau_1_err, tau_2_err),
                    ) = RKF78(
                        self._ode_l,
                        (t_j, Z_j, v_j, eta_j, tau_1_j, tau_2_j),
                        l_j,
                        l_k,
                        relTol=tol,
                        absTol=tol**2,
                        record=[],
                    )

                    l_j = l_k

                    # fmt: off
                    updBarData(
                        tag="", t=t_j, l=l_j, Z=Z_j, v=v_j, eta=eta_j, tau_1=tau_1_j,
                        tau_2=tau_2_j, t_err=t_err, l_err=0, Z_err=Z_err, v_err=v_err,
                        eta_err=eta_err, tau_1_err=tau_1_err, tau_2_err=tau_2_err
                    )
                    # fmt: on
            else:
                raise ValueError("Unknown domain")

        except ValueError as e:
            raise e
        finally:
            pass

        if progressQueue is not None:
            progressQueue.put(100)

        """
        sort the data points
        """

        data, error = [], []

        for bar_dataLine, bar_errorLine in zip(bar_data, bar_err):
            # fmt: off
            dtag, t, l, Z, v, p_1, p_2, eta, tau_1, tau_2 = bar_dataLine
            (etag, t_err, l_err, Z_err, v_err, p_1_err, p_2_err, eta_err,
             tau_1_err, tau_2_err) = bar_errorLine
            # fmt: on
            psi, psi_err = self.f_psi_Z(Z), abs(self.f_sigma_Z(Z) * Z_err)
            p_s, p_b = self.toPsPb(l, p_2, eta)

            T_1 = tau_1 * self.T_v
            T_2 = tau_2 * self.T_v
            T_1_err = tau_1_err * self.T_v
            T_2_err = tau_2_err * self.T_v
            # fmt: off
            data.append((dtag, t, l, psi, v, p_1, p_b, p_2, p_s, T_1, T_2, eta))
            error.append(
                (etag, t_err, l_err, psi_err, v_err, p_1_err, "---", p_2_err, "---",
                 T_1_err, T_2_err, eta_err)
            )
            # fmt: on

        """
        scale the records too
        """

        for t, (l, psi, v, p_1, p_2, eta, tau_1, tau_2) in record:
            if t in [line[1] for line in data]:
                continue
            p_s, p_b = self.toPsPb(l, p_2, eta)
            # fmt: off
            data.append(
                (
                    "*", t, l, psi, v, p_1, p_b, p_2, p_s, tau_1 * self.T_v,
                    tau_2 * self.T_v, eta
                )
            )
            # fmt: on

            errLine = ("L", *("--" for _ in range(11)))
            error.append(errLine)

        data, error = zip(
            *((a, b) for a, b in sorted(zip(data, error), key=lambda x: x[0][1]))
        )

        # calculate a pressure tracing.
        p_trace = []

        l_h = self.l_0 / self.chi_k  # physical length of the high pressure chamber
        l_l = self.l_1 / self.chi_k  # low pressure chamber

        for line in data:
            tag, t, l, psi, v, p_h, p_b, p, p_s, T_1, T_2, eta = line
            p_line = []
            for i in range(step):
                x = i / step * (l + l_l) + l_h
                p_x = self.toPx(l, p_h, p_b, p_s, x)
                p_line.append((x, p_x))

            p_line.append((l + l_l + l_h, p_s))
            p_trace.append((tag, psi, T_2, p_line))
            p_trace.append((tag, psi, T_1, [(0, p_h), (l_h * (1 - tol), p_h)]))

        try:
            if self.material is None:
                raise ValueError("Structural material not specified")

            structure = self.getStructural(data, step, tol)

        except Exception:
            import sys, traceback

            exc_type, exc_value, exc_traceback = sys.exc_info()
            errMsg = "".join(
                traceback.format_exception(exc_type, exc_value, exc_traceback)
            )
            print(errMsg)
            structure = [None, None, None]

        return data, error, p_trace, structure

    def getEff(self, vg):
        """
        te: thermal efficiency
        be: ballistic efficiency
        """
        te = (vg / self.v_j) ** 2
        be = te / self.phi
        return te, be

    def toPsPb(self, l, p, eta):
        """
        Convert average chamber pressure at a certain travel to
        shot base pressure, and breech face pressure

        eta: fraction of propellant gas in the low pressure chamber
        l: travel of the projectile
        p: average pressure

        Ps: pressure at shot
        Pb: pressure at breech
        """
        # Labda_g = l / self.l_1
        # labda_1_prime = (1 / 2) * (1 / self.chi_k + Labda_g) / (1 + Labda_g)
        # labda_2_prime = (1 / 3) * (1 / self.chi_k + Labda_g) / (1 + Labda_g)
        # factor_s = 1 + labda_2_prime * (self.omega * eta / (self.phi_1 * self.m))

        # factor_b = (self.phi_1 * self.m + labda_2_prime * self.omega * eta) / (
        #     self.phi_1 * self.m + labda_1_prime * self.omega * eta
        # )

        labda_1, labda_2 = 1 / 2, 1 / 3
        factor_s = 1 + labda_2 * (self.omega * eta / (self.phi_1 * self.m))

        factor_b = (self.phi_1 * self.m + labda_2 * self.omega * eta) / (
            self.phi_1 * self.m + labda_1 * self.omega * eta
        )

        return p / factor_s, p / factor_b

    def toPx(self, l, p_h, p_b, p_s, x):
        """
        | l_h | l_l |======l======
        """
        l_h = self.l_0 / self.chi_k  # physical length of the high pressure chamber
        l_l = self.l_1 / self.chi_k  # low pressure chamber

        if x < l_h:
            p_x = p_h
        else:
            r = (
                self.chi_k * (x - l_h)
                if x < (l_l + l_h)
                else (x - l_l - l_h) + self.l_1
            )
            k = (r / (self.l_1 + l)) ** 2

            p_x = p_s * k + p_b * (1 - k)

        return p_x

    def getStructural(self, data, step, tol):
        l_g = self.l_g
        chi_k = self.chi_k
        l_h = self.l_0 / chi_k  # physical length of high chamber
        l_l = self.l_1 / chi_k  # physical length of low chamber

        sigma = self.material.Y
        S = self.S

        r_s = 0.5 * self.caliber
        r_b = r_s * chi_k**0.5  # radius of breech
        S_b = S * chi_k  # area of breech

        x_probes = (
            [i / step * l_h for i in range(step)]
            + [l_h * (1 - tol)]
            + [i / step * l_l + l_h for i in range(step)]
            + [l_h + l_l * (1 - tol)]
            + [i / step * l_g + (l_h + l_l) for i in range(step)]
            + [l_g + l_h + l_l]
        )
        p_probes = [0] * len(x_probes)

        for line in data:
            tag, t, l, psi, v, p_h, p_b, p, p_s, T_1, T_2, eta = line

            for i, x in enumerate(x_probes):
                if (x - (l_h + l_l)) <= l:
                    p_x = self.toPx(l, p_h, p_b, p_s, x)
                    p_probes[i] = max(p_probes[i], p_x)
                else:
                    break

        for i, p in enumerate(p_probes):
            x = x_probes[i]
            p_probes[i] = p * self.ssf

        rho_probes = []
        V = 0
        if self.is_af:
            i, j = x_probes.index(l_h), x_probes.index(l_h + l_l)

            x_h, p_h = x_probes[:i], p_probes[:i]
            x_l, p_l = x_probes[i:j], p_probes[i:j]
            x_b, p_b = x_probes[j:], p_probes[j:]

            V_l, rho_l = Gun._Vrho_k(x_l, p_l, [S_b for _ in x_l], sigma, tol)
            V_h, rho_h = Gun._Vrho_k(
                x_h, p_h, [S_b for _ in x_h], sigma, tol, k_min=rho_l[0]
            )

            V_b, rho_b = Gun._Vrho_k(
                x_b,
                p_b,
                [S for _ in x_b],
                sigma,
                tol,
                p_ref=max(p_l),
                k_max=rho_l[-1] * chi_k**0.5,
            )

            V = V_h + V_l + V_b
            rho_probes = rho_h + rho_l + rho_b
        else:
            for p in p_probes:
                y = p / sigma

                if y > 3**-0.5:
                    raise ValueError(
                        f"Limit to conventional construction ({sigma * 3*1e-6:.3f} MPa)"
                        + " exceeded in section."
                    )
                rho = ((1 + y * (4 - 3 * y**2) ** 0.5) / (1 - 3 * y**2)) ** 0.5
                rho_probes.append(rho)

            for i in range(len(x_probes) - 1):
                x_0 = x_probes[i]
                x_1 = x_probes[i + 1]
                rho_0 = rho_probes[i]
                rho_1 = rho_probes[i + 1]
                dV = (rho_1**2 + rho_0**2 - 2) * 0.5 * S * (x_1 - x_0)
                if x_1 < (l_h + l_l):
                    V += dV * chi_k
                else:
                    V += dV

        tube_mass = V * self.material.rho

        hull = []
        for x, rho in zip(x_probes, rho_probes):
            if x < l_h + l_l:
                hull.append((x, r_b, rho * r_b))
            else:
                hull.append((x, r_s, rho * r_s))

        return tube_mass, None, hull


if __name__ == "__main__":
    """standard 7 port cylinder has d_0=e_1, port dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """
    from tabulate import tabulate
    from prop import GrainComp, Propellant

    compositions = GrainComp.readFile("data/propellants.csv")

    M17 = compositions["M17"]
    M1 = compositions["M1"]
    from prop import SimpleGeometry

    M17C = Propellant(M17, SimpleGeometry.CYLINDER, None, 25)
    M1C = Propellant(M1, SimpleGeometry.CYLINDER, None, 10)
    lf = 0.5
    print("DELTA/rho:", lf)
    test = Highlow(
        caliber=0.082,
        shotMass=5,
        propellant=M17C,
        grainSize=5e-3,
        chargeMass=0.5,
        chamberVolume=0.5 / M1C.rho_p / lf,
        expansionVolume=0.5 / M1C.rho_p / lf,
        startPressure=5e6,
        burstPressure=10e6,
        lengthGun=3.5,
        portArea=0.5 * (pi * 0.082**2 * 0.25),
        chambrage=1,
    )
    record = []

    print("\nnumerical: time")
    print(
        tabulate(
            test.integrate(10, 1e-3, dom=DOMAIN_TIME)[0],
            headers=(
                "tag",
                "t",
                "l",
                "psi",
                "v",
                "p1",
                "pb",
                "p2",
                "ps",
                "T1",
                "T2",
                "eta",
            ),
        )
    )

    # input()

    print("\nnumerical: length")
    print(
        tabulate(
            test.integrate(9, 1e-3, dom=DOMAIN_LENG)[0],
            headers=(
                "tag",
                "t",
                "l",
                "psi",
                "v",
                "p1",
                "pb",
                "p2",
                "ps",
                "T1",
                "T2",
                "eta",
            ),
        )
    )
    # print(test.getEff(942))
