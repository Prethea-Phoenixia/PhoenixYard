from highlow import Highlow

# from pso import pso
from _coordesc import coordesc
from math import pi, log, inf
from num import cubic, RKF78, dekker, gss

# from gun import
from gun import POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_EXIT
from highlow import POINT_PEAK_HIGH, POINT_PEAK_BLEED


class ConstrainedHighlow:
    def __init__(
        self,
        caliber,
        propellant,
        shotMass,
        burstPressure,
        startPressure,
        dragCoefficient,
        chambrage,
        tol,
        designHighPressure,
        designLowPressure,
        designVelocity,
        minWeb=1e-6,
        maxLength=1e3,
        maxEV=1,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        control=POINT_PEAK_AVG,
        **_,
    ):
        # cache the constants:
        self.S = (0.5 * caliber) ** 2 * pi
        self.propellant = propellant
        self.m = shotMass
        self.p_0_e = burstPressure
        self.p_0_s = startPressure

        self.p_d_h = designHighPressure
        self.p_d_l = designLowPressure
        self.v_d = designVelocity

        self.chi_k = chambrage

        self.phi_1 = 1 / (1 - dragCoefficient)  # drag work coefficient

        self.tol = tol

        self.maxLength = maxLength
        self.maxEV = maxEV

        self.minWeb = minWeb
        self.ambientRho = ambientRho
        self.ambientP = ambientP
        self.ambientGamma = ambientGamma
        self.control = control

        theta = self.theta
        gamma = theta + 1
        self.phi_2 = 0.15  # for small ports 1.5mm to 2mm in size

        self.cfpr = (2 / (gamma + 1)) ** (gamma / theta)
        self.K_0 = gamma**0.5 * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * theta))

    def __getattr__(self, attrName):
        if "propellant" in vars(self) and not (
            attrName.startswith("__") and attrName.endswith("__")
        ):
            try:
                return getattr(self.propellant, attrName)
            except AttributeError:
                raise AttributeError(
                    "%r object has no attribute %r"
                    % (self.__class__.__name__, attrName)
                )
        else:
            raise AttributeError

    def solve(
        self,
        loadFraction,
        chargeMassRatio,
        # portAreaRatio,
        portArea,
        lengthGun=None,
        known_bore=False,
        suppress=False,
        progressQueue=None,
        **_,
    ):
        if any((chargeMassRatio <= 0, loadFraction <= 0, loadFraction > 1)):
            raise ValueError("Invalid parameters to solve constrained design problem")

        if progressQueue is not None:
            progressQueue.put(1)

        m = self.m
        rho_p = self.rho_p
        theta = self.theta
        f = self.f
        chi = self.chi
        # chi_k = self.chi_k
        labda = self.labda
        mu = self.mu
        S = self.S
        maxLF = self.maxLF
        phi_1 = self.phi_1
        p_d_l = self.p_d_l
        p_d_h = self.p_d_h

        u_1, n = self.u_1, self.n
        alpha = self.alpha
        Z_b = self.Z_b

        f_psi_Z = self.f_psi_Z
        f_sigma_Z = self.f_sigma_Z

        p_0_e = self.p_0_e
        p_0_s = self.p_0_s

        cfpr = self.cfpr
        phi_2 = self.phi_2
        K_0 = self.K_0

        Z_b = self.Z_b

        tol = self.tol

        p_a = self.ambientP
        if self.ambientRho != 0:
            c_a = (self.ambientGamma * self.ambientP / self.ambientRho) ** 0.5
        else:
            c_a = 0
        k_1 = self.ambientGamma
        control = self.control

        if loadFraction > maxLF:
            raise ValueError("Specified Load Fraction Violates Geometrical Constraint")

        omega = m * chargeMassRatio
        V_0 = omega / (rho_p * loadFraction)

        labda_1 = 1 / 2
        labda_2 = 1 / 3

        psi_0 = (
            p_0_e
            * (V_0 - omega / rho_p)
            / (f * omega - p_0_e * (omega / rho_p - alpha * omega))
        )

        Zs = cubic(chi * mu, chi * labda, chi, -psi_0)
        # pick a valid solution between 0 and 1
        Zs = sorted(
            Z for Z in Zs if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to overcome"
                + " port open pressure, or has burnt to post fracture."
            )

        def _f_p_1(Z, eta, tau_1, psi=None):
            psi = psi if psi else f_psi_Z(Z)
            V_psi = V_0 - omega / rho_p * (1 - psi) - alpha * omega * (psi - eta)
            return f * omega * tau_1 / V_psi * (psi - eta)

        # S_j_bar = portAreaRatio
        # S_j = S * S_j_bar
        S_j = portArea
        # S_j_bar = S_j / S

        phi = phi_1 + labda_2 * omega / m

        def f_e_1(e_1):
            Z_0, t_0, eta_0, tau_1_0 = Zs[0], 0, 0, 1

            def _ode_Zi(Z, t, eta, tau_1, _):
                psi = f_psi_Z(Z)
                p_1 = _f_p_1(Z, eta, tau_1, psi)
                dt = 1 / (u_1 / e_1 * p_1**n)  # dt / dZ
                deta = (phi_2 * K_0 * p_1 * S_j / ((f * tau_1) ** 0.5 * omega)) * dt
                dpsi = f_sigma_Z(Z)  # dpsi / dZ
                dtau_1 = ((1 - tau_1) * dpsi - theta * tau_1 * deta) / (psi - eta)

                return dt, deta, dtau_1

            dt_0, deta_0, dtau_1_0 = _ode_Zi(Z_0, t_0, eta_0, tau_1_0, None)

            dZ = Z_0 * tol
            Z_0 += dZ
            t_0 += dt_0 * dZ
            eta_0 += deta_0 * dZ
            tau_1_0 += dtau_1_0 * dZ

            record_0 = [[Z_0, [t_0, eta_0, tau_1_0]]]

            def abort_Zi(x, ys, record):
                Z, _, eta, tau_1 = x, *ys
                p_1 = _f_p_1(Z, eta, tau_1)
                if len(record) < 1:
                    return False
                o_x, o_ys = record[-1]
                oZ, _, oeta, otau_1 = o_x, *o_ys
                op_1 = _f_p_1(oZ, oeta, otau_1)

                return p_1 < op_1

            def g(Z):
                i = record_0.index([v for v in record_0 if v[0] <= Z][-1])

                x = record_0[i][0]
                ys = record_0[i][1]

                r = []
                Z, (t, eta, tau_1), _ = RKF78(
                    _ode_Zi,
                    ys,
                    x,
                    Z,
                    relTol=tol,
                    absTol=tol**2,
                    abortFunc=abort_Zi,
                    record=r,
                )

                if Z not in [line[0] for line in r]:
                    r.append([Z, (t, eta, tau_1)])

                xs = [v[0] for v in record_0]
                record_0.extend(v for v in r if v[0] not in xs)
                record_0.sort()

                return _f_p_1(Z, eta, tau_1)

            Z_p_i = (
                sum(
                    gss(
                        lambda Z: g(Z),
                        Z_0,
                        Z_b,
                        y_rel_tol=0.5 * tol,
                        findMin=False,
                    )
                )
                * 0.5
            )

            return g(Z_p_i) - p_d_h

        dp_bar_probe = f_e_1(self.minWeb)
        probeWeb = self.minWeb

        if dp_bar_probe < 0:
            raise ValueError(
                "Design pressure cannot be achieved by varying" + " web down to minimum"
            )

        while dp_bar_probe > 0:
            probeWeb *= 2
            dp_bar_probe = f_e_1(probeWeb)

        def fr(x):
            progressQueue.put(round(x * 50))

        e_1, _ = dekker(
            f_e_1,
            probeWeb,  # >0
            0.5 * probeWeb,  # ?0
            y_abs_tol=p_d_h * self.tol,
            f_report=fr if progressQueue is not None else None,
        )  # this is the e_1 that satisifies the pressure specification.

        # print("solved e_1 at ", e_1, "with a pressure deviation of ", f_e_1(e_1))

        def f_ev(expansionVolume):
            V_1 = expansionVolume

            def _f_p_2(l, eta, tau_2):
                l_star = (V_1 - alpha * omega * eta) / S
                p_avg = f * omega * tau_2 * eta / (S * (l_star + l))

                if control == POINT_PEAK_AVG:
                    return p_avg

                factor_s = 1 + labda_2 * (omega * eta / (phi_1 * m))

                if control == POINT_PEAK_SHOT:
                    return p_avg / factor_s

                factor_b = (phi_1 * m + labda_2 * omega * eta) / (
                    phi_1 * m + labda_1 * omega * eta
                )
                if control == POINT_PEAK_BLEED:
                    return p_avg / factor_b

            Z_0, t_0, eta_0, tau_1_0, tau_2_0 = Zs[0], 0, 0, 1, 1 + theta

            def _ode_Zs(Z, t, eta, tau_1, tau_2, _):
                psi = f_psi_Z(Z)
                p_1 = _f_p_1(Z, eta, tau_1, psi)

                l_star = (V_1 - alpha * omega * eta) / S
                p_2 = f * omega * tau_2 * eta / (S * l_star)

                dt = 1 / (u_1 / e_1 * p_1**n)  # dt / dZ

                pr = p_2 / p_1
                if pr <= cfpr:
                    deta = (phi_2 * K_0 * p_1 * S_j / ((f * tau_1) ** 0.5 * omega)) * dt
                else:
                    gamma = theta + 1
                    deta = (
                        (phi_2 * p_1 * S_j)
                        / ((f * tau_1) ** 0.5 * omega)
                        * (
                            (2 * gamma / theta)
                            * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                        )
                        ** 0.5
                    ) * dt

                dpsi = f_sigma_Z(Z)  # dpsi / dZ
                dtau_1 = ((1 - tau_1) * dpsi - theta * tau_1 * deta) / (psi - eta)

                if eta == 0:
                    dtau_2 = None
                else:
                    dtau_2 = (((1 + theta) * tau_1 - tau_2) * deta) / eta

                return dt, deta, dtau_1, dtau_2

            dt_0, deta_0, dtau_1_0, _ = _ode_Zs(Z_0, t_0, eta_0, tau_1_0, tau_2_0, None)

            dZ = Z_0 * tol
            Z_0 += dZ
            t_0 += dt_0 * dZ
            eta_0 += deta_0 * dZ
            tau_1_0 += dtau_1_0 * dZ

            record_0 = [[Z_0, [t_0, eta_0, tau_1_0, tau_2_0]]]

            def abort(x, ys, record):
                Z, _, eta, tau_1, tau_2 = x, *ys
                p_1 = _f_p_1(Z, eta, tau_1)

                l_star = (V_1 - alpha * omega * eta) / S
                p_2 = f * omega * tau_2 * eta / (S * l_star)

                if len(record) < 1:
                    return False

                return p_2 > p_0_s or (p_1 - p_2 < tol * p_1)

            def g(Z):
                i = record_0.index([v for v in record_0 if v[0] <= Z][-1])

                x = record_0[i][0]
                ys = record_0[i][1]

                r = []
                Z, (t, eta, tau_1, tau_2), _ = RKF78(
                    _ode_Zs,
                    ys,
                    x,
                    Z,
                    relTol=tol,
                    absTol=tol**2,
                    abortFunc=abort,
                    record=r,
                )

                if Z not in [line[0] for line in r]:
                    r.append([Z, (t, eta, tau_1)])

                xs = [v[0] for v in record_0]
                record_0.extend(v for v in r if v[0] not in xs)
                record_0.sort()

                l_star = (V_1 - alpha * omega * eta) / S
                p_2 = f * omega * tau_2 * eta / (S * l_star)

                return p_2

            Z_sm, (t, eta, tau_1, tau_2), _ = RKF78(
                _ode_Zs,
                record_0[0][1],
                record_0[0][0],
                Z_b,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort,
                record=record_0,
            )
            p_1_sm = _f_p_1(Z_sm, eta, tau_1)
            p_2_sm = _f_p_2(0, eta, tau_2)

            if abs(p_1_sm - p_2_sm) < tol * p_1_sm and p_2_sm < p_0_s:
                raise ValueError(
                    "Pressure levels have coalesced between high and low chamber "
                    + "before shot has started, "
                    + f"at ({p_2_sm * 1e-6:.6f} MPa)."
                )

            elif p_2_sm < p_0_s:
                raise ValueError(
                    f"Maximum pressure developed in low-chamber ({p_2_sm * 1e-6:.6f} MPa) "
                    + "not enough to start the shot."
                )

            Z_1, _ = dekker(lambda x: g(x) - p_0_s, Z_0, Z_sm, y_abs_tol=p_0_s * tol)

            print("Z_1 according to opt", Z_1)
            record_1 = [[Z_0, (t_0, eta_0, tau_1_0, tau_2_0)]]
            t_1, eta_1, tau_1_1, tau_2_1 = RKF78(
                _ode_Zs,
                (t_0, eta_0, tau_1_0, tau_2_0),
                Z_0,
                Z_1,
                relTol=tol,
                absTol=tol**2,
                record=record_1,
            )[1]

            def abort_Z(x, ys, record):
                Z, _, l, _, eta, tau_1, tau_2 = x, *ys
                p_h = _f_p_1(Z, eta, tau_1)
                p_l = _f_p_2(l, eta, tau_2)

                o_x, o_ys = record[-1]
                oZ, _, ol, _, oeta, otau_1, otau_2 = o_x, *o_ys
                op_h = _f_p_1(oZ, oeta, otau_1)
                op_l = _f_p_2(ol, oeta, otau_2)

                return p_h < op_h and p_l < op_l

            def _ode_Z(Z, t, l, v, eta, tau_1, tau_2, _):
                psi = f_psi_Z(Z)
                p_1 = _f_p_1(Z, eta, tau_1, psi)

                l_star = (V_1 - alpha * omega * eta) / S
                p_2 = f * omega * tau_2 * eta / (S * (l_star + l))

                dt = 1 / (u_1 / e_1 * p_1**n)  # dt / dZ

                pr = p_2 / p_1
                if pr <= cfpr:
                    deta = (
                        (phi_2 * K_0 * p_1 * S_j) / ((f * tau_1) ** 0.5 * omega)
                    ) * dt
                else:
                    gamma = theta + 1
                    deta = (
                        (phi_2 * p_1 * S_j)
                        / ((f * tau_1) ** 0.5 * omega)
                        * (
                            (2 * gamma / theta)
                            * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                        )
                        ** 0.5
                    ) * dt

                dpsi = f_sigma_Z(Z)  # dpsi / dZ
                dtau_1 = ((1 - tau_1) * dpsi - theta * tau_1 * deta) / (psi - eta)

                if c_a != 0 and v > 0:
                    k = k_1  # gamma
                    v_r = v / c_a
                    p_d = (
                        0.25 * k * (k + 1) * v_r**2
                        + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
                    ) * p_a
                else:
                    p_d = 0

                dv = S / (phi * m) * (p_2 - p_d) * dt
                dl = v * dt

                dtau_2 = (
                    ((1 + theta) * tau_1 - tau_2) * deta
                    - (theta * phi * m) / (f * omega) * v * dv
                ) / eta

                return dt, dl, dv, deta, dtau_1, dtau_2

            record_2 = [[Z_1, (t_1, 0, 0, eta_1, tau_1_1, tau_2_1)]]
            Z_j, (t_j, l_j, v_j, eta_j, tau_1_j, tau_2_j), _ = RKF78(
                _ode_Z,
                (t_1, 0, 0, eta_1, tau_1_1, tau_2_1),
                Z_1,
                Z_b,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort_Z,
                record=record_2,
            )

            def _f_p_2_Z(Z):
                i = record_2.index([v for v in record_2 if v[0] <= Z][-1])
                x, ys = record_2[i]

                r = []
                _, (t, l, v, eta, tau_1, tau_2), _ = RKF78(
                    dFunc=_ode_Z,
                    iniVal=ys,
                    x_0=x,
                    x_1=Z,
                    relTol=tol,
                    absTol=tol**2,
                    record=r,
                )
                xs = [v[0] for v in record_2]
                record_2.extend(v for v in r if v[0] not in xs)
                record_2.sort()
                return _f_p_2(l, eta, tau_2), (Z, t, l, v, eta, tau_1, tau_2)

            def _f_p_1_Z(Z):
                if Z < Z_1:
                    i = record_1.index([v for v in record_1 if v[0] <= Z][-1])
                    x, ys = record_1[i]

                    r = []
                    _, (t, eta, tau_1, tau_2), _ = RKF78(
                        dFunc=_ode_Zs,
                        iniVal=ys,
                        x_0=x,
                        x_1=Z,
                        relTol=tol,
                        absTol=tol**2,
                        record=r,
                    )
                    xs = [v[0] for v in record_1]
                    record_1.extend(v for v in r if v[0] not in xs)
                    record_1.sort()

                else:
                    i = record_2.index([v for v in record_2 if v[0] <= Z][-1])
                    x, ys = record_2[i]

                    r = []
                    _, (t, l, v, eta, tau_1, tau_2), _ = RKF78(
                        dFunc=_ode_Z,
                        iniVal=ys,
                        x_0=x,
                        x_1=Z,
                        relTol=tol,
                        absTol=tol**2,
                        record=r,
                    )
                    xs = [v[0] for v in record_2]
                    record_2.extend(v for v in r if v[0] not in xs)
                    record_2.sort()

                return _f_p_1(Z, eta, tau_1)

            Z_p_1 = (
                sum(
                    gss(
                        _f_p_1_Z,
                        Z_0,
                        Z_j,
                        y_rel_tol=0.5 * tol,
                        findMin=False,
                    )
                )
                * 0.5
            )

            Z_p_2 = (
                sum(
                    gss(
                        lambda Z: _f_p_2_Z(Z)[0],
                        Z_1,
                        Z_j,
                        y_rel_tol=0.5 * tol,
                        findMin=False,
                    )
                )
                * 0.5
            )

            p_p_h = _f_p_1_Z(Z_p_1)
            p_p_l, vals_2 = _f_p_2_Z(Z_p_2)

            return p_p_h, p_p_l, vals_2

        probeEV = self.maxEV
        lastEV = probeEV
        while True:
            try:
                dp_bar_probe = f_ev(probeEV)[1]
                break
            except ValueError:
                lastEV = probeEV
                probeEV *= 0.5

        # actual maximum lies somewhere between probeEV and 2x probeEV

        def findBound(func, x_probe, x_bound, tol, record=None):
            if record is None:
                record = []

            try:
                record.append((x_bound, func(x_bound)))
                return x_bound
            except ValueError:
                x_valid = x_probe

                delta = x_bound - x_probe
                up = x_bound > x_probe

                while abs(2 * delta) > tol and (
                    x_probe <= x_bound if up else x_probe >= x_bound
                ):
                    try:
                        record.append((x_probe, func(x_probe)))
                        x_valid = x_probe
                    except ValueError:
                        delta *= 0.5
                    finally:
                        x_probe = x_valid + delta

                return x_valid

        evUpperBound = findBound(f_ev, probeEV, lastEV, tol * V_0)

        p_l_min = f_ev(evUpperBound)[1]

        if p_l_min > p_d_l:
            raise ValueError(
                "Minimum pressure level in low chamber exceeds design pressure at "
                + f"{p_l_min * 1e-3:.2f} MPa and "
                + f"bounding volume of < {evUpperBound * 1e-3:.2f} L.",
            )

        probeEV = evUpperBound
        lastEV = probeEV
        while True:
            try:
                dp_bar_probe = f_ev(probeEV)[1]
                lastEV = probeEV
                probeEV *= 0.5
            except ValueError:
                break

        evLowerBound = findBound(f_ev, probeEV, lastEV, tol * V_0)
        p_l_max = f_ev(evLowerBound)[1]

        if p_l_max < p_d_l:
            raise ValueError(
                "Maximum pressure level in low chamber subceeds design pressure at "
                + f"{p_l_max * 1e-3:.2f} MPa and "
                + f"bounding volume of > {evLowerBound * 1e-3:.2f} L.",
            )

        def fr(x):
            progressQueue.put(round(x * 50) + 50)

        V_1, _ = dekker(
            lambda ev: f_ev(ev)[1],
            evLowerBound,
            evUpperBound,
            y=p_d_l,
            y_rel_tol=tol,
            f_report=fr if progressQueue is not None else None,
        )

        # print("expansion volumed solved to be ", V_1)

        p_h_act, p_l_act, (Z_i, t_i, l_i, v_i, eta_i, tau_1_i, tau_2_i) = f_ev(V_1)

        if abs(p_h_act - p_d_h) > tol * p_d_h:
            raise ValueError(
                "High and low pressure chamber conditions are insufficiently "
                + "decorrelated."
            )

        print("high pressure actual", p_h_act, "low pressure actual", p_l_act)

        if known_bore:
            if progressQueue is not None:
                progressQueue.put(100)
            return e_1, V_1, lengthGun

        # l_1 = V_1 / S

        def _f_p_2(l, eta, tau_2):
            l_star = (V_1 - alpha * omega * eta) / S
            p_avg = f * omega * tau_2 * eta / (S * (l_star + l))

            return p_avg

            if control == POINT_PEAK_AVG:
                return p_avg

            factor_s = 1 + labda_2 * (omega * eta / (phi_1 * m))

            if control == POINT_PEAK_SHOT:
                return p_avg / factor_s

            factor_b = (phi_1 * m + labda_2 * omega * eta) / (
                phi_1 * m + labda_1 * omega * eta
            )
            if control == POINT_PEAK_BLEED:
                return p_avg / factor_b

        def _ode_v(v, t, Z, l, eta, tau_1, tau_2, _):
            psi = f_psi_Z(Z)
            p_1 = _f_p_1(Z, eta, tau_1, psi)

            l_star = (V_1 - alpha * omega * eta) / S
            p_2 = f * omega * tau_2 * eta / (S * (l_star + l))

            if c_a != 0 and v > 0:
                k = k_1  # gamma
                v_r = v / c_a
                p_d = (
                    0.25 * k * (k + 1) * v_r**2
                    + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
                ) * p_a
            else:
                p_d = 0

            # dv = S / (phi * m) * (p_2 - p_d)  # dv / dt
            dt = (phi * m) / (S * (p_2 - p_d))  # dt / dv
            # dl = v # dl / dt
            dl = v * dt  # dl / dv

            if Z <= Z_b:
                dZ = u_1 / e_1 * p_1**n * dt  # dZ / dv
            else:
                dZ = 0

            pr = p_2 / p_1
            if pr <= cfpr:
                deta = (phi_2 * K_0 * p_1 * S_j) / ((f * tau_1) ** 0.5 * omega) * dt
            else:
                gamma = theta + 1
                deta = (
                    (phi_2 * p_1 * S_j)
                    / ((f * tau_1) ** 0.5 * omega)
                    * (
                        (2 * gamma / theta)
                        * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                    )
                    ** 0.5
                ) * dt

            dpsi = f_sigma_Z(Z) * dZ  # dpsi / dv
            dtau_1 = ((1 - tau_1) * dpsi - theta * tau_1 * deta) / (
                psi - eta
            )  # dtau_1 / dv

            dtau_2 = (
                ((1 + theta) * tau_1 - tau_2) * deta
                - (theta * phi * m) / (f * omega) * v
            ) / eta  # dtau_2 / dv

            return dt, dZ, dl, deta, dtau_1, dtau_2

        v_d = self.v_d
        if v_i > v_d and not suppress:
            raise ValueError(
                f"Design velocity exceeded ({v_i:.4g} m/s > {v_d:.4g} m/s) before peak pressure."
            )

        l_d = self.maxLength

        def abort_v(x, ys, record):
            _, _, l, _, _, _ = ys
            return l > l_d

        vtzlett_record = [v_i, (t_i, Z_i, l_i, eta_i, tau_1_i, tau_2_i)]
        try:
            v_g, (t_g, Z_g, l_g, eta_g, tau_1_g, tau_2_g), _ = RKF78(
                _ode_v,
                (t_i, Z_i, l_i, eta_i, tau_1_i, tau_2_i),
                v_i,
                v_d,
                relTol=tol,
                absTol=tol**2,
                record=vtzlett_record,
                abortFunc=abort_v,
            )

            p_1_g, p_2_g = _f_p_1(Z_g, eta_g, tau_1_g), _f_p_2(l_g, eta_g, tau_2_g)

        except ValueError:
            v_m, (t_m, Z_m, l_m, eta_m, tau_1_m, tau_2_m) = vtzlett_record[-1]
            p_1_m, p_2_m = _f_p_1(Z_m, eta_m, tau_1_m), _f_p_2(l_m, eta_m, tau_2_m)

            raise ValueError(
                "Integration appears to be approaching asymptote, "
                + f"last calculated to v = {v_m:.4g} m/s, "
                + f"x = {l_m:.4g} m, p = {p_1_m * 1e-6:.4g} MPa (high), {p_2_m * 1e-6:.4g} MPa (low) . "
                + "This indicates an excessive velocity target relative to pressure developed."
            )

        if l_g > l_d:
            raise ValueError(
                "Solution requires excessive tube length, last calculated to "
                + f"v = {v_g:.4g} m/s, x = {l_g:.4g} m, "
                + f"p = {p_1_g * 1e-6:.4g} MPa (high), {p_2_g * 1e-6:.4g} (low)."
            )

        if abs(v_g - v_d) > (tol * v_d):
            raise ValueError(
                "Velocity target is not met, last calculated to "
                + f"v = {v_g:.4g} m/s ({(v_g - v_d) / v_d:+.3g} %), x = {l_g:.4g} m, "
                + f"p = {p_1_g * 1e-6:.4g} MPa (high), {p_2_g * 1e-6:.4g} (low)."
            )

        if progressQueue is not None:
            progressQueue.put(100)

        return e_1, V_1, l_g


if __name__ == "__main__":
    from tabulate import tabulate
    from math import pi
    from prop import GrainComp, Propellant

    compositions = GrainComp.readFile("data/propellants.csv")

    M17 = compositions["M17"]
    M1 = compositions["M1"]
    from prop import SimpleGeometry

    M17C = Propellant(M17, SimpleGeometry.CYLINDER, None, 2.5)
    M1C = Propellant(M1, SimpleGeometry.CYLINDER, None, 10)
    test = ConstrainedHighlow(
        caliber=0.082,
        propellant=M1C,
        shotMass=5,
        burstPressure=10e6,
        startPressure=5e6,
        dragCoefficient=3e-3,
        chambrage=1,
        tol=1e-3,
        designHighPressure=50e6,
        designLowPressure=20e6,
        designVelocity=150,
        control=POINT_PEAK_BLEED,
    )

    # result = test.constrained(
    #     minWeb=1e-6,
    #     maxWeb=1e-2,
    #     minLF=0.1,
    #     minPortRatio=0.1,
    #     maxPortRatio=1 / 0.15,
    #     minLength=0.01,
    #     maxLength=1,
    #     minEV=0,
    #     maxEV=2e-3,
    #     control=POINT_PEAK_AVG,  # targeted pressure
    #     designHighPressure=80e6,
    #     designLowPressure=40e6,
    #     designVelocity=75,
    # )
    result = test.solve(
        loadFraction=0.5,
        chargeMassRatio=0.5 / 5,
        portArea=0.5 * pi * 0.082**2 * 0.25,
    )

    print(result)
