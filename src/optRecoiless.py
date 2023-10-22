from num import gss, RKF78, cubic, dekker
from prop import Propellant
from random import uniform
from math import pi, inf
from recoiless import Recoiless
from opt import N
from bisect import bisect
from gun import POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT


class ConstrainedRecoiless:
    def __init__(
        self,
        caliber,
        shotMass,
        propellant,
        startPressure,
        dragCoefficient,
        designPressure,
        designVelocity,
        nozzleExpansion,
        nozzleEfficiency,
        chambrage,
        **_,
    ):
        # constants for constrained designs

        if any(
            (
                caliber <= 0,
                shotMass <= 0,
                startPressure <= 0,
                dragCoefficient < 0,
                dragCoefficient >= 1,
                nozzleExpansion < 1,
                nozzleEfficiency > 1,
                nozzleEfficiency <= 0,
                chambrage < 1,
            )
        ):
            raise ValueError("Invalid parameters for constrained design")

        if any(
            (
                designPressure <= 0,
                # designPressure <= startPressure,
                designVelocity <= 0,
            )
        ):
            raise ValueError("Invalid design constraint")

        self.S = (0.5 * caliber) ** 2 * pi
        self.m = shotMass
        self.propellant = propellant
        self.p_0 = startPressure
        self.phi_1 = 1 / (1 - dragCoefficient)

        # design limits
        self.p_d = designPressure
        self.v_d = designVelocity

        self.chi_0 = nozzleEfficiency
        self.A_bar = nozzleExpansion

        self.chi_k = chambrage

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
        tol,
        minWeb=1e-6,
        containBurnout=True,
        maxLength=1e3,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        control=POINT_PEAK_AVG,
        **_,
    ):
        if any(
            (
                minWeb <= 0,
                tol <= 0,
                chargeMassRatio <= 0,
                loadFraction <= 0,
                loadFraction > 1,
                maxLength <= 0,
                ambientRho < 0,
                ambientP < 0,
                ambientGamma < 1,
            )
        ):
            raise ValueError(
                "Invalid parameters to solve constrained design problem"
            )
        """
        minWeb  : represents minimum possible grain size
        """

        m = self.m
        rho_p = self.rho_p
        theta = self.theta
        f = self.f
        chi = self.chi
        chi_k = self.chi_k
        labda = self.labda
        mu = self.mu
        S = self.S
        maxLF = self.maxLF
        phi_1 = self.phi_1
        p_0 = self.p_0
        v_d = self.v_d
        p_d = self.p_d
        u_1 = self.u_1
        n = self.n
        alpha = self.alpha
        Z_b = self.Z_b

        f_psi_Z = self.f_psi_Z
        f_sigma_Z = self.f_sigma_Z

        chi_0 = self.chi_0
        A_bar = self.A_bar

        Sx = S * chi_k

        if loadFraction > maxLF:
            raise ValueError(
                "Specified Load Fraction Violates Geometrical Constraint"
            )

        omega = m * chargeMassRatio
        V_0 = omega / (rho_p * loadFraction)
        Delta = omega / V_0
        l_0 = V_0 / S

        gamma = theta + 1

        phi = phi_1 + omega / (3 * m)

        S_j_bar = 1 / (Recoiless.getCf(gamma, A_bar, tol) * chi_0)
        if S_j_bar > chi_k:
            raise ValueError(
                "Achieving recoiless condition necessitates"
                + " a larger throat area than could be fit into"
                + " breech face."
            )
        S_j = S_j_bar * S

        K_0 = (2 / (gamma + 1)) ** (
            (gamma + 1) / (2 * (gamma - 1))
        ) * gamma**0.5

        phi_2 = 1
        C_A = (
            (0.5 * theta * phi * m / omega) ** 0.5 * K_0 * phi_2
        )  # flow rate value

        """
        it is impossible to account for the chambrage effect given unspecified
        barrel length, in our formulation
        """
        v_j = (2 * f * omega / (theta * phi * m)) ** 0.5

        if ambientRho != 0:
            c_a_bar = (ambientGamma * ambientP / ambientRho) ** 0.5 / v_j
            p_a_bar = ambientP / (f * Delta)
        else:
            c_a_bar = 0
            p_a_bar = 0

        gamma_1 = ambientGamma

        if v_j < v_d:
            raise ValueError(
                "Propellant load too low to achieve design velocity. "
                + " The 2nd ballistic limit for this loading conditions is"
                + " {:.4g} m/s,".format(v_j)
                + " and recoiless guns only achieve a part of that as well."
            )

        psi_0 = (1 / Delta - 1 / rho_p) / (f / p_0 + alpha - 1 / rho_p)

        Zs = cubic(
            a=chi * mu,
            b=chi * labda,
            c=chi,
            d=-psi_0,
        )
        # pick a valid solution between 0 and 1
        Zs = tuple(
            Z
            for Z in Zs
            if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to overcome"
                + " start pressure, or has burnt to post fracture."
            )
        Z_0 = Zs[0]

        # p_bar_0 = p_0 / (f * Delta)

        def _f_p_bar(Z, l_bar, v_bar, eta, tau):
            psi = f_psi_Z(Z)
            l_psi_bar = 1 - Delta * ((1 - psi) / rho_p + alpha * (psi - eta))
            p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

            if control == POINT_PEAK_AVG:
                return p_bar

            else:
                y = omega * eta
                m_dot = C_A * v_j * S_j * p_bar * Delta / (tau**0.5)
                vx = m_dot * (V_0 + S * l_bar * l_0) / (Sx * (omega - y))

                H_1, H_2 = (
                    vx / (v_j * v_bar) if v_bar != 0 else inf,
                    2 * phi_1 * m / (omega - y) + 1,
                )

                H = min(H_1, H_2)

                p_s_bar = p_bar / (
                    1 + (omega - y) / (3 * phi_1 * m) * (1 - 0.5 * H)
                )
                if control == POINT_PEAK_SHOT:
                    return p_s_bar
                elif control == POINT_PEAK_BREECH:
                    return (
                        p_s_bar * (1 + (omega - y) / (2 * phi_1 * m) * (1 - H))
                        if H == H_1
                        else 0
                    )

        p_bar_d = p_d / (f * Delta)  # convert to unitless
        l_bar_d = maxLength / l_0
        """
        if p_bar_d < p_bar_s:
            raise ValueError("Pressure constraint lower than starting value.")
        """
        """
        step 1, find grain size that satisifies design pressure
        """

        def abort_Z(x, ys, record):
            Z, t_bar, l_bar, v_bar, eta, tau = x, *ys
            p_bar = _f_p_bar(Z, l_bar, v_bar, eta, tau)
            return (p_bar > 2 * p_bar_d) or l_bar > l_bar_d

        def _f_p_e_1(e_1, tol=tol):
            """
            Find pressure maximum bracketed by the range of
            Z_0 < Z < Z_d
            l_bar < l_bar_d,
            p_bar < 2 * p_bar_d.
            """
            # print("e_1", e_1)
            B = (
                S**2
                * e_1**2
                / (f * phi * omega * m * u_1**2)
                * (f * Delta) ** (2 * (1 - n))
            )

            # integrate this to end of burn

            def _ode_Z(Z, t_bar, l_bar, v_bar, eta, tau, _):
                """burnout domain ode of internal ballistics"""
                psi = f_psi_Z(Z)
                dpsi = f_sigma_Z(Z)  # dpsi/dZ

                l_psi_bar = 1 - Delta * (
                    (1 - psi) / rho_p + alpha * (psi - eta)
                )
                p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

                if c_a_bar != 0 and v_bar > 0:
                    v_r = v_bar / c_a_bar
                    p_d_bar = (
                        +0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                        + gamma_1
                        * v_r
                        * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                    ) * p_a_bar
                else:
                    p_d_bar = 0

                if Z <= Z_b:
                    dt_bar = (2 * B / theta) ** 0.5 * p_bar**-n
                    dl_bar = v_bar * dt_bar
                    dv_bar = 0.5 * theta * (p_bar - p_d_bar) * dt_bar

                else:
                    # technically speaking it is undefined in this area
                    dt_bar = 0  # dt_bar/dZ
                    dl_bar = 0  # dl_bar/dZ
                    dv_bar = 0  # dv_bar/dZ

                deta = C_A * S_j_bar * p_bar / tau**0.5 * dt_bar  # deta / dZ
                dtau = (
                    (1 - tau) * (dpsi)
                    - 2 * v_bar * (dv_bar)
                    - theta * tau * (deta)
                ) / (psi - eta)

                return (dt_bar, dl_bar, dv_bar, deta, dtau)

            # stepVanished = False
            record = [[Z_0, [0, 0, 0, 0, 1]]]
            try:
                (
                    Z_j,
                    (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j),
                    _,
                ) = RKF78(
                    dFunc=_ode_Z,
                    iniVal=(0, 0, 0, 0, 1),
                    x_0=Z_0,
                    x_1=Z_b,
                    relTol=tol,
                    absTol=tol**2,
                    abortFunc=abort_Z,
                    record=record,
                    # debug=True,
                )

                if Z_j not in [line[0] for line in record]:
                    record.append(
                        [Z_j, [t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j]]
                    )
            except ValueError:
                Z_j, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j) = record[-1]
                # stepVanished = True

            p_bar_j = _f_p_bar(Z_j, l_bar_j, v_bar_j, eta_j, tau_j)

            """print(
                *[
                    (Z, _f_p_bar(Z, l_bar, v_bar, eta, tau) * f * Delta)
                    for (Z, (t_bar, l_bar, v_bar, eta, tau)) in record
                ],
                sep="\n",
            )"""

            # print(*record, sep="\n")

            """
            find the peak pressure point.
            """
            Z_i = Z_0
            peak = None
            if len(record) > 2:
                p_bars = [
                    _f_p_bar(Z, l_bar, v_bar, eta, tau)
                    for (Z, (t_bar, l_bar, v_bar, eta, tau)) in record
                ]

                for i, (l, c, r) in enumerate(
                    zip(p_bars[:-2], p_bars[1:-1], p_bars[2:])
                ):
                    if l < c and c > r:
                        peak = i + 1

            if peak is None:
                # no peak, so it suffice to compare the end points.
                if p_bar_j == 0:
                    p_bar_j = inf
                return (p_bar_j - p_bar_d, record[-1][0], *record[-1][-1])
            else:  # peak exist, must compare the peak and the two end points.
                Z_i = record[peak - 1][0]
                Z_j = record[peak + 1][0]

                def _f_p_Z(Z):
                    i = record.index([v for v in record if v[0] <= Z][-1])
                    x = record[i][0]
                    ys = record[i][1]

                    r = []
                    _, (t_bar, l_bar, v_bar, eta, tau), _ = RKF78(
                        dFunc=_ode_Z,
                        iniVal=ys,
                        x_0=x,
                        x_1=Z,
                        relTol=tol,
                        absTol=tol**2,
                        record=r,
                    )

                    xs = [v[0] for v in record]
                    record.extend(v for v in r if v[0] not in xs)
                    record.sort()
                    return _f_p_bar(Z, l_bar, v_bar, eta, tau)

                Z_1, Z_2 = gss(
                    _f_p_Z,
                    Z_i,
                    Z_j,
                    y_rel_tol=tol,
                    findMin=False,
                )
                Z_p = 0.5 * (Z_1 + Z_2)

                if abs(Z_p - Z_b) < tol:
                    Z_p = Z_b

                p_bar_p = _f_p_Z(Z_p)
                i = [line[0] for line in record].index(Z_p)

                p_bar_max = max((p_bar_p, p_bar_j))
                if p_bar_p == p_bar_max:
                    return (p_bar_p - p_bar_d, record[i][0], *record[i][-1])
                else:
                    return (p_bar_j - p_bar_d, record[-1][0], *record[-1][-1])

        dp_bar_probe, Z, *_ = _f_p_e_1(minWeb)
        probeWeb = minWeb

        if dp_bar_probe < 0:
            raise ValueError(
                "Design pressure cannot be achieved by varying web down to minimum. "
                + "Peak pressure found at phi = {:.4g} at {:.4g} MPa".format(
                    f_psi_Z(Z), (dp_bar_probe + p_bar_d) * 1e-6 * f * Delta
                )
            )

        while dp_bar_probe > 0:
            probeWeb *= 2
            dp_bar_probe = _f_p_e_1(probeWeb)[0]

        e_1, e_1_2 = dekker(
            lambda web: _f_p_e_1(web)[0],
            probeWeb,  # >0
            0.5 * probeWeb,  # ?0
            # x_tol=1e-14,
            y_abs_tol=p_bar_d * tol,
            # debug=True,
        )  # this is the e_1 that satisifies the pressure specification.

        """
        e_1 and e_2 brackets the true solution
        """

        dp_bar_i, *vals_1 = _f_p_e_1(e_1)
        dp_bar_j, *vals_2 = _f_p_e_1(e_1_2)

        (Z_i, t_bar_i, l_bar_i, v_bar_i, eta_i, tau_i) = vals_1

        if abs(dp_bar_i) > tol * p_bar_d:
            raise ValueError(
                "Design pressure is not met, current best solution peaked at "
                + "P = {:.4g}MPa ({:+.3g}%) ".format(
                    (dp_bar_i + p_bar_d) * (f * Delta * 1e-6),
                    dp_bar_i / p_bar_d * 100,
                )
                + "and the bracketing solution at "
                + "P = {:.4g}MPa ({:+.3g}%),".format(
                    (dp_bar_j + p_bar_d) * (f * Delta * 1e-6),
                    dp_bar_j / p_bar_d * 100,
                )
                + " with a web difference of {:.4g} mm.".format(
                    (e_1 - e_1_2) * 1e3
                )
                + " If the pressures are too disparate this may indicate desired"
                + " solution lies at the edge or outside of solution space."
            )

        if abs(Z_i - Z_b) < tol:
            """
            fudging the starting Z value for velocity integration to prevent
            driving integrator to 0 at the transition point.
            """
            Z_i = Z_b + tol

        """
        step 2, find the requisite muzzle length to achieve design velocity
        """
        v_bar_d = v_d / v_j

        if v_bar_i > v_bar_d and containBurnout:
            raise ValueError(
                "Design velocity exceeded ({:.4g} m/s > {:.4g} m/s) before peak pressure.".format(
                    v_bar_i * v_j, v_bar_d * v_j
                )
            )

        # TODO: find some way of making this cross constraint less troublesome.
        B = (
            S**2
            * e_1**2
            / (f * phi * omega * m * u_1**2)
            * (f * Delta) ** (2 * (1 - n))
        )

        def _ode_v(v_bar, t_bar, Z, l_bar, eta, tau, _):
            psi = f_psi_Z(Z)
            dpsi = f_sigma_Z(Z)  # dpsi/dZ

            l_psi_bar = 1 - Delta * ((1 - psi) / rho_p + alpha * (psi - eta))
            p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

            if c_a_bar != 0 and v_bar > 0:
                v_r = v_bar / c_a_bar
                p_d_bar = (
                    +0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                    + gamma_1
                    * v_r
                    * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                ) * p_a_bar
            else:
                p_d_bar = 0

            dt_bar = 2 / (theta * (p_bar - p_d_bar))

            if Z <= Z_b:
                dZ = dt_bar * (0.5 * theta / B) ** 0.5 * p_bar**n
            else:
                dZ = 0

            dl_bar = v_bar * dt_bar

            deta = C_A * S_j_bar * p_bar / tau**0.5 * dt_bar  # deta / dv_bar
            dtau = (
                (1 - tau) * (dpsi * dZ)  # dZ/dt_bar
                - 2 * v_bar  # dv_bar/dt_bar
                - theta * tau * (deta)  # deta/dt_bar
            ) / (
                psi - eta
            )  # dtau/dt_bar

            return (dt_bar, dZ, dl_bar, deta, dtau)

        def abort_v(x, ys, record):
            # v_bar = x
            t_bar, _, l_bar, _, _ = ys

            ot_bar, *_ = record[-1][-1]
            return l_bar > l_bar_d or t_bar < ot_bar

        try:
            vtzlet_record = [[v_bar_i, (t_bar_i, Z_i, l_bar_i, eta_i, tau_i)]]
            (
                v_bar_g,
                (t_bar_g, Z_g, l_bar_g, eta, tau),
                _,
            ) = RKF78(
                dFunc=_ode_v,
                iniVal=(t_bar_i, Z_i, l_bar_i, eta_i, tau_i),
                x_0=v_bar_i,
                x_1=v_bar_d,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort_v,
                record=vtzlet_record,  # debug
                # debug=True,
            )

        except ValueError:
            vmax = vtzlet_record[-1][0] * v_j
            lmax = vtzlet_record[-1][-1][2] * l_0
            raise ValueError(
                "Velocity asymptotically approaching {:.4g} m/s at {:.4g} m.".format(
                    vmax, lmax
                )
                + " This indicates an excessive velocity target relative to pressure developed."
            )

        if l_bar_g > l_bar_d:
            raise ValueError(
                "Solution requires excessive tube length. Last calculated at {:.4g} m/s at {:.4g}.".format(
                    v_bar_g * v_j, l_bar_g * l_0, maxLength
                )
            )

        elif t_bar_g < t_bar_i:
            raise ValueError(
                "Projectile is already decelerating at peak pressure point ({:.4g} m/s at {:.4g} m).".format(
                    v_bar_i * v_j, l_bar_i * l_0
                )
                + " This indicate excessively low propulsive effort compared to drag"
                + " in bore. Reduce the caliber, or specify more energetic propellant charge."
            )
        if abs(v_bar_g - v_bar_d) > (v_bar_d * tol):
            raise ValueError(
                "Velocity specification is not met: {:.4g} m/s ({:+.3g} %)".format(
                    v_bar_g * v_j, (v_bar_g - v_bar_d) / v_bar_d * 1e2
                )
            )

        return e_1, l_bar_g * l_0

    def findMinV(
        self,
        chargeMassRatio,
        tol,
        minWeb,
        maxLength,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        loadFraction=None,  # hint
        control=POINT_PEAK_AVG,
        **_,
    ):
        """
        find the minimum volume solution.
        """

        """
        Step 1, find a valid range of values for load fraction,
        using psuedo-bisection.

        high lf -> high web
        low lf -> low web
        """
        omega = self.m * chargeMassRatio
        rho_p = self.rho_p
        S = self.S
        solve = self.solve

        def f(lf):
            V_0 = omega / (rho_p * lf)
            l_0 = V_0 / S

            e_1, l_g = solve(
                loadFraction=lf,
                chargeMassRatio=chargeMassRatio,
                tol=tol,  # this is to ensure unimodality up to ~tol
                minWeb=minWeb,
                containBurnout=False,
                maxLength=maxLength,
                ambientP=ambientP,
                ambientRho=ambientRho,
                ambientGamma=ambientGamma,
                control=control,
            )
            return e_1, (l_g + l_0), l_g

        records = []

        for i in range(N):
            startProbe = (
                loadFraction
                if (i == 1 and loadFraction is not None)
                else uniform(tol, 1 - tol)
            )
            try:
                _, lt_i, lg_i = f(startProbe)
                records.append((startProbe, lt_i))
                break
            except ValueError:
                pass

        if i == N - 1:
            raise ValueError(
                "Unable to find any valid load"
                + " fraction with {:d} random samples.".format(N)
            )

        low = tol
        probe = startProbe
        delta_low = low - probe

        new_low = probe + delta_low

        while abs(2 * delta_low) > tol:
            try:
                _, lt_i, lg_i = f(new_low)
                records.append((new_low, lt_i))
                probe = new_low
            except ValueError:
                delta_low *= 0.5
            finally:
                new_low = probe + delta_low

        low = probe

        high = 1 - tol
        probe = startProbe
        delta_high = high - probe

        new_high = probe + delta_high
        while abs(2 * delta_high) > tol and new_high < 1:
            try:
                _, lt_i, lg_i = f(new_high)
                records.append((new_high, lt_i))
                probe = new_high
            except ValueError:
                delta_high *= 0.5
            finally:
                new_high = probe + delta_high

        high = probe

        if abs(high - low) < tol:
            raise ValueError("No range of values satisfying constraint.")

        if len(records) > 2:
            records.sort(key=lambda x: x[0])
            for l, m, h in zip(records[:-2], records[1:-1], records[2:]):
                if l[1] > m[1] and h[1] > m[1]:
                    low = l[0]
                    high = h[0]

        delta = high - low

        # Edge values are some times only semi-stable, i.e. when calling
        # f() with the same value will spuriously raise value errors. Therefore
        # we conservatively shrink the range by tolerance to avoid this issue.

        low += delta * tol
        high -= delta * tol

        """
        Step 2, gss to min.
        """
        lf_low, lf_high = gss(
            lambda lf: f(lf)[1], low, high, x_tol=tol, findMin=True
        )
        lf = 0.5 * (lf_high + lf_low)

        e_1, l_t, l_g = f(lf)

        return lf, e_1, l_g


if __name__ == "__main__":
    from prop import GrainComp, MultPerfGeometry, SimpleGeometry

    compositions = GrainComp.readFile("data/propellants.csv")
    S22 = compositions["ATK PRD(S)22"]
    M8 = compositions["M8"]
    M1 = compositions["M1"]

    S22S = Propellant(S22, SimpleGeometry.TUBE, 1, 2.5)
    M8S = Propellant(M8, SimpleGeometry.TUBE, 1, 2.5)
    M17P = Propellant(M1, MultPerfGeometry.SEVEN_PERF_CYLINDER, 1.0, 2.50)
    test = ConstrainedRecoiless(
        caliber=93e-3,
        shotMass=4,
        propellant=M17P,
        startPressure=10e6,
        dragCoefficient=5e-2,
        designPressure=5e6,
        designVelocity=100,
        nozzleExpansion=4,
        nozzleEfficiency=0.92,
        chambrage=1.5,
    )
    print(
        test.solve(
            loadFraction=0.5,
            chargeMassRatio=0.309 / 4,
            tol=1e-3,
            minWeb=1e-6,
            maxLength=100,
            control=POINT_PEAK_BREECH,
        )
    )
    """
    datas = []
    for i in range(1):
        datas.append(
            test.findMinV(
                chargeMassRatio=0.309 / 4,
                tol=1e-4,
                minWeb=100e-6,
                maxLength=100,
            )
        )

    from tabulate import tabulate

    means = [sum(x) / len(datas) for x in zip(*datas)]

    delta = []
    for line in datas:
        delta.append((v - m) / m for v, m in zip(line, means))

    print(tabulate(datas, headers=("load fract.", "web", "length")))
    print(*means)
    print(
        tabulate(
            delta, headers=("load fract.", "web", "length"), floatfmt=".3e"
        )
    )
    """
