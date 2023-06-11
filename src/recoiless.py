from math import pi, log, exp, inf
from num import gss, RKF78, cubic
from prop import GrainComp, Propellant

from gun import DOMAIN_TIME, DOMAIN_LENG
from gun import (
    POINT_START,
    POINT_PEAK,
    POINT_FRACTURE,
    POINT_BURNOUT,
    POINT_EXIT,
)


class Recoiless:
    def __init__(
        self,
        caliber,
        shotMass,
        propellant,
        grainSize,
        chargeMass,
        chamberVolume,
        startopenPressure,
        lengthGun,
        nozzleExpansion,
        dragCoe=0,
        nozzleEff=0.92,
    ):
        if any(
            (
                caliber <= 0,
                shotMass <= 0,
                chargeMass <= 0,
                grainSize <= 0,
                chamberVolume <= 0,
                lengthGun <= 0,
                nozzleExpansion < 1,
                nozzleEff > 1,
                dragCoe < 0,
                startopenPressure < 0,
            )
        ):
            raise ValueError("Invalid gun parameters")

        e_1 = 0.5 * grainSize
        self.chi_k = 1
        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.propellant = propellant
        self.omega = chargeMass
        self.V_0 = chamberVolume
        self.p_0 = startopenPressure
        self.l_g = lengthGun
        self.chi_0 = nozzleEff
        self.Delta = self.omega / self.V_0
        self.l_0 = self.V_0 / self.S
        Labda = self.l_g / self.l_0
        self.phi_1 = 1 + dragCoe  # drag work coefficient
        self.phi = self.phi_1 + self.omega / (3 * self.m)
        # chamberage effect is explicitly not accounted for.
        self.v_j = (
            2 * self.f * self.omega / (self.theta * self.phi * self.m)
        ) ** 0.5

        if self.p_0 == 0:
            raise NotImplementedError(
                "Current implementation use exponential burn rate and does not"
                + " allow for solving the case with 0 shot start pressure."
            )
        else:
            self.psi_0 = (1 / self.Delta - 1 / self.rho_p) / (
                self.f / self.p_0 + self.alpha - 1 / self.rho_p
            )
            if self.psi_0 <= 0:
                raise ValueError(
                    "Initial burnup fraction is solved to be negative."
                    + " In practice this implies a detonation of the gun breech"
                    + " will likely occur."
                    + " Suggest reducing load fraction."
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
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to overcome"
                + " start pressure, or has burnt to post fracture."
            )

        self.Z_0 = Zs[0]

        # additional calculation for recoiless weapons:
        gamma = self.theta + 1
        self.S_j_bar = self.getCf(gamma, 1) / (
            self.getCf(gamma, nozzleExpansion) * nozzleEff
        )  # the "recoiless condition."
        self.S_j = self.S_j_bar * self.S

        self.K_0 = (2 / (gamma + 1)) ** (
            (gamma + 1) / (2 * (gamma - 1))
        ) * gamma**0.5

        phi_2 = 1
        self.C_A = (
            (self.theta * self.phi * self.m / (2 * self.omega)) ** 0.5
            * self.K_0
            * phi_2
        )  # flow rate value

    def __getattr__(self, attrName):
        try:
            return getattr(self.propellant, attrName)
        except:
            raise AttributeError("object has no '%s'" % attrName)

    def _fp_bar(self, Z, l_bar, eta, tau):
        psi = self.f_psi_Z(Z)
        l_psi_bar = 1 - self.Delta * (
            (1 - psi) / self.rho_p - self.alpha * (psi - eta)
        )

        p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

        return p_bar

    def _ode_t(self, t_bar, Z, l_bar, v_bar, p_bar, eta, tau):
        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ

        l_psi_bar = 1 - self.Delta * (
            (1 - psi) / self.rho_p - self.alpha * (psi - eta)
        )

        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n
        else:
            dZ = 0  # dZ/dt_bar
        dl_bar = v_bar
        dv_bar = self.theta * 0.5 * p_bar

        deta = self.C_A * self.S_j_bar * p_bar / tau**0.5  # deta / dt_bar
        dtau = (
            (1 - tau) * (dpsi * dZ)
            - 2 * v_bar * dv_bar
            - self.theta * tau * deta
        ) / (
            psi - eta
        )  # dtau/dt_bar

        return (dZ, dl_bar, dv_bar, 0, deta, dtau)

    def _pf_t(self, x, *ys):
        t_bar, Z, l_bar, v_bar, p_bar, eta, tau = x, *ys
        p_bar_prime = self._fp_bar(Z, l_bar, eta, tau)
        return (Z, l_bar, v_bar, p_bar_prime, eta, tau)

    def _ode_l(self, l_bar, t_bar, Z, v_bar, p_bar, eta, tau):
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition.

        in general, d/dl_bar = d/dt_bar * dt_bar/dl_bar

        """

        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ

        l_psi_bar = 1 - self.Delta * (
            (1 - psi) / self.rho_p - self.alpha * (psi - eta)
        )

        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n / v_bar
        else:
            dZ = 0  # dZ /dl_bar
        dv_bar = self.theta * 0.5 * p_bar / v_bar  # dv_bar / dl_bar
        dt_bar = 1 / v_bar  # dt_bar / dl_bar

        dp_bar = tau / (l_bar + l_psi_bar) * (psi - eta) * dt_bar

        deta = (
            self.C_A * self.S_j_bar * p_bar / tau**0.5 * dt_bar
        )  # deta / dl_bar
        dtau = (
            (
                (1 - tau) * (dpsi * dZ / dt_bar)  # dpsi/dt_bar
                - 2 * v_bar * (dv_bar / dt_bar)  # dv_bar/dt_bar
                - self.theta * tau * (deta / dt_bar)
            )
            / (psi - eta)
            * dt_bar
        )  # dtau/dl_bar

        return (dt_bar, dZ, dv_bar, 0, deta, dtau)

    def _pf_l(self, x, *ys):
        l_bar, t_bar, Z, v_bar, p_bar, eta, tau = x, *ys
        p_bar_prime = self._fp_bar(Z, l_bar, eta, tau)
        return (t_bar, Z, v_bar, p_bar_prime, eta, tau)

    def _ode_Z(self, Z, t_bar, l_bar, v_bar, p_bar, eta, tau):
        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ

        l_psi_bar = 1 - self.Delta * (
            (1 - psi) / self.rho_p - self.alpha * (psi - eta)
        )

        if Z <= self.Z_b:
            dt_bar = (2 * self.B / self.theta) ** 0.5 * p_bar**-self.n
            dl_bar = v_bar * (2 * self.B / self.theta) ** 0.5 * p_bar**-self.n
            dv_bar = (self.B * self.theta * 0.5) ** 0.5 * p_bar ** (1 - self.n)

        else:
            # technically speaking it is undefined in this area
            dt_bar = 0  # dt_bar/dZ
            dl_bar = 0  # dl_bar/dZ
            dv_bar = 0  # dv_bar/dZ

        deta = (
            self.C_A * self.S_j_bar * p_bar / tau**0.5 * dt_bar
        )  # deta / dZ
        dtau = (
            (
                (1 - tau) * (dpsi / dt_bar)
                - 2 * v_bar * (dv_bar / dt_bar)
                - self.theta * tau * (deta / dt_bar)
            )
            / (psi - eta)
            * dt_bar
        )

        return (dt_bar, dl_bar, dv_bar, 0, deta, dtau)

    def _pf_Z(self, x, *ys):
        Z, t_bar, l_bar, v_bar, p_bar, eta, tau = x, *ys
        p_bar_prime = self._fp_bar(Z, l_bar, eta, tau)
        return (t_bar, l_bar, v_bar, p_bar_prime, eta, tau)

    def _T(self, psi, l, p):
        """
        given pressure and travel, return temperature
        using the Nobel- Abel EOS
        """
        l_psi = self.l_0 * (
            1
            - self.Delta / self.rho_p
            - self.Delta * (self.alpha * -1 / self.rho_p) * psi
        )

        return self.S * p * (l + l_psi) / (self.omega * psi * self.R)

    def integrate(self, steps=10, tol=1e-5, dom=DOMAIN_TIME, record=None):
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
        minTol = 1e-16  # based on experience

        if any((steps < 0, tol < 0)):
            raise ValueError("Invalid integration specification")

        l_g_bar = self.l_g / self.l_0

        tScale = self.l_0 / self.v_j
        pScale = self.f * self.Delta

        p_bar_0 = self.p_0 / pScale
        Z_b = self.Z_b
        Z_0 = self.Z_0

        bar_data = []
        bar_err = []

        def updBarData(
            tag,
            t_bar,
            l_bar,
            Z,
            v_bar,
            p_bar,
            t_bar_err,
            l_bar_err,
            Z_err,
            v_bar_err,
            p_bar_err,
        ):
            bar_data.append((tag, t_bar, l_bar, Z, v_bar, p_bar))
            """
            Worst case scenario for pressure deviation is when:
            Z higher than actual, l lower than actual, v lower than actual
            more burnt propellant, less volume, lower speed of projectile.
            """
            bar_err.append(
                (
                    "L",
                    t_bar_err,
                    l_bar_err,
                    Z_err,
                    v_bar_err,
                    p_bar_err,
                )
            )

        updBarData(
            tag=POINT_START,
            t_bar=0,
            l_bar=0,
            Z=Z_0,
            v_bar=0,
            p_bar=p_bar_0,
            t_bar_err=0,
            l_bar_err=0,
            Z_err=0,
            v_bar_err=0,
            p_bar_err=0,
        )

        if record is not None:
            record.append((0, (0, self.psi_0, 0, p_bar_0)))

        Z_i = Z_0
        Z_j = Z_b
        N = 1
        Delta_Z = Z_b - Z_0

        t_bar_i, l_bar_i, v_bar_i, p_bar_i, eta_i, tau_i = (
            0,
            0,
            0,
            p_bar_0,
            0,
            1,
        )

        isBurnOutContained = True

        """
        Instead of letting the integrator handle the heavy lifting, we
        partition Z and integrate upwards until either barrel exit or
        burnout has been achieved. This seek-and-validate inspired
        from bisection puts the group of value denoted by subscript
        i within or on the muzzle, with the propellant either still
        burning or right on the burnout point..
        """
        ztlvpet_record = []
        p_max = 1e9  # 1GPa
        p_bar_max = p_max / pScale

        def abort(x, ys, o_ys):
            t_bar, l_bar, v_bar, p_bar, eta, tau = ys
            return l_bar > l_g_bar or p_bar > p_bar_max

        while Z_i < Z_b:  # terminates if burnout is achieved
            ztlvpet_record_i = []
            if Z_j == Z_i:
                raise ValueError(
                    "Numerical accuracy exhausted in search of exit/burnout point."
                )
            try:
                if Z_j > Z_b:
                    Z_j = Z_b
                t_bar_j, l_bar_j, v_bar_j, p_bar_j, eta_j, tau_j = RKF78(
                    self._ode_Z,
                    (t_bar_i, l_bar_i, v_bar_i, p_bar_i, eta_i, tau_i),
                    Z_i,
                    Z_j,
                    relTol=tol,
                    absTol=tol,
                    minTol=minTol,
                    parFunc=self._pf_Z,
                    abortFunc=abort,
                    record=ztlvpet_record_i,
                )[1]

            except ValueError as e:
                raise ValueError(
                    "Unable to integrate due to ill defined system, requiring"
                    + " vanishingly step size."
                )

            if p_bar_j > p_bar_max:
                raise ValueError(
                    "Nobel-Abel EoS is generally accurate enough below 600MPa. However,"
                    + " Unreasonably high pressure (>{:.0f} MPa) was encountered.".format(
                        p_max / 1e6
                    )
                )

            if any(
                (
                    t_bar_i == t_bar_j,
                    l_bar_i == l_bar_j,
                    v_bar_i == v_bar_j,
                    p_bar_i == p_bar_j,
                )
            ):
                raise ValueError(
                    "Numerical integration stalled in search of exit/burnout point."
                )

            if l_bar_j >= l_g_bar:
                if abs(l_bar_i - l_g_bar) / (l_g_bar) > tol or l_bar_i == 0:
                    N *= 2
                    Z_j = Z_i + Delta_Z / N
                else:
                    isBurnOutContained = False
                    break  # l_bar_i is solved to within a tol of l_bar_g

            else:
                ztlvpet_record.extend(ztlvpet_record_i)
                (
                    t_bar_i,
                    l_bar_i,
                    v_bar_i,
                    p_bar_i,
                    eta_i,
                    tau_i,
                ) = (
                    t_bar_j,
                    l_bar_j,
                    v_bar_j,
                    p_bar_j,
                    eta_j,
                    tau_j,
                )
                Z_i = Z_j
                """
                this way the group of values denoted by _i is always updated
                as a group.
                """
                Z_j += Delta_Z / N

        if t_bar_i == 0:
            raise ValueError("exit/burnout point found to be at the origin.")

        """
        Cludge code to force the SoE past the discontinuity at Z = Z_b, since 
        we wrote the SoE to be be piecewise continous from (0, Z_b] and (Z_b, +inf)
        it is necessary to do this to prevent the RKF integrator coming up with 
        irreducible error estimates and driving the step size to 0 around Z = Z_b
        """
        if isBurnOutContained:
            Z_i = Z_b + tol

        if record is not None:
            record.extend(
                (
                    t_bar,
                    (
                        l_bar,
                        self.f_psi_Z(Z),
                        v_bar,
                        p_bar,
                    ),
                )
                for (
                    Z,
                    (t_bar, l_bar, v_bar, p_bar, eta, tau),
                ) in ztlvpet_record
            )

            # print(*record, sep="\n")
        """
        Subscript e indicate exit condition.
        At this point, since its guaranteed that point i will be further
        towards the chamber of the firearm than point e, we do not have
        to worry about the dependency of the correct behaviour of this ODE
        on the positive direction of integration.
        """

        ltzvpet_record = []
        (
            _,
            (t_bar_e, Z_e, v_bar_e, p_bar_e, eta_e, tau_e),
            (t_bar_err, Z_err, v_bar_err, p_bar_err, eta_err, tau_err),
        ) = RKF78(
            self._ode_l,
            (t_bar_i, Z_i, v_bar_i, p_bar_i, eta_i, tau_i),
            l_bar_i,
            l_g_bar,
            relTol=tol,
            absTol=tol,
            minTol=minTol,
            record=ltzvpet_record,
            parFunc=self._pf_l,
        )

        if record is not None:
            record.extend(
                (
                    t_bar,
                    (
                        l_bar,
                        self.f_psi_Z(Z),
                        v_bar,
                        p_bar,
                    ),
                )
                for (
                    l_bar,
                    (t_bar, Z, v_bar, p_bar, eta, tau),
                ) in ltzvpet_record
            )
            # print(*record, sep="\n")
            # record.sort(key=lambda line: line[0])

        updBarData(
            tag=POINT_EXIT,
            t_bar=t_bar_e,
            l_bar=l_g_bar,
            Z=Z_e,
            v_bar=v_bar_e,
            p_bar=p_bar_e,
            t_bar_err=t_bar_err,
            l_bar_err=0,
            Z_err=Z_err,
            v_bar_err=v_bar_err,
            p_bar_err=p_bar_err,
        )

        t_bar_f = None
        if Z_b > 1.0 and Z_e >= 1.0:  # fracture point exist and is contained
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile
            movement to charge fracture
            """
            (
                _,
                (t_bar_f, l_bar_f, v_bar_f, p_bar_f, eta_f, tau_f),
                (
                    t_bar_err_f,
                    l_bar_err_f,
                    v_bar_err_f,
                    p_bar_err_f,
                    eta_err_f,
                    tau_err_f,
                ),
            ) = RKF78(
                self._ode_Z,
                (0, 0, 0, p_bar_0, 0, 1),
                Z_0,
                1,
                relTol=tol,
                absTol=tol,
                minTol=minTol,
                parFunc=self._pf_Z,
            )

            updBarData(
                tag=POINT_FRACTURE,
                t_bar=t_bar_f,
                l_bar=l_bar_f,
                Z=1,
                v_bar=v_bar_f,
                p_bar=p_bar_f,
                t_bar_err=t_bar_err_f,
                l_bar_err=l_bar_err_f,
                Z_err=0,
                v_bar_err=v_bar_err_f,
                p_bar_err=p_bar_err_f,
            )

        t_bar_b = None
        if isBurnOutContained:
            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to Z_b, from onset of projectile
            movement to charge burnout.
            """

            (
                _,
                (t_bar_b, l_bar_b, v_bar_b, p_bar_b, eta_b, tau_b),
                (
                    t_bar_err_b,
                    l_bar_err_b,
                    v_bar_err_b,
                    p_bar_err_b,
                    eta_err_b,
                    tau_err_b,
                ),
            ) = RKF78(
                self._ode_Z,
                (0, 0, 0, p_bar_0, 0, 1),
                Z_0,
                Z_b,
                relTol=tol,
                absTol=tol,
                minTol=minTol,
                parFunc=self._pf_Z,
            )

            updBarData(
                tag=POINT_BURNOUT,
                t_bar=t_bar_b,
                l_bar=l_bar_b,
                Z=Z_b,
                v_bar=v_bar_b,
                p_bar=p_bar_b,
                t_bar_err=t_bar_err_b,
                l_bar_err=l_bar_err_b,
                Z_err=0,
                v_bar_err=v_bar_err_b,
                p_bar_err=p_bar_err_b,
            )

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.

        we hereby simply solve p golden section searching it from origin
        to point e, i.e. inside the barrel.
        """

        def f(t):
            _, _, _, p_bar, _, _ = RKF78(
                self._ode_t,
                (Z_0, 0, 0, p_bar_0, 0, 1),
                0,
                t,
                relTol=tol,
                absTol=tol,
                minTol=minTol,
                parFunc=self._pf_t,
            )[1]
            return p_bar

        """
            tolerance is specified a bit differently for gold section search
            GSS tol is the length between the upper bound and lower bound
            of the maxima/minima, thus by our definition of tolerance (one sided),
            we take the median value.
        """

        t_bar_tol = tol * min(
            t for t in (t_bar_e, t_bar_b, t_bar_f) if t is not None
        )

        t_bar_p_1, t_bar_p_2 = gss(
            f,
            0,
            t_bar_e if t_bar_b is None else t_bar_b,
            tol=t_bar_tol,
            findMin=False,
        )

        t_bar_p = 0.5 * (t_bar_p_1 + t_bar_p_2)

        if t_bar_f is not None and abs(t_bar_p - t_bar_f) / t_bar_p < tol:
            # peak pressure occurs sufficiently close to fracture.
            Z_p = 1
            l_bar_p = l_bar_f
            v_bar_p = v_bar_f
            t_bar_p = t_bar_f
            p_bar_p = p_bar_f

            eta_p = eta_f
            tau_p = tau_f

            Z_err_p = 0
            l_bar_err_p = l_bar_err_f
            v_bar_err_p = v_bar_err_f
            t_bar_err_p = t_bar_err_f
            p_bar_err_p = p_bar_err_f

            eta_err_p = eta_err_f
            tau_err_p = tau_err_f

        elif t_bar_b is not None and abs(t_bar_p - t_bar_b) / t_bar_p < tol:
            # peak pressure occurs sufficiently close to burnout.
            Z_p = Z_b
            l_bar_p = l_bar_b
            v_bar_p = v_bar_b
            t_bar_p = t_bar_b
            p_bar_p = p_bar_b

            eta_p = eta_b
            tau_p = tau_b

            Z_err_p = 0
            l_bar_err_p = l_bar_err_b
            v_bar_err_p = v_bar_err_b
            t_bar_err_p = t_bar_err_b
            p_bar_err_p = p_bar_err_b

            eta_err_p = eta_err_b
            tau_err_p = tau_err_b

        else:
            (
                _,
                (Z_p, l_bar_p, v_bar_p, p_bar_p, eta_p, tau_p),
                (
                    Z_err_p,
                    l_bar_err_p,
                    v_bar_err_p,
                    p_bar_err_p,
                    eta_err_p,
                    tau_err_p,
                ),
            ) = RKF78(
                self._ode_t,
                (Z_0, 0, 0, p_bar_0, 0, 1),
                0,
                t_bar_p,
                relTol=tol,
                absTol=tol,
                minTol=minTol,
                parFunc=self._pf_t,
            )
            t_bar_err_p = 0.5 * t_bar_tol

        updBarData(
            tag=POINT_PEAK,
            t_bar=t_bar_p,
            l_bar=l_bar_p,
            Z=Z_p,
            v_bar=v_bar_p,
            p_bar=p_bar_p,
            t_bar_err=t_bar_err_p,
            l_bar_err=l_bar_err_p,
            Z_err=Z_err_p,
            v_bar_err=v_bar_err_p,
            p_bar_err=p_bar_err_p,
        )

        """
        populate data for output purposes
        """
        try:
            if dom == DOMAIN_TIME:
                (Z_j, l_bar_j, v_bar_j, t_bar_j, p_bar_j, eta_j, tau_j) = (
                    Z_0,
                    0,
                    0,
                    0,
                    p_bar_0,
                    0,
                    1,
                )
                for j in range(steps):
                    t_bar_k = t_bar_e / (steps + 1) * (j + 1)
                    (
                        _,
                        (Z_j, l_bar_j, v_bar_j, p_bar_j, eta_j, tau_j),
                        (
                            Z_err,
                            l_bar_err,
                            v_bar_err,
                            p_bar_err,
                            eta_err,
                            tau_err,
                        ),
                    ) = RKF78(
                        self._ode_t,
                        (Z_j, l_bar_j, v_bar_j, p_bar_j, eta_j, tau_j),
                        t_bar_j,
                        t_bar_k,
                        relTol=tol,
                        absTol=tol,
                        minTol=minTol,
                        parFunc=self._pf_t,
                    )
                    t_bar_j = t_bar_k

                    updBarData(
                        tag="",
                        t_bar=t_bar_j,
                        l_bar=l_bar_j,
                        Z=Z_j,
                        v_bar=v_bar_j,
                        p_bar=p_bar_j,
                        t_bar_err=0,
                        l_bar_err=l_bar_err,
                        Z_err=Z_err,
                        v_bar_err=v_bar_err,
                        p_bar_err=p_bar_err,
                    )

            else:
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
                Z_j, l_bar_j, v_bar_j, p_bar_j, eta_j, tau_j = RKF78(
                    self._ode_t,
                    (Z_0, 0, 0, p_bar_0, 0, 1),
                    0,
                    t_bar_j,
                    relTol=tol,
                    absTol=tol,
                    minTol=minTol,
                    parFunc=self._pf_t,
                )[1]

                for j in range(steps):
                    l_bar_k = l_g_bar / (steps + 1) * (j + 1)

                    (
                        _,
                        (t_bar_j, Z_j, v_bar_j, p_bar_j, eta_j, tau_j),
                        (
                            t_bar_err,
                            Z_err,
                            v_bar_err,
                            p_bar_err,
                            eta_err,
                            tau_err,
                        ),
                    ) = RKF78(
                        self._ode_l,
                        (t_bar_j, Z_j, v_bar_j, p_bar_j, eta_j, tau_j),
                        l_bar_j,
                        l_bar_k,
                        relTol=tol,
                        absTol=tol,
                        minTol=minTol,
                        parFunc=self._pf_l,
                    )

                    l_bar_j = l_bar_k

                    updBarData(
                        tag="",
                        t_bar=t_bar_j,
                        l_bar=l_bar_j,
                        Z=Z_j,
                        v_bar=v_bar_j,
                        p_bar=p_bar_j,
                        t_bar_err=t_bar_err,
                        l_bar_err=0,
                        Z_err=Z_err,
                        v_bar_err=v_bar_err,
                        p_bar_err=p_bar_err,
                    )

        except ValueError as e:
            print(e)
        finally:
            pass

        """
        sort the data points
        """

        bar_data, bar_err = zip(
            *(
                (a, b)
                for a, b in sorted(
                    zip(bar_data, bar_err), key=lambda x: x[0][1]
                )
            )
        )

        data = []
        error = []

        for bar_dataLine, bar_errorLine in zip(bar_data, bar_err):
            dtag, t_bar, l_bar, Z, v_bar, p_bar = bar_dataLine
            (
                etag,
                t_bar_err,
                l_bar_err,
                Z_err,
                v_bar_err,
                p_bar_err,
            ) = bar_errorLine

            t, t_err = (t_bar * tScale, t_bar_err * tScale)
            l, l_err = (l_bar * self.l_0, l_bar_err * self.l_0)
            psi, psi_err = (self.f_psi_Z(Z), abs(self.f_sigma_Z(Z) * Z_err))
            v, v_err = v_bar * self.v_j, v_bar_err * self.v_j

            p, p_err = p_bar * pScale, p_bar_err * pScale

            T = self._T(psi, l, p)
            """
            Since for known propellants 1/alpha, or the reciprocal of 
            the covolume is smaller than the density, therefore the free
            volume decrease monotically, the number of molecules increase
            monotically.

            Therefore, the worst case error for temperature is achieved by: 
            psi lower than actual, length longer than actual, pressure higher than actual
            According to the Nobel-Abel EoS
            """
            T_err = max(
                self._T(psi - psi_err, l + l_err, p + p_err) - T,
                T - self._T(psi + psi_err, l - l_err, p - p_err),
            )
            data.append((dtag, t, l, psi, v, p, T))
            error.append((etag, t_err, l_err, psi_err, v_err, p_err, T_err))

        """
        scale the records too
        """
        if record is not None:
            for i, (t_bar, (l_bar, psi, v_bar, p_bar)) in enumerate(record):
                record[i] = (
                    t_bar * tScale,
                    (
                        l_bar * self.l_0,
                        psi,
                        v_bar * self.v_j,
                        p_bar * pScale,
                    ),
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

    def toPb(self, p, L):
        """
        Convert average chamber pressure at a certain travel to the shot bass
        pressure.
        """
        theta_0 = self.V_0 / (self.V_0 + self.S * L)
        epsilon_prime = self.omega / (self.phi_1 * self.m)
        factor = 1 + epsilon_prime / 3 * (
            1 - 1.5 * theta_0**3 * (1 - self.chi_k**-2)
        )

        return p / factor

    def toPt(self, p, L):
        """
        Convert the average pressure at a certain shot travel in the gun to the
        chamber base pressure

        A_0: breech face area,      L_0: chamber size
        A_1: barrel cross section   L_1: shot travel

        A_0 = A_1 * chi_k <- chamber expansion factor
        """

        theta_0 = self.V_0 / (self.V_0 + self.S * L)
        epsilon_prime = self.omega / (self.phi_1 * self.m)
        factor = 1 - epsilon_prime / 6 * (
            1 - 3 * theta_0**2 * (1 - theta_0) * (1 - self.chi_k**-2)
        )  # p/p_t

        return p / factor

    @staticmethod
    def getCf(gamma, Sr, tol=1e-5):
        """
        takes the adiabatic index and area ration between constriction throat
        and the exit, calculate the thrust factor Cf
        See Hunt (1953) appendix I.A01-A03
        Sr = Se/St
        Vr = V/Vt
        """
        Vr_old = 0
        Vr = 1
        while abs(Vr - Vr_old) / Vr > tol:
            Vr_old = Vr
            Vr = (
                (gamma + 1)
                / (gamma - 1)
                * (1 - 2 / (gamma + 1) * (Vr * Sr) ** (1 - gamma))
            ) ** 0.5

        Cf = (2 / (gamma + 1)) ** (gamma / (gamma - 1)) * (
            gamma * Vr + Sr ** (1 - gamma) * Vr ** (-gamma)
        )

        return Cf


if __name__ == "__main__":
    """standard 7 hole cylinder has d_0=e_1, hole dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """
    from tabulate import tabulate

    compositions = GrainComp.readFile("data/propellants.csv")

    M17 = compositions["M17"]
    from prop import SimpleGeometry

    M17SHC = Propellant(M17, SimpleGeometry.SPHERE, 2, 2.5)

    lf = 0.5
    print("DELTA:", lf * M17SHC.maxLF)
    test = Recoiless(
        caliber=0.082,
        shotMass=1.0,
        propellant=M17SHC,
        grainSize=1e-3,
        chargeMass=1,
        chamberVolume=1.0 / M17SHC.rho_p / M17SHC.maxLF / lf,
        startopenPressure=30e6,
        lengthGun=3.5,
        nozzleExpansion=2,
    )
    record = []

    print("\nnumerical: time")
    print(
        tabulate(
            test.integrate(100, 1e-3, dom=DOMAIN_TIME)[0],
            headers=("tag", "t", "l", "phi", "v", "p", "T"),
        )
    )
    print("\nnumerical: length")
    print(
        tabulate(
            test.integrate(100, 1e-3, dom="length")[0],
            headers=("tag", "t", "l", "phi", "v", "p", "T"),
        )
    )
    print(test.getEff(942))
