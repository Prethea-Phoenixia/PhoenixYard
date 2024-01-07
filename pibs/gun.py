from math import pi, log, inf, exp
from num import gss, RKF78, cubic, intg, secant, dekker
from prop import GrainComp, Propellant


DOMAIN_TIME = "DOMAIN_TIME"
DOMAIN_LENG = "DOMAIN_LENG"

POINT_START = "SHOT_START"
POINT_PEAK_AVG = "PEAK_AVG_P"
POINT_PEAK_BREECH = "PEAK_BREECH_P"
POINT_PEAK_SHOT = "PEAK_SHOT_P"
POINT_FRACTURE = "FRACTURE"
POINT_BURNOUT = "BURNOUT"
POINT_EXIT = "SHOT_EXIT"

SOL_LAGRANGE = "SOL_LAGRANGE"
SOL_PIDDUCK = "SOL_PIDDUCK"
SOL_MAMONTOV = "SOL_MAMONTOV"

minTol = 1e-14  # based on experience


def pidduck(wpm, k, tol):
    """
    Pidduck's limiting solution to the Lagrange problem.
    wpm : w/(phi_1 * m), charge mass to equivalent corrected (fictitious) shot
          weight
    k   : adiabatic index of the gas, in practice this is not a great influence
    tol : numerical tolerance

    Pidduck's solution is reduced to that of M.A.Mamontov's solution at k -> 1,
    however numerical difficulty necessitates taking the limit.
    """
    if k < 1:
        raise ValueError("Invalid adiabatic index passed", k)

    def f(Omega, x):
        if k == 1:
            return exp(-Omega * x**2)
        else:
            return (1 - Omega * x**2) ** (1 / (k - 1))

    def g(Omega, x):
        return f(Omega, x) * x**2

    def f_Omega(Omega):
        if Omega == 0:
            return -inf

        I, _ = intg(lambda x: f(Omega, x), 0, 1, tol)

        if k == 1:
            return I - 0.5 * wpm * exp(-Omega) / Omega
        else:
            return I - 0.5 * ((k - 1) / k) * wpm * (
                (1 - Omega) ** (k / (k - 1)) / Omega
            )

    a, b = dekker(f_Omega, 0, 1, tol, y_abs_tol=minTol)
    Omega = 0.5 * (a + b)

    if k == 1:
        labda_1 = (exp(Omega) - 1) / wpm
    else:
        labda_1 = ((1 - Omega) ** (k / (1 - k)) - 1) / wpm
    # Pidduck's solution

    I_u, _ = intg(lambda x: g(Omega, x), 0, 1, tol)
    I_l, _ = intg(lambda x: f(Omega, x), 0, 1, tol)
    labda_2 = I_u / I_l

    return labda_1, labda_2


class Gun:
    def __init__(
        self,
        caliber,
        shotMass,
        propellant,
        grainSize,
        chargeMass,
        chamberVolume,
        startPressure,
        lengthGun,
        chambrage,
        structuralMaterial=None,
        structuralSafetyFactor=1.1,
        dragCoefficient=0,
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
                lengthGun <= 0,
                chambrage < 1,
                dragCoefficient < 0,
                dragCoefficient >= 1,
                structuralSafetyFactor <= 1,
            )
        ):
            raise ValueError("Invalid gun parameters")

        if chargeMass > (propellant.maxLF * propellant.rho_p * chamberVolume):
            raise ValueError(
                "Specified Load Fraction Violates Geometrical Constraint"
            )

        self.propellant = propellant
        self.caliber = caliber

        self.e_1 = 0.5 * grainSize
        self.S = (0.5 * caliber) ** 2 * pi
        self.m = shotMass
        self.omega = chargeMass
        self.V_0 = chamberVolume
        self.p_0 = startPressure
        self.l_g = lengthGun
        self.chi_k = chambrage  # ration of l_0 / l_chamber
        self.l_0 = self.V_0 / self.S
        self.l_c = self.l_0 / self.chi_k
        self.Delta = self.omega / self.V_0

        self.phi_1 = 1 / (1 - dragCoefficient)  # drag work coefficient
        self.material = structuralMaterial
        self.ssf = structuralSafetyFactor
        self.is_af = autofrettage

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
                    + " In practice this implies a detonation of the chmaber"
                    + " will likely occur."
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
                "Propellant either could not develop enough pressure to"
                + " overcome start pressure, or has burnt to post fracture."
            )

        self.Z_0 = Zs[0]  # pick the smallest solution

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

    def _f_p_bar(self, Z, l_bar, v_bar):
        psi = self.f_psi_Z(Z)

        l_psi_bar = 1 - self.Delta * (
            (1 - psi) / self.rho_p + (self.alpha * psi)
        )

        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

        return p_bar

    def _ode_t(self, t_bar, Z, l_bar, v_bar, _):
        psi = self.f_psi_Z(Z)

        l_psi_bar = (
            1
            - self.Delta / self.rho_p
            - self.Delta * (self.alpha - 1 / self.rho_p) * psi
        )
        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                0.25 * k * (k + 1) * v_r**2
                + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar
        else:
            p_d_bar = 0

        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n
        else:
            dZ = 0  # dZ/dt_bar

        dl_bar = v_bar
        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar)

        return (dZ, dl_bar, dv_bar)

    def _ode_l(self, l_bar, t_bar, Z, v_bar, _):
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition."""

        psi = self.f_psi_Z(Z)

        l_psi_bar = (
            1
            - self.Delta / self.rho_p
            - self.Delta * (self.alpha - 1 / self.rho_p) * psi
        )
        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                0.25 * k * (k + 1) * v_r**2
                + k * v_r * (1 + (0.25 * (k + 1) * v_r) ** 2) ** 0.5
            ) * self.p_a_bar
        else:
            p_d_bar = 0

        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n / v_bar
        else:
            dZ = 0  # dZ /dl_bar

        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar) / v_bar  # dv_bar/dl_bar
        dt_bar = 1 / v_bar  # dt_bar / dl_bar

        return (dt_bar, dZ, dv_bar)

    def _ode_Z(self, Z, t_bar, l_bar, v_bar, _):
        psi = self.f_psi_Z(Z)

        l_psi_bar = (
            1
            - self.Delta / self.rho_p
            - self.Delta * (self.alpha - 1 / self.rho_p) * psi
        )
        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                0.25 * k * (k + 1) * v_r**2
                + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar

        else:
            p_d_bar = 0

        if Z <= self.Z_b:
            dt_bar = (
                2 * self.B / self.theta
            ) ** 0.5 * p_bar**-self.n  # dt_bar/dZ
            dl_bar = v_bar * dt_bar  # dv_bar/dZ
            dv_bar = 0.5 * self.theta * (p_bar - p_d_bar) * dt_bar
        else:
            # technically speaking it is undefined in this area
            dt_bar = 0  # dt_bar/dZ
            dl_bar = 0  # dl_bar/dZ
            dv_bar = 0  # dv_bar/dZ

        return (dt_bar, dl_bar, dv_bar)

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

    def _dPdZ(self, Z, l_bar, v_bar):
        psi = self.f_psi_Z(Z)
        p_bar = self._f_p_bar(Z, l_bar, v_bar)
        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n
        else:
            return 0

        l_psi_bar = 1 - self.Delta * (
            (1 - psi) / self.rho_p + (self.alpha * psi)
        )
        dp_bar = (
            (
                (1 + p_bar * self.Delta * (self.alpha - 1 / self.rho_p))
                * self.f_sigma_Z(Z)  # dpsi/dZ
                * dZ  # dZ/dt_bar
                - p_bar * v_bar * (1 + self.theta)
            )
            / (l_bar + l_psi_bar)  # dp_bar/dt_bar
            / dZ
        )  # dp_bar/dZ

        return dp_bar

    def integrate(
        self,
        step=10,
        tol=1e-5,
        dom=DOMAIN_TIME,
        sol=SOL_PIDDUCK,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        **_,
    ):
        """
        Runs a full numerical solution for the gun in the specified domain
        sampled evenly at specified number of step, using a scaled numerical
        tolerance as specified.

        tolerance is meant to be interpreted as the maximum relative deviation
        each component is allowed to have, at each step of integration point.

        Through significant trials and errors, it was determined that for this
        particular system, the error due to compounding does not appear to be
        significant, usually on the order of 1e-16 - 1e-14 as compared to much
        larger for component errors.

        Partial analytical approximation to the Lagrangian problem:
        d rho / dx = 0: Lagrange
          d S / dx = 0: Pidduck
          d T / dx = 0: M.A.Mamontov
        All solutions assumes gas velocity increasing linearlly from 0
        at breech face and shot velocity at shot base.

        """
        record = []

        if any((step < 0, tol < 0)):
            raise ValueError("Invalid integration specification")

        if any((ambientP < 0, ambientRho < 0, ambientGamma < 1)):
            raise ValueError("Invalid ambient condition")

        if sol == SOL_LAGRANGE:
            labda_1, labda_2 = 0.5, 1 / 3
        elif sol == SOL_PIDDUCK:
            labda_1, labda_2 = pidduck(
                self.omega / (self.phi_1 * self.m), self.theta + 1, tol
            )
        elif sol == SOL_MAMONTOV:
            labda_1, labda_2 = pidduck(
                self.omega / (self.phi_1 * self.m), 1, tol
            )
        else:
            raise ValueError("Unknown Solution")

        self.labda_1 = labda_1
        self.labda_2 = labda_2
        # labda_2 = 1 / 3

        Labda = self.l_g / self.l_0
        cc = (
            1 - (1 - 1 / self.chi_k) * log(Labda + 1) / Labda
        )  # chambrage correction factor

        self.phi = self.phi_1 + labda_2 * self.omega / self.m * cc  # per ref.

        self.B = (
            self.S**2
            * self.e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_1**2)
            * (self.f * self.Delta) ** (2 * (1 - self.n))
        )

        self.v_j = (
            2 * self.f * self.omega / (self.theta * self.phi * self.m)
        ) ** 0.5

        tScale = self.l_0 / self.v_j
        pScale = self.f * self.Delta

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
            t_bar_err,
            l_bar_err,
            Z_err,
            v_bar_err,
        ):
            p_bar = self._f_p_bar(Z, l_bar, v_bar)
            bar_data.append((tag, t_bar, l_bar, Z, v_bar, p_bar))

            p_bar_err = abs(self._dPdZ(Z, l_bar, v_bar) * Z_err)

            bar_err.append(
                ("L", t_bar_err, l_bar_err, Z_err, v_bar_err, p_bar_err)
            )

        updBarData(
            tag=POINT_START,
            t_bar=0,
            l_bar=0,
            Z=Z_0,
            v_bar=0,
            t_bar_err=0,
            l_bar_err=0,
            Z_err=0,
            v_bar_err=0,
        )

        record.append((0, (0, self.psi_0, 0, p_bar_0)))

        Z_i = Z_0
        Z_j = Z_b
        N = 1
        Delta_Z = Z_b - Z_0

        t_bar_i, l_bar_i, v_bar_i = 0, 0, 0

        isBurnOutContained = True

        """
        Instead of letting the integrator handle the heavy lifting, we
        partition Z and integrate upwards until either barrel exit or
        burnout has been achieved. This seek-and-validate inspired
        from bisection puts the group of value denoted by subscript
        i within or on the muzzle, with the propellant either still
        burning or right on the burnout point..
        """
        ztlv_record = [(Z_0, (0, 0, 0))]
        p_max = 1e9  # 1GPa
        p_bar_max = p_max / pScale

        def abort(x, ys, record):
            Z = x
            t_bar, l_bar, v_bar = ys

            p_bar = self._f_p_bar(Z, l_bar, v_bar)

            return l_bar > l_g_bar or p_bar > p_bar_max or v_bar < 0

        while Z_i < Z_b:  # terminates if burnout is achieved
            ztlv_record_i = []
            if Z_j == Z_i:
                raise ValueError(
                    "Numerical accuracy exhausted in search of exit/burnout point."
                )
            try:
                if Z_j > Z_b:
                    Z_j = Z_b
                Z, (t_bar_j, l_bar_j, v_bar_j), _ = RKF78(
                    self._ode_Z,
                    (t_bar_i, l_bar_i, v_bar_i),
                    Z_i,
                    Z_j,
                    relTol=tol,
                    absTol=tol**2,
                    minTol=minTol,
                    abortFunc=abort,
                    record=ztlv_record_i,
                )

                p_bar_j = self._f_p_bar(Z_j, l_bar_j, v_bar_j)

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
                ztlv_record.extend(ztlv_record_i)
                if v_bar_j <= 0:
                    Z, t_bar, l_bar, v_bar = (
                        ztlv_record[-1][0],
                        *ztlv_record[-1][1],
                    )

                    raise ValueError(
                        "Squib load condition detected: Shot stopped in bore.\n"
                        + "Shot is last calculated at:.0f} mm at {:.0f} mm/s after {:.0f} ms".format(
                            l_bar * self.l_0 * 1e3,
                            v_bar * self.v_j * 1e3,
                            t_bar * tScale * 1e3,
                        )
                    )

                if any(v < 0 for v in (t_bar_j, l_bar_j, p_bar_j)):
                    raise ValueError(
                        "Numerical Integration diverged: negative"
                        + " values encountered in results.\n"
                        + "{:.0f} ms, {:.0f} mm, {:.0f} m/s, {:.0f} MPa".format(
                            t_bar_j * tScale * 1e3,
                            l_bar_j * self.l_0 * 1e3,
                            v_bar_j * self.v_j,
                            p_bar_j * pScale / 1e6,
                        )
                    )

                # if Z != Z_j:
                if p_bar_j > p_bar_max:
                    raise ValueError(
                        "Nobel-Abel EoS is generally accurate enough below"
                        + " 600MPa. However, Unreasonably high pressure "
                        + "(>{:.0f} MPa) was encountered.".format(p_max / 1e6),
                    )  # in practice most press-related spikes are captured here

                (t_bar_i, l_bar_i, v_bar_i) = (t_bar_j, l_bar_j, v_bar_j)

                Z_i = Z_j
                """
                this way the group of values denoted by _i is always updated
                as a group.
                """
                Z_j += Delta_Z / N

        if t_bar_i == 0:
            raise ValueError("burnout point found to be at the origin.")

        """
        Cludge code to force the SoE past the discontinuity at Z = Z_b, since
        we wrote the SoE to be be piecewise continous from (0, Z_b] and (Z_b,
        +inf) it is necessary to do this to prevent the RKF integrator coming
        up with irreducible error estimates and driving the step size to 0
        around Z = Z_b
        """
        if isBurnOutContained:
            Z_i = Z_b + tol

        record.extend(
            (
                t_bar,
                (
                    l_bar,
                    self.f_psi_Z(Z),
                    v_bar,
                    self._f_p_bar(Z, l_bar, v_bar),
                ),
            )
            for (Z, (t_bar, l_bar, v_bar)) in ztlv_record
            if t_bar != 0
        )

        """
        Subscript e indicate exit condition.
        At this point, since its guaranteed that point i will be further
        towards the chamber of the firearm than point e, we do not have
        to worry about the dependency of the correct behaviour of this ODE
        on the positive direction of integration.
        """

        ltzv_record = []
        (
            l_bar,
            (t_bar_e, Z_e, v_bar_e),
            (t_bar_err, Z_err, v_bar_err),
        ) = RKF78(
            self._ode_l,
            (t_bar_i, Z_i, v_bar_i),
            l_bar_i,
            l_g_bar,
            relTol=tol,
            absTol=tol**2,
            minTol=minTol,
            record=ltzv_record,
        )

        if l_bar != l_g_bar:
            if v_bar_e <= 0:
                raise ValueError(
                    "Squib load condition detected post burnout:"
                    + " Round stopped in bore at {:}".format(l_bar * self.l_0)
                )

        record.extend(
            (
                t_bar,
                (
                    l_bar,
                    self.f_psi_Z(Z),
                    v_bar,
                    self._f_p_bar(Z, l_bar, v_bar),
                ),
            )
            for (l_bar, (t_bar, Z, v_bar)) in ltzv_record
        )

        updBarData(
            tag=POINT_EXIT,
            t_bar=t_bar_e,
            l_bar=l_g_bar,
            Z=Z_e,
            v_bar=v_bar_e,
            t_bar_err=t_bar_err,
            l_bar_err=0,
            Z_err=Z_err,
            v_bar_err=v_bar_err,
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
                (t_bar_f, l_bar_f, v_bar_f),
                (t_bar_err_f, l_bar_err_f, v_bar_err_f),
            ) = RKF78(
                self._ode_Z,
                (0, 0, 0),
                Z_0,
                1,
                relTol=tol,
                absTol=tol**2,
                minTol=minTol,
            )

            updBarData(
                tag=POINT_FRACTURE,
                t_bar=t_bar_f,
                l_bar=l_bar_f,
                Z=1,
                v_bar=v_bar_f,
                t_bar_err=t_bar_err_f,
                l_bar_err=l_bar_err_f,
                Z_err=0,
                v_bar_err=v_bar_err_f,
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
                (t_bar_b, l_bar_b, v_bar_b),
                (t_bar_err_b, l_bar_err_b, v_bar_err_b),
            ) = RKF78(
                self._ode_Z,
                (0, 0, 0),
                Z_0,
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
                t_bar_err=t_bar_err_b,
                l_bar_err=l_bar_err_b,
                Z_err=0,
                v_bar_err=v_bar_err_b,
            )

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.

        we hereby simply solve p golden section searching it from origin
        to point e, i.e. inside the barrel.
        """

        def findPeak(f, tag):
            """
            tolerance is specified a bit differently for gold section search
            GSS tol is the length between the upper bound and lower bound
            of the maxima/minima, thus by our definition of tolerance (one
            sided), we take the median value.
            """
            t_bar_tol = tol * min(
                t for t in (t_bar_e, t_bar_b, t_bar_f) if t is not None
            )

            t_bar_1, t_bar_2 = gss(
                f,
                0,
                t_bar_e if t_bar_b is None else t_bar_b,
                x_tol=t_bar_tol,
                findMin=False,
            )

            t_bar = 0.5 * (t_bar_1 + t_bar_2)

            (
                _,
                (Z, l_bar, v_bar),
                (Z_err, l_bar_err, v_bar_err),
            ) = RKF78(
                self._ode_t,
                (Z_0, 0, 0),
                0,
                t_bar,
                relTol=tol,
                absTol=tol**2,
                minTol=minTol,
            )
            t_bar_err = 0.5 * t_bar_tol

            updBarData(
                tag=tag,
                t_bar=t_bar,
                l_bar=l_bar,
                Z=Z,
                v_bar=v_bar,
                t_bar_err=t_bar_err,
                l_bar_err=l_bar_err,
                Z_err=Z_err,
                v_bar_err=v_bar_err,
            )

        def f(t, m="a"):
            Z, l_bar, v_bar = RKF78(
                self._ode_t,
                (Z_0, 0, 0),
                0,
                t,
                relTol=tol,
                absTol=tol**2,
                minTol=minTol,
            )[1]

            p_bar = self._f_p_bar(Z, l_bar, v_bar)

            if m == "a":
                return p_bar

            Ps, Pb = self.toPsPb(l_bar * self.l_0, p_bar * pScale)
            if m == "s":
                return Ps / pScale
            elif m == "b":
                return Pb / pScale

        findPeak(lambda x: f(x, "a"), POINT_PEAK_AVG)
        findPeak(lambda x: f(x, "s"), POINT_PEAK_SHOT)
        findPeak(lambda x: f(x, "b"), POINT_PEAK_BREECH)

        """
        populate data for output purposes
        """
        try:
            if dom == DOMAIN_TIME:
                (Z_j, l_bar_j, v_bar_j, t_bar_j) = (Z_0, 0, 0, 0)
                for j in range(step):
                    t_bar_k = t_bar_e / (step + 1) * (j + 1)
                    (
                        _,
                        (Z_j, l_bar_j, v_bar_j),
                        (Z_err, l_bar_err, v_bar_err),
                    ) = RKF78(
                        self._ode_t,
                        (Z_j, l_bar_j, v_bar_j),
                        t_bar_j,
                        t_bar_k,
                        relTol=tol,
                        absTol=tol**2,
                        minTol=minTol,
                    )
                    t_bar_j = t_bar_k

                    updBarData(
                        tag="",
                        t_bar=t_bar_j,
                        l_bar=l_bar_j,
                        Z=Z_j,
                        v_bar=v_bar_j,
                        t_bar_err=0,
                        l_bar_err=l_bar_err,
                        Z_err=Z_err,
                        v_bar_err=v_bar_err,
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
                Z_j, l_bar_j, v_bar_j = RKF78(
                    self._ode_t,
                    (Z_0, 0, 0),
                    0,
                    t_bar_j,
                    relTol=tol,
                    absTol=tol**2,
                    minTol=minTol,
                )[1]

                for j in range(step):
                    l_bar_k = l_g_bar / (step + 1) * (j + 1)

                    (
                        _,
                        (t_bar_j, Z_j, v_bar_j),
                        (t_bar_err, Z_err, v_bar_err),
                    ) = RKF78(
                        self._ode_l,
                        (t_bar_j, Z_j, v_bar_j),
                        l_bar_j,
                        l_bar_k,
                        relTol=tol,
                        absTol=tol**2,
                        minTol=minTol,
                    )
                    l_bar_j = l_bar_k

                    updBarData(
                        tag="",
                        t_bar=t_bar_j,
                        l_bar=l_bar_j,
                        Z=Z_j,
                        v_bar=v_bar_j,
                        t_bar_err=t_bar_err,
                        l_bar_err=0,
                        Z_err=Z_err,
                        v_bar_err=v_bar_err,
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

        data = []
        error = []

        for bar_dataLine, bar_errorLine in zip(bar_data, bar_err):
            dtag, t_bar, l_bar, Z, v_bar, p_bar = bar_dataLine
            # fmt:off
            (
                etag, t_bar_err, l_bar_err, Z_err, v_bar_err, p_bar_err
            ) = bar_errorLine
            # fmt:on

            t, t_err = (t_bar * tScale, t_bar_err * tScale)
            l, l_err = (l_bar * self.l_0, l_bar_err * self.l_0)
            psi, psi_err = (self.f_psi_Z(Z), abs(self.f_sigma_Z(Z) * Z_err))
            v, v_err = v_bar * self.v_j, v_bar_err * self.v_j

            p, p_err = p_bar * pScale, p_bar_err * pScale

            ps, pb = self.toPsPb(l, p)

            T = self._T(psi, l, p)

            data.append((dtag, t, l, psi, v, pb, p, ps, T))
            error.append(
                (
                    etag,
                    t_err,
                    l_err,
                    psi_err,
                    v_err,
                    "---",
                    p_err,
                    "---",
                    "---",
                )
            )

        """
        scale the records too
        """

        for t_bar, (l_bar, psi, v_bar, p_bar) in record:
            t = t_bar * tScale
            if t in [line[1] for line in data]:
                continue

            l = l_bar * self.l_0
            p = p_bar * pScale
            ps, pb = self.toPsPb(l, p)
            data.append(
                (
                    "*",
                    t_bar * tScale,
                    l,
                    psi,
                    v_bar * self.v_j,
                    pb,
                    p,
                    ps,
                    self._T(psi, l, p),
                ),
            )
            errLine = ("L", *("--" for _ in range(8)))
            error.append(errLine)

        data, error = zip(
            *(
                (a, b)
                for a, b in sorted(zip(data, error), key=lambda x: x[0][1])
            )
        )

        # calculate a pressure and flow velcoity tracing.
        p_trace = []
        l_c = self.l_c

        for line in data:
            tag, t, l, psi, v, p_b, p, p_s, T = line
            p_line = []
            for i in range(step):
                x = i / step * (l + l_c)
                p_x, _ = self.toPxU(l, p_s, p_b, v, x)
                p_line.append((x, p_x))

            p_line.append((l + l_c, p_s))
            p_trace.append((tag, psi, T, p_line))

        try:
            if self.material is None:
                raise ValueError("Structural material not specified")

            structure = self.getStructural(data, step, tol)

        except Exception as e:
            import sys, traceback

            exc_type, exc_value, exc_traceback = sys.exc_info()
            errMsg = "".join(
                traceback.format_exception(exc_type, exc_value, exc_traceback)
            )
            print(errMsg)
            structure = [None, None, None, None]

        return data, error, p_trace, structure

    def getEff(self, vg):
        """
        te: thermal efficiency
        be: ballistic efficiency
        """
        te = (vg / self.v_j) ** 2
        be = te / self.phi
        return te, be

    def toPsPb(self, l, p):
        """
        Convert average chamber pressure at a certain travel to
        shot base pressure, and breech face pressure

        l: travel of the projectile
        p: average pressure

        Ps: pressure at shot
        Pb: pressure at breech
        """
        Labda_g = l / self.l_0
        labda_1_prime = (
            self.labda_1 * (1 / self.chi_k + Labda_g) / (1 + Labda_g)
        )
        labda_2_prime = (
            self.labda_2 * (1 / self.chi_k + Labda_g) / (1 + Labda_g)
        )

        factor_s = 1 + labda_2_prime * (
            self.omega / (self.phi_1 * self.m)
        )  # factor_b = P/P_b = phi / phi_1

        factor_b = (self.phi_1 * self.m + labda_2_prime * self.omega) / (
            self.phi_1 * self.m + labda_1_prime * self.omega
        )

        return p / factor_s, p / factor_b

    def toPxU(self, l, p_s, p_b, v, x):
        """
        Convert the average chamber to pressure and gas flow speed
        at arbitrary point x for a projectile travel of l and average pressure
        of p, **assuming the Lagrangian distribution**.

        Note that with the current state of research, only characteristic point
        values are available for other distributions, use toPsPb() instead for that.

        l: projectile travel
        p_s: pressure of shot
        p_b: pressure of breech
        x: probe point, start from the breech face.
        """
        L_1 = l
        L_0 = self.l_0 / self.chi_k  # physical length of the chamber.

        A_1 = self.S
        A_0 = A_1 * self.chi_k

        # p_0 = p_s * (
        #     1
        #     + self.labda_1  # 0.5
        #     * (self.omega / (self.phi_1 * self.m))  # epsilon_0
        #     * (1 - (self.l_0 / (self.l_0 + l)) ** 2)  # (1- theta_0**2)
        # )

        r = self.chi_k * x if x < L_0 else (x - L_0) + self.l_0
        p_x = p_s * (
            1
            + self.labda_1  # 0.5
            * (self.omega / (self.phi_1 * self.m))  # epsilon_0
            * (1 - (r / (self.l_0 + l)) ** 2)  # (1- theta_0**2))
        )

        if x < L_0:
            u = A_1 * x * v / (self.V_0 + A_1 * L_1)
        else:
            u = (A_1 * x + (A_0 - A_1) * L_0) * v / (self.V_0 + A_1 * L_1)

        return p_x, u

    def getStructural(self, data, step, tol):
        # step 1. calculate the barrel mass
        r = 0.5 * self.caliber
        l_c = self.l_c
        l_g = self.l_g
        chi_k = self.chi_k
        sigma = self.material.Y
        S = self.S

        r_b = r * chi_k**0.5  # radius of breech
        S_b = S * chi_k  # area of breech

        x_probes = (
            [i / step * l_c for i in range(step)]
            + [l_c * (1 - tol)]
            + [i / step * l_g + l_c for i in range(step)]
            + [l_g + l_c]
        )
        p_probes = [0] * len(x_probes)

        for line in data:
            tag, t, l, psi, v, p_b, p, p_s, T = line
            for i, x in enumerate(x_probes):
                if (x - l_c) <= l:
                    p_x, _ = self.toPxU(l, p_s, p_b, v, x)
                    p_probes[i] = max(p_probes[i], p_x)
                else:
                    break

        for i, p in enumerate(p_probes):
            # x = x_probes[i]
            p_probes[i] = p * self.ssf

        if self.is_af:
            """
            m : r_a / r_i
            k : r_o / r_i
            n : p_vM_max / sigma

            1 < m < k

            The point of optimum autofrettage, or the minimum autofrettage
            necessary to contain the working pressure, is
            """
            i = x_probes.index(l_c)
            x_c, p_c = x_probes[:i], p_probes[:i]  # c for chamber
            x_b, p_b = x_probes[i:], p_probes[i:]  # b for barrel

            V_c, rho_c = Gun._Vrho_k(
                x_c,
                p_c,
                [S * chi_k for _ in x_c],
                sigma,
                tol,
            )

            V_b, rho_b = Gun._Vrho_k(
                x_b,
                p_b,
                [S for _ in x_b],
                sigma,
                tol,
                p_ref=max(p_c),
                k_max=rho_c[-1] * chi_k**0.5,
            )

            V = V_c + V_b
            rho_probes = rho_c + rho_b

        else:
            """
            The yield criterion chosen here is the fourth strength
            theory (von Mises) as it is generally accepted to be the most
            applicable for this application.

            The limiting stress points circumferentially along the circum-
            ference of the barrel.

            P_4 = sigma_e * (rho^2-1)/ (3*rho**4 + 1) ** 0.5
            lim (x->inf) (x^2-1)/sqrt(3*x**4+1) = 1/sqrt(3)

            the inverse of (x^2-1)/sqrt(3*x**4+1) is:
            sqrt(
                [-sqrt(-x**2*(3*x**2-4)) - 1]/(3 * x**2 - 1)
            )
            (x < -1 or x > 1)
            and
            sqrt(
                [sqrt(-x**2*(3*x**2-4)) - 1]/(3 * x**2 - 1)
            )
            (x from -1 to 1)
            """
            rho_probes = []
            V = 0

            for p in p_probes:
                y = p / sigma
                if y > 3**-0.5:
                    raise ValueError(
                        f"Limit to conventional construction ({sigma * 3*1e-6:.3f} MPa)"
                        + " exceeded in section."
                    )
                rho = (
                    (1 + y * (4 - 3 * y**2) ** 0.5) / (1 - 3 * y**2)
                ) ** 0.5

                rho_probes.append(rho)

            for i in range(len(x_probes) - 1):
                x_0 = x_probes[i]
                x_1 = x_probes[i + 1]
                rho_0 = rho_probes[i]
                rho_1 = rho_probes[i + 1]
                dV = (rho_1**2 + rho_0**2 - 2) * 0.5 * S * (x_1 - x_0)
                if x_1 <= l_c:
                    V += dV * chi_k
                else:
                    V += dV

        P__sigma = p_probes[0] / sigma
        R2__rb = rho_probes[0]

        R1__R2 = max((1 - 1 / R2__rb * (P__sigma) ** 0.5) ** 0.5, R2__rb**-1)
        R1__rb = R1__R2 * R2__rb
        L__rb = 0.5 * (
            P__sigma**0.5
            * (R1__rb**2 / P__sigma - R1__R2**2 / (1 - R1__R2**2)) ** -0.5
        )

        L = L__rb * r_b * 2  # *2 to account for screw interruption
        R1 = R1__rb * r_b

        bore = []
        for x, rho in zip(x_probes, rho_probes):
            if x < l_c:
                bore.append((x, rho * r_b))
            else:
                bore.append((x, rho * r))

        bore_mass = (
            V + (R2__rb**2 - R1__rb**2) * S_b * L
        ) * self.material.rho

        breech = [(-L, R1), (0.0, R1)]

        breech_mass = R1__rb**2 * S_b * L * self.material.rho

        return bore_mass, bore, breech_mass, breech

    @staticmethod
    def _Vrho_k(
        x_s, p_s, S_s, sigma, tol, k_max=None, k_min=None, index=0, p_ref=None
    ):
        def f(m):
            rho_s = []
            V = 0
            for p in p_s:
                sigma_max = Gun._sigma_vM(m, p, m, sigma)
                """
                the limit as k -> +inf for the stress is:

                lim sigma_tr =
                 k-> +inf
                  sigma * [1 - (1 + 2 ln(m))/m**2 ] + 2p/m**2
                """
                if sigma > sigma_max:
                    rho = m
                else:
                    rho, _ = secant(
                        lambda k: Gun._sigma_vM(k, p, m, sigma),
                        m + tol,
                        m + 2 * tol,
                        y=sigma,
                        x_min=m,
                        y_rel_tol=tol**2,
                    )

                rho_s.append(rho)

            for i in range(len(x_s) - 1):
                x_0, x_1 = x_s[i], x_s[i + 1]
                S_0, S_1 = S_s[i], S_s[i + 1]
                rho_0, rho_1 = rho_s[i], rho_s[i + 1]
                dV = (
                    0.5
                    * ((rho_0**2 - 1) * S_0 + (rho_1**2 - 1) * S_1)
                    * (x_1 - x_0)
                )
                V += dV
            return V, rho_s

        if p_ref is not None:
            p_max = p_ref
        else:
            p_max = max(p_s)

        def sigma_min(m):
            return (
                sigma * (1 - (1 + 2 * log(m)) / m**2) + 2 * p_max / m**2
            ) * (3**0.5 * 0.5)

        m_opt = exp(max(p_s) / sigma * 3**0.5 * 0.5)
        # optimal autofrettage for this pressure

        if sigma_min(m_opt) > sigma:
            """if the minimum junction stress at the optimal autofrettage
            fraction cannot be achieved down to material yield even as
            the thickness goes to infinity, raise an error and abort
            calculation"""
            raise ValueError(
                "Plastic-elastic junction stress exceeds material "
                + f"yield ({sigma * 1e-6:.3f} MPa) for autofrettaged construction."
            )

        elif sigma_min(1) > sigma:
            """if the minimum junction stress at an autofrettage fraction
            of 1 exceeds material yield, implies a certain amount of
            autofrettaging is required to contain the pressure"""
            m_min, _ = dekker(
                sigma_min,
                1,
                m_opt,
                y=sigma,
                y_rel_tol=tol,
            )

            """safety, fudge code to ensure a valid minimum autofrettage
            fraction is found.
            """
            while sigma_min(m_min) > sigma:
                m_min *= 1 + tol**2

        else:
            m_min = 1 + tol**2

        m_max = m_opt
        if k_max is not None:
            """another constraint is the continuation criteria, i.e. the
            barrel base should not require a thickness jump from the front of
            breech, within reason.

            A maximum ratio corresponds to a minimum autofrettage ratio.
            manually setting f(m_opt) to -1 since this is always true,
            but the limitations to numerical accuracy can cause the result
            to float around ~ +/- 10 * tol

            Since in all use cases k_max is set using a seciton of the chamber
            with higher pressure ratings, it is almost guaranteed that the
            corresponding m_min must exist. In the case this is not true, the
            resulting error raised will cause the program to gracefully default
            to finding the minimum autofrettaged mass.
            """
            try:
                m_k, _ = dekker(
                    lambda m: f(m)[1][index] if m != m_opt else -1,
                    m_min,
                    m_opt,
                    y=k_max,
                    y_rel_tol=tol,
                )
                m_min = max(m_k, m_min)
            except ValueError:
                pass

        if k_min is not None:
            try:
                m_k, _ = dekker(
                    lambda m: f(m)[1][index] if m != m_opt else -1,
                    m_min,
                    m_opt,
                    y=k_min,
                    y_rel_tol=tol,
                )
                m_max = min(m_k, m_max)
            except ValueError:
                pass

        # print("m_min", m_min, "m_max", m_max)

        m_best, _ = gss(
            lambda m: f(m)[0], m_min, m_max, y_rel_tol=tol, findMin=True
        )

        # print("m_best", m_best)

        return f(m_best)

    @staticmethod
    def _sigma_vM(k, p, m, sigma):
        """
        k   : probing point, radius ratio
        p   : working pressure (from within)
        m   : autofrettaged (plastically yielded region) rim radius over
            : barrel radius.

        Calculate the von Misses stress at point radius ratio k.
        When supplied with k = m, the result is for at plastic-elastic
        juncture. This is the limiting stress point for an auto-
        frettaged gun barrel under internal pressure loading.

        In general, for increasing m up until m=k, the ability of a barrel to
        tolerate stress increase.
        """
        sigma_tr = (
            sigma
            * (k / m) ** 2
            * ((m / k) ** 2 - (1 - (m / k) ** 2 + 2 * log(m)) / (k**2 - 1))
            + 2 * p / (k**2 - 1) * (k / m) ** 2
        )
        return (
            sigma_tr * 3**0.5 * 0.5
        )  # convert Tresca to von Misses equivalent stress


if __name__ == "__main__":
    """standard 7 hole cylinder has d_0=e_1, hole dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """
    from tabulate import tabulate

    compositions = GrainComp.readFile("data/propellants.csv")

    M17 = compositions["M17"]

    from prop import SimpleGeometry, MultPerfGeometry
    from material import Material

    M17SHC = Propellant(M17, MultPerfGeometry.SEVEN_PERF_CYLINDER, 2, 2.5)

    lf = 0.5
    print("DELTA/rho:", lf)
    cm = 0.1
    test = Gun(
        caliber=0.05,
        shotMass=1.0,
        propellant=M17SHC,
        grainSize=6.66e-3,
        chargeMass=cm,
        chamberVolume=cm / M17SHC.rho_p / lf,
        startPressure=30e6,
        lengthGun=3.5,
        chambrage=1.5,
        dragCoefficient=0.05,
        # autofrettage=False,
    )
    """

    print("\nnumerical: time")
    print(
        tabulate(
            test.integrate(100, 1e-6, dom=DOMAIN_TIME)[0],
            headers=("tag", "t", "l", "phi", "v", "p", "T", "eta"),
        )
    )
    print("\nnumerical: length")
    print(
        tabulate(
            test.integrate(100, 1e-6, dom=DOMAIN_LENG, sol=SOL_MAMONTOV)[0],
            headers=("tag", "t", "l", "phi", "v", "p", "T", "eta"),
        )
    )
    """
    result = test.integrate(100, 1e-3, dom=DOMAIN_TIME, sol=SOL_LAGRANGE)
    print(
        tabulate(
            result[0],
            headers=("tag", "t", "l", "phi", "v", "pb", "p", "ps", "T", "eta"),
        )
    )

    # print(result[2])
    """
    # density:  lbs/in^3 -> kg/m^3, multiply by 27680
    # covolume: in^3/lbs -> m^3/kg, divide by 27680
    # force:    ft-lbs per lb ->J/kg multiply by 2.98907
    # burn rate coefficient:
    # mm/s/(MPa)**exponent -> m/s/Pa**exponent * 1e-3 /(1e6)**exponent
    # cm/s/(MPa)**exponent -> m/s/Pa**exponent * 1e-2 /(1e6)**exponent
    """
