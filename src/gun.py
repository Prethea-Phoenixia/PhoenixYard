import enum

import csv

from math import pi, log, exp
from num import *


class Geometry(enum.Enum):
    """table 1-4 from ref[1] page 33"""

    def __new__(cls, *args, **kwds):
        value = len(cls.__members__) + 1
        obj = object.__new__(cls)
        obj._value_ = value
        return obj

    def __init__(self, desc, A, B, C, rhoDiv, nHole):
        self.desc = desc
        self.A = A
        self.B = B
        self.C = C

        self.rhoDiv = rhoDiv  # rho divided by (e_1+d_0/2)
        self.nHole = nHole

    SEVEN_PERF_CYLINDER = (
        "7 Perf Cylinder",
        1,
        7,
        0,
        0.2956,
        7,
    )
    SEVEN_PERF_ROSETTE = (
        "7 Perf Rosette",
        2,
        8,
        12 * 3**0.5 / pi,
        0.1547,
        7,
    )
    NINETEEN_PERF_ROSETTE = (
        "19 Perf Rosette",
        3,
        21,
        36 * 3**0.5 / pi,
        0.1547,
        19,
    )
    NINETEEN_PERF_CYLINDER = (
        "19 Perf Cylinder",
        1,
        19,
        0,
        0.3559,
        19,
    )
    NINETEEN_PERF_HEXAGON = (
        "19 Perf Hexagon",
        18 / pi,
        19,
        18 * (3 * 3**0.5 - 1) / pi,
        0.1864,
        19,
    )
    NINETEEN_PERF_ROUNDED_HEXAGON = (
        "19 Perf Rounded Hexagon",
        3**0.5 + 12 / pi,
        19,
        3 - 3**0.5 + 12 * (4 * 3**0.5 - 1) / pi,
        0.1977,
        19,
    )


class GrainComp:
    """
    detonation velocity and pressure exponent are related to:

    r_dot = u_1 * p^n

    tested velocity for n=1 roughly equals
    Nitrocellulose, 6e-7 - 9e-7 m/(MPa*s)
    Nitroglycerin, 7e-7 - 9e-7 m/(MPa*s)
    """

    def __init__(
        self,
        name,
        desc,
        propellantForce,
        covolume,
        density,
        redAdbIndex,
        detVel,
        pressureExp,
        linBurnCoe,
        flameTemp,
    ):
        self.name = name
        self.desc = desc
        self.f = propellantForce
        """
        Propellant force is related to the flame temperature
        by:
            f = R * T_1
            R = R_0/M
            R_0 = 8.314... J/(K*mol)
        where T_1 is the temeprature propellants develop when
        burning in an iso-volume chamber, ignoring losses.
        M is the molar mass of the gas developed. (kg/mol)
        """
        self.alpha = covolume
        self.rho_p = density
        self.theta = redAdbIndex
        self.u_1 = detVel
        self.n = pressureExp
        self.u_0 = linBurnCoe
        self.T_1 = flameTemp
        self.R = self.f / self.T_1

    def readFile(fileName):
        composition = []
        with open(fileName, newline="") as csvfile:
            propReader = csv.reader(
                csvfile,
                delimiter=",",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )

            skipFirstLine = True

            for prop in propReader:
                if skipFirstLine:
                    skipFirstLine = False
                    continue
                (
                    name,
                    desc,
                    adb,
                    density,
                    propellantForce,
                    covolume,
                    pressureExp,
                    burnRateCoe,
                    linBurnCoe,
                    flameTemp,
                ) = prop

                redAdbIndex = float(adb) - 1

                newComp = GrainComp(
                    name,
                    desc,
                    float(propellantForce),
                    float(covolume),
                    float(density),
                    float(redAdbIndex),
                    float(burnRateCoe),
                    float(pressureExp),
                    float(linBurnCoe),
                    float(flameTemp),
                )

                composition.append(newComp)

        return {i.name: i for i in composition}


class Propellant:
    # assumed to be multi-holed propellants
    def __init__(self, composition, propGeom, arcThick, grainPR, grainLDR):
        """
        LDR: length ot diameter ratio
        PR: perf diameter to arc thickness ratio.
        """
        self.composition = composition
        self.geometry = propGeom

        self.e_1 = arcThick / 2
        self.d_0 = arcThick * grainPR

        if propGeom == Geometry.SEVEN_PERF_CYLINDER:
            self.D_0 = (
                8 * self.e_1 + 3 * self.d_0
            )  # enforce geometry constraints
            b, a = self.D_0, 0

        elif propGeom == Geometry.SEVEN_PERF_ROSETTE:
            self.D_0 = 8 * self.e_1 + 3 * self.d_0
            b, a = self.d_0 + 4 * self.e_1, self.d_0 + 2 * self.e_1

        elif propGeom == Geometry.NINETEEN_PERF_ROSETTE:
            self.D_0 = 12 * self.e_1 + 5 * self.d_0
            b, a = self.d_0 + 4 * self.e_1, self.d_0 + 2 * self.e_1

        elif propGeom == Geometry.NINETEEN_PERF_CYLINDER:
            self.D_0 = 12 * self.e_1 + 5 * self.d_0
            b, a = self.D_0, 0

        elif propGeom == Geometry.NINETEEN_PERF_HEXAGON:
            self.D_0 = 12 * self.e_1 + 5 * self.d_0
            b, a = self.d_0 + 2 * self.e_1, self.d_0 + 2 * self.e_1

        elif propGeom == Geometry.NINETEEN_PERF_ROUNDED_HEXAGON:
            self.D_0 = 12 * self.e_1 + 5 * self.d_0
            b, a = self.d_0 + 2 * self.e_1, self.d_0 + 2 * self.e_1

        else:
            raise ValueError(
                "unhandled propellant geometry {}".format(propGeom)
            )

        grainLength = grainLDR * self.D_0
        # grainLength = self.D_0 * grainLDR
        self.c = grainLength / 2

        A, B, C = propGeom.A, propGeom.B, propGeom.C
        self.rho = propGeom.rhoDiv * (self.e_1 + self.d_0 / 2)

        Pi_1 = (A * b + B * self.d_0) / (2 * self.c)
        Q_1 = (C * a**2 + A * b**2 - B * self.d_0**2) / (2 * self.c) ** 2

        S_T = Q_1 * pi * self.c**2
        n = propGeom.nHole

        """
        condition for progressive burning:
        (n-1) > (D_0 + n*d_0)/c
        """

        self.maxLF = S_T / (S_T + n * pi * 0.25 * self.d_0**2)

        """
        maximum load factor, the volume fraction of propellant
        within the geometry of the grain defined. In reality
        this will be a lot lower for any real weapons since
        this both cause undesirably high pressure spike at start
        of shot, and in addition is physically impossible given
        realistic grain packing behaviours
        """

        beta = self.e_1 / self.c

        # first phase of burning rate parameters, Z from 0 to 1
        self.chi = (Q_1 + 2 * Pi_1) / Q_1 * beta
        self.labda = (
            (n - 1 - 2 * Pi_1) / (Q_1 + 2 * Pi_1) * beta
        )  # deliberate misspell to prevent issues with python lambda keyword
        self.mu = -(n - 1) / (Q_1 + 2 * Pi_1) * beta**2

        self.Z_b = (
            self.e_1 + self.rho
        ) / self.e_1  # second phase burning Z upper limit

        psi_s = self.chi * (1 + self.labda + self.mu)

        # second phase of burning rate parameters, Z from 1 to Z_b
        self.chi_s = (1 - psi_s * self.Z_b**2) / (self.Z_b - self.Z_b**2)
        self.labda_s = psi_s / self.chi_s - 1

    def __getattr__(self, attrName):
        try:
            return getattr(self.composition, attrName)
        except:
            raise AttributeError("object has no attribute '%s'" % attrName)


class Gun:
    def __init__(
        self,
        caliber,
        shotMass,
        propellant,
        chargeMass,
        chamberVolume,
        startPressure,
        lengthGun,
        chamberExpansion,
    ):
        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.propellant = propellant
        self.omega = chargeMass
        self.V_0 = chamberVolume
        self.p_0 = startPressure
        self.l_g = lengthGun
        self.chi_k = chamberExpansion  # ration of l_0 / l_chamber
        self.Delta = self.omega / self.V_0
        self.l_0 = self.V_0 / self.S
        Labda = self.l_g / self.l_0
        self.phi = 1 + self.omega / (3 * self.m) * (
            1 - (1 - 1 / self.chi_k) * 2.303 * log(Labda + 1) / Labda
        )  # extra work factor, chamberage effect averaged over entire length
        self.v_j = (
            2 * self.f * self.omega / (self.theta * self.phi * self.m)
        ) ** 0.5

        if self.p_0 == 0:
            raise ValueError(
                "Currently solving for the case with 0 starting pressure"
                + " is not implemented."
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
        # this will overwrite the definition of Geometry.B
        self.B = (
            self.S**2
            * self.e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_1**2)
            * (self.f * self.Delta) ** (2 * (1 - self.n))
        )

        Zs = cubic(
            self.chi * self.mu, self.chi * self.labda, self.chi, -self.psi_0
        )
        # pick a valid solution between 0 and 1
        Zs = tuple(Z for Z in Zs if (Z > 0 and Z < 1))
        if len(Zs) < 1:
            raise ValueError(
                "Propellant has burnt to fracture before start"
                + " of shot movement. Suggest reducing starting"
                + " pressure or specifying higher load fraction."
            )

        self.Z_0 = Zs[0]

    def _fpsi(self, Z):
        if Z < 0:
            return 0
        elif Z < 1.0:
            return self.chi * Z * (1 + self.labda * Z + self.mu * Z**2)
        elif Z < self.Z_b:
            return self.chi_s * Z * (1 + self.labda_s * Z)
        else:
            return 1.0

    def _fp_bar(self, Z, l_bar, v_bar):
        psi = self._fpsi(Z)
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

    def _ode_t(self, t_bar, Z, l_bar, v_bar):
        """time domain ode of internal ballistics"""
        p_bar = self._fp_bar(Z, l_bar, v_bar)
        if Z < self.Z_b:
            dZ = (self.theta / (2 * self.B)) ** 0.5 * p_bar**self.n  # dt_bar
        else:
            dZ = 0
        dl_bar = v_bar  # over dt_bar
        dv_bar = self.theta * 0.5 * p_bar  # over dt_bar

        return (dZ, dl_bar, dv_bar)

    def _ode_l(self, l_bar, t_bar, Z, v_bar):
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition."""
        p_bar = self._fp_bar(Z, l_bar, v_bar)
        if Z < self.Z_b:
            dZ = (self.theta / (2 * self.B)) ** 0.5 * p_bar**self.n / v_bar
        else:
            dZ = 0
        dv_bar = self.theta * 0.5 * p_bar / v_bar
        dt_bar = 1 / v_bar
        return (dt_bar, dZ, dv_bar)

    def _ode_Z(self, Z, t_bar, l_bar, v_bar):
        """burnout domain ode of internal ballistics"""
        p_bar = self._fp_bar(Z, l_bar, v_bar)
        if Z < self.Z_b:
            dt_bar = ((2 * self.B) / self.theta) ** 0.5 * p_bar**-self.n
            dl_bar = (
                v_bar * ((2 * self.B) / self.theta) ** 0.5 * p_bar**-self.n
            )
            dv_bar = (self.B * self.theta / 2) ** 0.5 * p_bar ** (1 - self.n)
        else:
            # technically speaking it is undefined in this area
            dt_bar = 0
            dl_bar = 0
            dv_bar = 0

        return (dt_bar, dl_bar, dv_bar)

    def _ode_v(self, v_bar, t_bar, Z, l_bar):
        p_bar = self._fp_bar(Z, l_bar, v_bar)
        if Z < self.Z_b:
            dZ = (2 / (self.B * self.theta)) ** 0.5 * p_bar ** (self.n - 1)
        else:
            dZ = 0
        dl_bar = 2 * v_bar / (self.theta * p_bar)
        dt_bar = 2 / (self.theta * p_bar)

        return (dt_bar, dZ, dl_bar)

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

    def integrate(self, steps=10, tol=1e-5, dom="time"):
        """
        Runs a full numerical solution for the gun in the specified domain sampled
        evenly at specified number of steps, using a scaled numerical tolerance as
        specified.

        tolerance is meant to be interpreted as in the worst case, each unitless result
        is allowed to be tol distance away from the real solution as specified.

        Through significant trials and errors, it was determined that for this particular
        system, the error due to compounding does not appear to be significant.
        """

        l_g_bar = self.l_g / self.l_0
        bar_data = []

        bar_data.append(
            (
                "SHOT START",
                0,
                0,
                self.Z_0,
                0,
                self.p_0 / (self.f * self.Delta),
            )
        )

        N = 1
        Delta_Z = self.Z_b - self.Z_0
        Z_i = self.Z_0

        Z_j = Z_i + Delta_Z / N
        t_bar_i, l_bar_i, v_bar_i = 0, 0, 0

        """
        Instead of letting the integrator handle the heavy lifting, we
        partition Z and integrate upwards until either barrel exit or
        burnout has been achieved. This seek-and-validate inspired
        from bisection puts the group of value denoted by subscript
        i within or on the muzzle, with the propellant either still
        burning or right on the burnout point..
        """
        while Z_i < self.Z_b:  # terminates if burnout is achieved
            # print(Z_i, Z_j, N)

            if Z_j == Z_i:
                raise ValueError(
                    "Numerical accuracy exhausted in search of exit/burnout point."
                )
            try:
                t_bar_j, l_bar_j, v_bar_j = RKF45OverTuple(
                    self._ode_Z,
                    (t_bar_i, l_bar_i, v_bar_i),
                    Z_i,
                    Z_j,
                    tol=tol,
                    termAbv=(None, l_g_bar, None),
                )

            except ValueError as e:
                raise ValueError(
                    "Unable to integrate system due to ill defined system, requiring"
                    + " vanishingly step size. Reducing tolerance limit is generally"
                    + " not useful. This is usually caused by excessive pressure spike."
                    + " Reduced propellant load or energy content,"
                    + " lower chamber load fraction, longer burn times via greater"
                    + " geometrical size of the propellant grain is suggested."
                )

            if any(
                (t_bar_i == t_bar_j, l_bar_i == l_bar_j, v_bar_i == v_bar_j)
            ):
                raise ValueError(
                    "Numerical integration stalled in search of exit/burnout point."
                )

            if l_bar_j > l_g_bar:
                if abs(l_bar_i - l_g_bar) > tol**0.5 or l_bar_i == 0:
                    N *= 2
                    Z_j = Z_i + Delta_Z / N
                else:
                    break  # l_bar_i is solved to within a tol of l_bar_g
            else:
                t_bar_i, l_bar_i, v_bar_i = t_bar_j, l_bar_j, v_bar_j
                Z_i = Z_j
                """
                this way the group of values denoted by _i is always updated
                as a group.
                """
                Z_j += Delta_Z / N
                if Z_j > self.Z_b:
                    Z_j = self.Z_b

        if t_bar_i == 0:
            raise ValueError("exit/burnout point found to be at the origin.")

        """
        # verify the cumulative error is small enough
        t_bar_v, l_bar_v, v_bar_v = RKF45OverTuple(
            self._ode_Z,
            (0, 0, 0),
            self.Z_0,
            Z_i,
            tol=tol,
            termAbv=(None, l_g_bar, None),
        )
        print(t_bar_v - t_bar_i, l_bar_v - l_bar_i, v_bar_v - v_bar_i)
        """
        """
        Subscript e indicate exit condition.
        At this point, since its guaranteed that point i will be further
        towards the chamber of the firearm than point e, we do not have
        to worry about the dependency of the correct behaviour of this ODE
        on the positive direction of integration.
        """
        t_bar_e, Z_e, v_bar_e = RKF45OverTuple(
            self._ode_l,
            (t_bar_i, Z_i, v_bar_i),
            l_bar_i,
            l_g_bar,
            tol=tol,
        )
        bar_data.append(
            (
                "SHOT EXIT",
                t_bar_e,
                l_g_bar,
                Z_e,
                v_bar_e,
                self._fp_bar(Z_e, l_g_bar, v_bar_e),
            )
        )

        t_bar_b = None
        if Z_e > 1:  # fracture point is contained:
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile
            movement to charge fracture
            """
            t_bar_f, l_bar_f, v_bar_f = RKF45OverTuple(
                self._ode_Z,
                (0, 0, 0),
                self.Z_0,
                1,
                tol=tol,
            )

            bar_data.append(
                (
                    "FRACTURE",
                    t_bar_f,
                    l_bar_f,
                    1,
                    v_bar_f,
                    self._fp_bar(1, l_bar_f, v_bar_f),
                )
            )

        if Z_e >= self.Z_b:
            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to Z_b, from onset of projectile
            movement to charge burnout.
            """

            t_bar_b, l_bar_b, v_bar_b = RKF45OverTuple(
                self._ode_Z,
                (0, 0, 0),
                self.Z_0,
                self.Z_b,
                tol=tol,
            )
            bar_data.append(
                (
                    "BURNOUT",
                    t_bar_b,
                    l_bar_b,
                    self.Z_b,
                    v_bar_b,
                    self._fp_bar(self.Z_b, l_bar_b, v_bar_b),
                )
            )

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.

        we hereby simply solve p golden section searching it from origin
        to point e, i.e. inside the barrel.
        """

        def f(t_bar):
            Z, l_bar, v_bar = RKF45OverTuple(
                self._ode_t,
                (self.Z_0, 0, 0),
                0,
                t_bar,
                tol=tol,
            )
            return self._fp_bar(Z, l_bar, v_bar)

        # tolerance is specified a bit differently for gold section search
        t_bar_p_1, t_bar_p_2 = gss(
            f,
            0,
            t_bar_e if t_bar_b is None else t_bar_b,
            tol=tol,
            findMin=False,
        )
        t_bar_p = (t_bar_p_1 + t_bar_p_2) / 2

        Z_p, l_bar_p, v_bar_p = RKF45OverTuple(
            self._ode_t, (self.Z_0, 0, 0), 0, t_bar_p, tol=tol
        )

        bar_data.append(
            (
                "PEAK PRESSURE",
                t_bar_p,
                l_bar_p,
                Z_p,
                v_bar_p,
                self._fp_bar(Z_p, l_bar_p, v_bar_p),
            )
        )
        """
        populate data for output purposes
        """

        bar_sample_data = []
        try:
            if dom == "time":
                for j in range(1, steps):
                    t_bar_j = t_bar_e / steps * j

                    Z_j, l_bar_j, v_bar_j = RKF45OverTuple(
                        self._ode_t,
                        (self.Z_0, 0, 0),
                        0,
                        t_bar_j,
                        tol=tol,
                    )

                    bar_sample_data.append(
                        (
                            "",
                            t_bar_j,
                            l_bar_j,
                            Z_j,
                            v_bar_j,
                            self._fp_bar(Z_j, l_bar_j, v_bar_j),
                        )
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
                t_bar_s = 0.5 * t_bar_i
                Z_s, l_bar_s, v_bar_s = RKF45OverTuple(
                    self._ode_t,
                    (self.Z_0, 0, 0),
                    0,
                    t_bar_s,
                    tol=tol,
                )

                for i in range(1, steps):
                    l_bar_j = l_g_bar / steps * i

                    t_bar_j, Z_j, v_bar_j = RKF45OverTuple(
                        self._ode_l,
                        (t_bar_s, Z_s, v_bar_s),
                        l_bar_s,
                        l_bar_j,
                        tol=tol,
                    )
                    bar_sample_data.append(
                        (
                            "",
                            t_bar_j,
                            l_bar_j,
                            Z_j,
                            v_bar_j,
                            self._fp_bar(Z_j, l_bar_j, v_bar_j),
                        )
                    )
        except ValueError:
            bar_sample_data = []
        finally:
            bar_data.extend(bar_sample_data)

        bar_data.sort(key=lambda x: x[1])  # sort by scaled time
        data = list(
            (
                tag,
                t_bar * self.l_0 / self.v_j,
                l_bar * self.l_0,
                self._fpsi(Z),
                v_bar * self.v_j,
                p_bar * self.f * self.Delta,
            )
            for (tag, t_bar, l_bar, Z, v_bar, p_bar) in bar_data
        )
        """
        t_bar_v, Z_v, l_bar_v = RKF45OverTuple(
            self._ode_v, (0, self.Z_0, 0), 0, v_bar_e, tol=tol
        )
        print(t_bar_v - t_bar_e, Z_v - self.Z_b, l_bar_v - l_g_bar)
        """
        data = list(
            (tag, t, l, psi, v, p, self._T(psi, l, p))
            for (tag, t, l, psi, v, p) in data
        )

        return data

    def analyze(self, it=250, tol=1e-5):
        """run the psi-bar analytical solution on the defined weapon

        key assumptions:
            >burn rate scales linearlly with pressure
            >propellant volume fraction (ψ, psi) is a polynomial
             function of linear burn fraction (Z),
             (this involves rewriting the psi-z polynomial to
             match the starting and end points of the cubic equation,
             which is not exactly rigorous, but is necessary
             for analytically solving the SOE.)
            >l_psi does not vary much for the entire burn duration
             (since with current propellant the increase in free
             chamber volume as a consequence of burning is out-
             weighed by the effect of covolume, this is approxima-
             tely true for moderate chamber loading fractions.)

        when u_1 refers to the linear burn rate model, it is substituted
        for u_0, deviating from the source's use.

        alternate chi, labda for quadratic form equation we use here
        to approximate the more rigorous cubic polynomial form equations.
        values are chosen such that the shot start and fracture points
        match exactly.

        Since this is an approximation, we do not specify an accuracy and
        instead supplies an iteration number. This allows it to be used
        where the calculation finishing and returning is required, such
        as say numerial optimisation routines
        """
        Z_0 = self.Z_0
        psi_s = self.chi * (1 + self.labda + self.mu)

        labda_prime = (psi_s * Z_0 - self.psi_0) / (
            self.psi_0 - Z_0**2 * psi_s
        )
        chi_prime = psi_s / (1 + labda_prime)

        # the equivalent of B used for linear burn rates.
        B_0 = (
            self.S**2
            * self.e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_0**2)
        )

        # labda_prime = self.labda
        # chi_prime = self.chi
        sigma_0 = (1 + 4 * labda_prime / chi_prime * self.psi_0) ** 0.5

        """
        # only applicable in the quadratic form case
        if labda_prime == 0:
            Z_0 = self.psi_0 / chi_prime
        else:
            Z_0 = (sigma_0 - 1) / (2 * labda_prime)
        """
        K_1 = chi_prime * sigma_0
        B_1 = B_0 * self.theta * 0.5 - chi_prime * labda_prime

        v_k = self.S * self.e_1 / (self.u_0 * self.phi * self.m)
        gamma = B_1 * self.psi_0 / K_1**2

        # x = Z - Z_0
        def _v(x):
            return v_k * x

        def _psi(x):
            """
            In effect:
            ψ(x) = ψ_0 + χ * σ_0 * x + λ * χ * x**2
            ψ(x) = ψ_0 + K_1 * x + λ * χ * x**2
            """
            return self.psi_0 + K_1 * x + labda_prime * chi_prime * x**2

        def _l_psi_avg(x_i, x_j):
            psi_i, psi_j = _psi(x_i), _psi(x_j)
            psi_avg = (psi_i + psi_j) / 2
            return self.l_0 * (
                1
                - self.Delta / self.rho_p
                - self.Delta * (self.alpha - 1 / self.rho_p) * psi_avg
            )

        def _Z_x(x):
            # distinct from Z which should be straightforwardly related to x
            beta = B_1 / K_1 * x
            b = (1 + 4 * gamma) ** 0.5
            return (1 - 2 * beta / (b + 1)) ** (0.5 * (b + 1) / b) * (
                1 + 2 * beta / (b - 1)
            ) ** (0.5 * (b - 1) / b)

        B_1_prime = -B_1
        gamma_prime = B_1_prime * self.psi_0 / K_1**2

        def _Z_prime_x(x):
            beta_prime = B_1_prime / K_1 * x
            b_prime = (1 - 4 * gamma_prime) ** 0.5
            return (1 + 2 * beta_prime / (b_prime + 1)) ** (
                0.5 * (b_prime + 1) / b_prime
            ) * (1 - 2 * beta_prime / (b_prime - 1)) ** (
                0.5 * (b_prime - 1) / b_prime
            )

        def _l(x_i, x_j, l_i):
            """
            given a bunch of parameters, calculate l_j
            """
            l_psi_avg = _l_psi_avg(x_i, x_j)
            if B_1 > 0:
                Z_x_i = _Z_x(x_i)
                Z_x_j = _Z_x(x_j)

                return (l_i + l_psi_avg) * (Z_x_j / Z_x_i) ** (
                    -B_0 / B_1
                ) - l_psi_avg

            elif B_1 == 0:
                xi_prime = K_1 / self.psi_0 * x_i
                xj_prime = K_1 / self.psi_0 * x_j

                return (
                    exp(xj_prime - xi_prime) * (1 + xi_prime) / (1 + xj_prime)
                ) ** (B * self.psi_0 / K_1**2) * (l_i + l_psi_avg) - l_psi_avg

            else:
                Z_prime_i = _Z_prime_x(x_i)
                Z_prime_j = _Z_prime_x(x_j)

                return (l_i + l_psi_avg) * (Z_prime_j / Z_prime_i) ** (
                    -B_0 / B_1
                ) - l_psi_avg

        def propagate(x, it):
            # propagate l from x= 0 to x_k (1-Z_0)
            # simulatneousely, also propagate a time.
            if x == 0:
                return 0, 0
            l_i = 0
            t = 0
            step = x / it
            x_i = 0
            for i in range(it):  # 8192
                x_j = x_i + step
                l_j = _l(x_i, x_j, l_i)
                Delta_l = l_j - l_i
                t += 2 * Delta_l / (_v(x_i) + _v(x_j))
                l_i = l_j
                x_i = x_j

            return l_i, t

        def _p(x, l):
            l_psi = _l_psi_avg(x, x)
            return (
                self.f
                * self.omega
                * (self.psi_0 + K_1 * x - B_1 * x**2)
                / (self.S * (l + l_psi))
            )

        delta_1 = 1 / (self.alpha - 1 / self.rho_p)
        x_k = 1 - Z_0  # upper limit of x_m, also known as x_s

        # values reflceting fracture point
        l_k, t_k = propagate(x_k, it=it)
        # l_k = _l(x_k, _l_psi_avg(x_k, 0))

        p_k = _p(x_k, l_k)
        v_k_1 = _v(x_k)

        """
            xi valid range  (0,1)
            psi valid range (0,1)
            Z valid range (0,Z_b)
        """

        psi_k = _psi(x_k)
        xi_k = 1 / self.Z_b

        """
            In the reference, the following terms are chi_s and labda_s.
            Issue is, this definition IS NOT the same as the definition
            self.chi_s and self.labda_s.from the previous chapter.
        """

        chi_k = (psi_k / xi_k - xi_k) / (1 - xi_k)
        labda_k = 1 - 1 / chi_k

        # print(-chi_k * labda_k, "xi**2", chi_k, "xi")

        def _psi_fracture(xi):
            """xi as in greek letter ξ
            ψ(ξ) = χ_K * ξ - λ_k * χ_K * ξ^2
            """
            return chi_k * xi * (1 - labda_k * xi)

        I_s = (self.e_1 + self.rho) / (self.u_0)
        xi_0 = Z_0 / self.Z_b

        v_kk = self.S * I_s * (xi_k - xi_0) / (self.m * self.phi)

        def _v_fracture(xi):
            return v_kk * (xi - xi_0) / (xi_k - xi_0)
            # return self.S * I_s / (self.phi * self.m) * (xi - xi_k) + v_kk

        Labda_1 = (
            l_k / self.l_0
        )  # reference is wrong, this is the implied definition
        # print(l_k)
        B_2 = self.S**2 * I_s**2 / (self.f * self.omega * self.phi * self.m)
        B_2_bar = B_2 / (chi_k * labda_k)
        xi_k_bar = labda_k * xi_k

        def _Labda(xi_i, xi_j, Labda_i):  # Labda = l / l_0
            xi_i_bar = labda_k * xi_i
            xi_j_bar = labda_k * xi_j

            Labda_psi_avg = _Labda_psi_avg(xi_i, xi_j)

            r = (
                (1 - (1 + 0.5 * B_2_bar * self.theta) * xi_j_bar)
                / (1 - (1 + 0.5 * B_2_bar * self.theta) * xi_i_bar)
            ) ** (-B_2_bar / (1 + 0.5 * B_2_bar * self.theta))

            Labda_j = r * (Labda_i + Labda_psi_avg) - Labda_psi_avg
            return Labda_j

        def _Labda_psi_avg(xi_i, xi_j):
            psi_i, psi_j = _psi_fracture(xi_i), _psi_fracture(xi_j)
            psi_avg = (psi_i + psi_j) / 2

            return (
                1
                - self.Delta / self.rho_p
                - self.Delta * (self.alpha - 1 / self.rho_p) * psi_avg
            )

        def _p_fracture(xi, Labda):
            """
            Equation is wrong for reference work
            """
            psi = _psi_fracture(xi)
            # v = _v_fracture(xi)
            Labda_psi = _Labda_psi_avg(xi, xi)

            return (
                self.f
                * self.Delta
                * (psi - 0.5 * B_2 * self.theta * (xi - xi_0) ** 2)
                / (Labda + Labda_psi)
            )
            """
            return (
                self.f * self.omega * psi
                - self.theta * self.phi * self.m * v**2 * 0.5
            ) / (self.S * (Labda_psi + Labda) * self.l_0)
            """

        def propagate_fracture(xi, it):
            # propagate Labda in the fracture regime
            Labda_i = Labda_1
            t = t_k
            step = (xi - xi_k) / it
            xi_i = xi_k
            for i in range(it):
                xi_j = xi_i + step
                Labda_j = _Labda(xi_i, xi_j, Labda_i)
                Delta_Labda = Labda_j - Labda_i
                xi_i = xi_j
                t += (
                    2
                    * Delta_Labda
                    * self.l_0
                    / (_v_fracture(xi_i) + _v_fracture(xi_j))
                )

            return Labda_j, t

        # values reflecting end of burn
        Labda_2, t_k_2 = propagate_fracture(1, it=it)
        # Labda_2 = _Labda(1, _Labda_psi_avg(1, xi_k))
        l_k_2 = Labda_2 * self.l_0
        p_k_2 = _p_fracture(1, Labda_2)
        v_k_2 = _v_fracture(1)
        psi_k_2 = _psi_fracture(1)

        """
            Defining equations for post burnout, adiabatic
            expansion phase.
        """
        # l_1 is the the limit for l_psi_avg as psi -> 0
        l_1 = self.l_0 * (1 - self.alpha * self.Delta)

        def _p_adb(l):
            return p_k_2 * ((l_1 + l_k) / (l_1 + l)) ** (self.theta + 1)

        def _v_adb(l):
            return (
                self.v_j
                * (
                    1
                    - ((l_1 + l_k) / (l_1 + l)) ** self.theta
                    * (1 - (v_k_2 / self.v_j) ** 2)
                )
                ** 0.5
            )

        def _t_adb(l, it):
            t = t_k_2
            step = (l - l_k_2) / it
            l_i = l_k_2
            t = t_k_2
            for i in range(it):  # 8192
                l_j = l_i + step
                t += 2 * step / (_v_adb(l_i) + _v_adb(l_j))
                l_i = l_j
            return t

        """
            Solving for muzzle exit condition
            """

        l_e = self.l_g
        Labda_e = l_e / self.l_0

        data = []
        if l_e >= l_k_2:  # shot exit happened after burnout
            v_e = _v_adb(l_e)
            p_e = _p_adb(l_e)
            psi_e = 1
            t_e = _t_adb(l_e, it=it)

        elif l_e >= l_k:  # shot exit happened after fracture
            xi_e = bisect(
                lambda xi: propagate_fracture(xi, it=it)[0] - Labda_e,
                # lambda xi: _Labda(xi, _Labda_psi_avg(xi, xi_k)) - Labda_e,
                xi_k,
                1,
                tol=tol,
            )[0]

            _, t_e = propagate_fracture(xi_e, it=it)
            v_e = _v_fracture(xi_e)
            p_e = _p_fracture(xi_e, Labda_e)
            psi_e = _psi_fracture(xi_e)

        else:  # shot exit happend before fracture
            x_e = bisect(
                lambda x: propagate(x, it=it)[0] - l_e,
                # lambda x: _l(x, _l_psi_avg(x, 0)) - l_e,
                0,
                x_k,
                tol=tol,
            )[0]
            _, t_e = propagate(x_e, it=it)
            v_e = _v(x_e)
            p_e = _p(x_e, l_e)
            psi_e = _psi(x_e)

        def findPeak(x):
            if x < x_k:  # pre burnout
                l, t = propagate(x, it=it)
                p = _p(x, l)
                psi = _psi(x)
                v = _v(x)
            else:
                xi = x / self.Z_b
                Labda, t = propagate_fracture(xi, it=it)
                l = Labda * self.l_0
                p = _p_fracture(xi, Labda)
                psi = _psi_fracture(xi)
                v = _v_fracture(xi)

            return (t, l, psi, v, p)

        """
        This should be the only non-deterministic part of this routine.
        """
        x_p_1, x_p_2 = gss(
            lambda x: findPeak(x)[4], 0, self.Z_b - Z_0, tol=tol, findMin=False
        )
        if x_p_2 == self.Z_b - Z_0:  # peak @ burnout
            t_m, l_m, psi_m, v_m, p_m = t_k_2, l_k_2, psi_k_2, v_k_2, p_k_2
        else:
            t_m, l_m, psi_m, v_m, p_m = findPeak((x_p_1 + x_p_2) * 0.5)

        data.append(("SHOT EXIT", t_e, l_e, psi_e, v_e, p_e))
        data.append(("SHOT START", 0.0, 0.0, self.psi_0, 0, self.p_0))
        data.append(("PEAK PRESSURE", t_m, l_m, psi_m, v_m, p_m))
        data.append(("FRACTURE", t_k, l_k, psi_k, v_k_1, p_k))
        data.append(("BURNOUT", t_k_2, l_k_2, psi_k_2, v_k_2, p_k_2))

        data = list(
            (tag, t, l, psi, v, p, self._T(psi, l, p))
            for (tag, t, l, psi, v, p) in data
        )

        data.sort(key=lambda x: x[1].real)

        return data

    def getBP(self, abortLength, abortVel, tol=1e-3):
        """routine, meant exclusively to integrate the ODE to the point
        where burnout has occured, irregardless of barrel length specified.
        Then, the peak-finding routine is ran to find the peak pressure
        point.

        """
        l_g_bar = self.l_g / self.l_0
        try:
            t_bar_b, l_bar_b, v_bar_b = RKF45OverTuple(
                self._ode_Z,
                (0, 0, 0),
                self.Z_0,
                self.Z_b,
                tol=tol,
                termAbv=(None, abortLength / self.l_0, abortVel / self.v_j),
            )

            if any(
                (
                    l_bar_b > abortLength / self.l_0,
                    v_bar_b > abortVel / self.v_j,
                )
            ):
                raise ValueError(
                    "Burnout cannot be found within the abort range specified"
                )

        except ValueError as e:
            raise e

        def f(t_bar):
            Z, l_bar, v_bar = RKF45OverTuple(
                self._ode_t,
                (self.Z_0, 0, 0),
                0,
                t_bar,
                tol=tol,
            )
            return self._fp_bar(Z, l_bar, v_bar)

        # tolerance is specified a bit differently for gold section search
        t_bar_p_1, t_bar_p_2 = gss(
            f,
            0,
            t_bar_b,
            tol=tol,
            findMin=False,
        )
        t_bar_p = (t_bar_p_1 + t_bar_p_2) / 2

        Z_p, l_bar_p, v_bar_p = RKF45OverTuple(
            self._ode_t, (self.Z_0, 0, 0), 0, t_bar_p, tol=tol
        )
        p_bar_p = self._fp_bar(Z_p, l_bar_p, v_bar_p)
        # returns the minimum barrel length required to contain the burnup.
        # the velocity at the burnup point
        # peak pressure
        return (
            l_bar_b * self.l_0,
            v_bar_b * self.v_j,
            p_bar_p * self.f * self.Delta,
        )

    @classmethod
    def constrained(
        cls,
        caliber,
        shotMass,
        propComp,
        propGeom,
        chargeMass,
        grainPR,
        grainLDR,
        startPressure,
        chamberExpansion,
        designedPress,
        designedVel,
        lengthGunMin,
        lengthGunMax,
        loadFractionMin=0.2,
        loadFractionMax=0.7,
        grainArcMin=0.1e-3,
        grainArcMax=1e-3,
        tol=1e-3,
        step=5,
    ):
        """does constrained design, finding the design space, with
        requisite grain size and
        length of gun to achieve the specified gun peak pressure and
        shot velocity, subjected to the constraint given.
        """

        if loadFractionMax >= 0.99 or loadFractionMin <= 0.01:
            raise ValueError(
                "The specified load fraction is impractical. (<=0.01/>=0.99)"
            )

        designSpace = []
        xs = []
        for i in range(step + 1):
            loadFraction = (
                i * (loadFractionMax - loadFractionMin) / step + loadFractionMin
            )
            xs.append(loadFraction)

            def f_a(a):
                print(loadFraction, a)
                prop = Propellant(propComp, propGeom, a, grainPR, grainLDR)
                cv = chargeMass / (prop.rho_p * prop.maxLF * loadFraction)
                gun = Gun(
                    caliber,
                    shotMass,
                    prop,
                    chargeMass,
                    cv,
                    startPressure,
                    1,  # this doesn't matter now.
                    chamberExpansion,
                )
                try:
                    l_b, v_b, p_p = gun.getBP(
                        abortLength=lengthGunMax,
                        abortVel=designedVel,
                        tol=tol,
                    )
                except ValueError as e:
                    """this approach is only possible since we have confidence
                    in that the valid solutions are islands surrounded by invalid
                    solutions. If instead the solution space was any other shape
                    this would have been cause for a grievious error."""

                    return 1e99

                return abs(p_p - designedPress)

            try:
                a_1, a_2 = gss(
                    f_a,
                    grainArcMin,
                    grainArcMax,
                    tol=tol * grainArcMin,
                    findMin=True,
                )
                a = 0.5 * (a_1 + a_2)
                print("solved to be", a_1, a_2)

                DeltaP = f_a(a)
                if DeltaP < designedPress * tol:
                    designSpace.append(a)
                else:
                    raise ValueError

            except ValueError as e:
                designSpace.append("x")

        from tabulate import tabulate

        print("designed vel: ", designedVel)
        print("designedPress: ", designedPress)
        print(tabulate([designSpace], headers=xs))

    # values
    def getEff(self, vg):
        """
        te: thermal efficiency
        be: ballistic efficiency
        """
        te = (vg / self.v_j) ** 2
        be = te / self.phi
        return te, be

    def getErr(self, tol):
        """
        returns the worst case maximum deviation given specified tolerance
        in domain of:
        time, length, charge burnt volume ratio, velocity, pressure.
        """
        t_err = (self.l_0 / self.v_j) * tol
        l_err = self.l_0 * tol
        v_err = self.v_j * tol

        """
        technically speaking we should use self.Z_0 here, but that must be solved
        """
        Zs = (0, -self.labda / (3 * self.mu), 1, self.Z_b)

        def dPsi(Z):
            """
            returns dPsi/dZ
            """
            if Z < 1:
                return (
                    self.chi
                    + 2 * self.labda * self.chi * Z
                    + 3 * self.mu * self.chi * Z**2
                )
            else:
                return self.chi_s + 2 * Z * self.labda_s * self.chi_s

        psi_err = max(abs(dPsi(Z) * tol) for Z in Zs)

        return (t_err, l_err, psi_err, v_err)

    def __getattr__(self, attrName):
        try:
            return getattr(self.propellant, attrName)
        except:
            raise AttributeError("object has no '%s'" % attrName)


if __name__ == "__main__":
    """standard 7 hole cylinder has d_0=e_1, hole dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """
    from tabulate import tabulate

    compositions = GrainComp.readFile("data/propellants.csv")
    M17 = compositions["M17"]

    M17SHC = Propellant(M17, Geometry.SEVEN_PERF_ROSETTE, 1.5e-3, 0, 2.5)

    # print(1 / M17SHC.rho_p / M17SHC.maxLF / 1)
    lf = 0.5
    print("DELTA:", lf * M17SHC.maxLF)
    test = Gun(
        0.035,
        1.0,
        M17SHC,
        0.8,
        0.8 / M17SHC.rho_p / M17SHC.maxLF / lf,
        30000000.0,
        3,
        1.1,
    )
    try:
        print("\nnumerical: time")
        print(
            tabulate(
                test.integrate(10, 1e-5, dom="time"),
                headers=("tag", "t", "l", "phi", "v", "p", "T"),
            )
        )
        print("\nnumerical: length")
        print(
            tabulate(
                test.integrate(10, 1e-5, dom="length"),
                headers=("tag", "t", "l", "phi", "v", "p", "T"),
            )
        )

    except ValueError as e:
        print(e)
        pass

    print("\nErrors")
    print(
        tabulate(
            [test.getErr(1e-5)],
            headers=("t", "l", "phi", "v", "p"),
        )
    )

    print("\nanalytical:")
    print(
        tabulate(
            test.analyze(it=100, tol=1e-5),
            headers=("tag", "t", "l", "phi", "v", "p"),
        )
    )
    print(test.getBP(tol=1e-3, abortVel=1200, abortLength=3.5))
    test = Gun.constrained(
        0.035,
        1.0,
        M17,
        Geometry.SEVEN_PERF_ROSETTE,
        0.8,
        0,
        2.5,
        30e6,
        1.1,
        300e6,
        2000,
        1,
        35,
    )

    # print(*test.integrate(10, 1e-5, dom="length"), sep="\n")
    # density:  lbs/in^3 -> kg/m^3, multiply by 27680
    # covolume: in^3/lbs -> m^3/kg, divide by 27680
    # force:    ft-lbs per lb ->J/kg multiply by 2.98907
    # burn rate coefficient:
    # mm/s/(MPa)**exponent -> m/s/Pa**exponent * 1e-3 /(1e6)**exponent
