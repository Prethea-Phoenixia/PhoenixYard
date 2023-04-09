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
    ):
        self.name = name
        self.desc = desc
        self.f = propellantForce
        self.alpha = covolume
        self.rho_p = density
        self.theta = redAdbIndex
        self.u_1 = detVel
        self.n = pressureExp
        self.u_0 = linBurnCoe

    def readFile(fileName):
        composition = []
        with open(fileName, newline="") as csvfile:
            propReader = csv.reader(
                csvfile,
                delimiter=",",
                quotechar="|",
                quoting=csv.QUOTE_NONNUMERIC,
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
                ) = prop

                redAdbIndex = adb - 1

                newComp = GrainComp(
                    name,
                    desc,
                    propellantForce,
                    covolume,
                    density,
                    redAdbIndex,
                    burnRateCoe,
                    pressureExp,
                    linBurnCoe,
                )

                composition.append(newComp)

        return {i.name: i for i in composition}


class Propellant:
    # assumed to be multi-holed propellants
    def __init__(self, composition, propGeom, arcThick, holeDia, grainLDR):
        self.composition = composition
        self.geometry = propGeom

        self.e_1 = arcThick / 2
        self.d_0 = holeDia

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

        grainLength = self.D_0 * grainLDR
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
        squeezePressure,
        lengthGun,
        chamberExpansion,
    ):
        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.propellant = propellant
        self.omega = chargeMass
        self.V_0 = chamberVolume
        self.p_0 = squeezePressure
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

    def integrate(self, steps=10, tol=1e-5, dom="time"):
        """
        Runs a full numerical solution for the gun in the specified domain sampled
        evenly at specified number of steps, using a scaled numerical tolerance as
        specified.

        tolerance is interpreted as, the sum of all error in the unitless result
        is allowed to lie specified tol distance away from the real solution.

        tolerance specification is halved due to the following simplistic
        error analysis:
            |--integrating in Z to fracture: 0.5
            |--integrating in Z to burnout: 1
            |--integrating in l to exit: 0.5+0.5
            |--numerically solving peak pressure point in time: (0.5)+0.5
            |--populating in tiem: 1
            |--populating in distance: 0.5+0.5
        worst case: 1
        """

        B = (
            self.S**2
            * self.e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_1**2)
            * (self.f * self.Delta) ** (2 * (1 - self.n))
        )

        def _fpsi(Z):
            if Z < 0:
                return 0
            elif Z < 1.0:
                return self.chi * Z * (1 + self.labda * Z + self.mu * Z**2)
            elif Z < self.Z_b:
                return self.chi_s * Z * (1 + self.labda_s * Z)
            else:
                return 1.0

        def _fp_bar(Z, l_bar, v_bar):
            psi = _fpsi(Z)
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

        def _ode_t(t_bar, Z, l_bar, v_bar):
            """time domain ode of internal ballistics"""
            p_bar = _fp_bar(Z, l_bar, v_bar)
            if Z < self.Z_b:
                dZ = (self.theta / (2 * B)) ** 0.5 * p_bar**self.n  # dt_bar
            else:
                dZ = 0
            dl_bar = v_bar  # over dt_bar
            dv_bar = self.theta * 0.5 * p_bar  # over dt_bar

            return (dZ, dl_bar, dv_bar)

        def _ode_l(l_bar, t_bar, Z, v_bar):
            """length domain ode of internal ballistics
            the 1/v_bar pose a starting problem that prevent us from using it from
            initial condition."""
            p_bar = _fp_bar(Z, l_bar, v_bar)
            if Z < self.Z_b:
                dZ = (self.theta / (2 * B)) ** 0.5 * p_bar**self.n / v_bar
            else:
                dZ = 0
            dv_bar = self.theta * 0.5 * p_bar / v_bar
            dt_bar = 1 / v_bar
            return (dt_bar, dZ, dv_bar)

        def _ode_Z(Z, t_bar, l_bar, v_bar):
            """burnout domain ode of internal ballistics"""
            p_bar = _fp_bar(Z, l_bar, v_bar)
            if Z < self.Z_b:
                dt_bar = ((2 * B) / self.theta) ** 0.5 * p_bar**-self.n
                dl_bar = (
                    v_bar * ((2 * B) / self.theta) ** 0.5 * p_bar**-self.n
                )
                dv_bar = (B * self.theta / 2) ** 0.5 * p_bar ** (1 - self.n)
            else:
                # technically speaking it is undefined in this area
                dt_bar = 0
                dl_bar = 0
                dv_bar = 0

            return (dt_bar, dl_bar, dv_bar)

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

        Z_0 = Zs[0]

        l_g_bar = self.l_g / self.l_0

        bar_data = []

        bar_data.append(
            (
                "SHOT START",
                0,
                0,
                Z_0,
                0,
                self.p_0 / (self.f * self.Delta),
            )
        )

        Z_i = 1
        t_bar_i = None
        try:
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile
            movement to charge fracture
            """
            t_bar_f, l_bar_f, v_bar_f = RKF45OverTuple(
                _ode_Z, (0, 0, 0), Z_0, 1, 0.5 * tol
            )

            bar_data.append(
                (
                    "FRACTURE",
                    t_bar_f,
                    l_bar_f,
                    1,
                    v_bar_f,
                    _fp_bar(1, l_bar_f, v_bar_f),
                )
            )

            Z_i, t_bar_i, l_bar_i, v_bar_i = 1, t_bar_f, l_bar_f, v_bar_f

            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to Z_b, from onset of projectile
            movement to charge burnout.

            This is only done when fracture has happened inside the barrel.
            """
            if l_bar_f < l_g_bar:
                t_bar_b, l_bar_b, v_bar_b = RKF45OverTuple(
                    _ode_Z, (0, 0, 0), Z_0, self.Z_b, tol
                )
                bar_data.append(
                    (
                        "BURNOUT",
                        t_bar_b,
                        l_bar_b,
                        self.Z_b,
                        v_bar_b,
                        _fp_bar(self.Z_b, l_bar_b, v_bar_b),
                    )
                )

        except ValueError:
            """
            If program execution has branched here, this indicate that
            the fracture happens so far out that the RKF integrator has
            reached its cycle count limit.

            In this case, we try to salvage the calculation by repeatedly
            golden section searching the interval, until a valid position
            has been found.
            """
            while t_bar_i is None:
                Z_i = (Z_i - Z_0) * 0.618 + Z_0
                try:
                    t_bar_i, l_bar_i, v_bar_i = RKF45OverTuple(
                        _ode_Z, (0, 0, 0), Z_0, Z_i, 0.5 * tol
                    )
                except ValueError:
                    pass

        """
        Subscript e indicate exit condition
        integrate forth or backwards to the barrel exit condition, along
        length. Instead of Zb, Zf < Zb is used to allow for correct
        handling of burning in the reverse direction. 
        """
        t_bar_e, Z_e, v_bar_e = RKF45OverTuple(
            _ode_l, (t_bar_i, Z_i, v_bar_i), l_bar_i, l_g_bar, 0.5 * tol
        )

        bar_data.append(
            (
                "SHOT EXIT",
                t_bar_e,
                l_g_bar,
                Z_e,
                v_bar_e,
                _fp_bar(Z_e, l_g_bar, v_bar_e),
            )
        )

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.
        """

        def f(t_bar):
            Z, l_bar, v_bar = RKF45OverTuple(
                _ode_t, (Z_0, 0, 0), 0, t_bar, 0.5 * tol
            )
            return _fp_bar(Z, l_bar, v_bar)

        # tolerance is specified a bit differently for gold section search
        t_bar_p_1, t_bar_p_2 = gss(f, 0, t_bar_e, tol=tol, findMin=False)
        t_bar_p = (t_bar_p_1 + t_bar_p_2) / 2
        """
        so the worst case scenario for the error is basically:
                t_bar_p_1               t_bar_p_1
                   RKF                     RKF
        +--(tol/2)--|--(tol/2)--+--(tol/2)--|--(tol/2)--+
                    +-------gss|(tol)-------+
                            t_bar_p
        """
        Z_p, l_bar_p, v_bar_p = RKF45OverTuple(
            _ode_t, (Z_0, 0, 0), 0, t_bar_p, tol
        )

        bar_data.append(
            (
                "PEAK PRESSURE",
                t_bar_p,
                l_bar_p,
                Z_p,
                v_bar_p,
                _fp_bar(Z_p, l_bar_p, v_bar_p),
            )
        )
        """
        populate data for output purposes
        """

        if dom == "time":
            for i in range(1, steps):
                t_bar_i = t_bar_e / steps * i

                Z_i, l_bar_i, v_bar_i = RKF45OverTuple(
                    _ode_t,
                    (Z_0, 0, 0),
                    0,
                    t_bar_i,
                    tol,
                )

                bar_data.append(
                    (
                        "",
                        t_bar_i,
                        l_bar_i,
                        Z_i,
                        v_bar_i,
                        _fp_bar(Z_i, l_bar_i, v_bar_i),
                    )
                )
        else:
            """
            Due to the compounding error for integrating from the fracture
            point, it is necessary to use reduced tolerance to ensure the
            final result lies within acceptable error bar
            """
            for i in range(1, steps):
                l_bar_i = l_g_bar / steps * i

                t_bar_i, Z_i, v_bar_i = RKF45OverTuple(
                    _ode_l,
                    (t_bar_f, 1, v_bar_f),
                    l_bar_f,
                    l_bar_i,
                    0.5 * tol,
                )
                bar_data.append(
                    (
                        "",
                        t_bar_i,
                        l_bar_i,
                        Z_i,
                        v_bar_i,
                        _fp_bar(Z_i, l_bar_i, v_bar_i),
                    )
                )

        bar_data.sort(key=lambda x: x[1])  # sort by scaled time
        data = list(
            (
                tag,
                t_bar * self.l_0 / self.v_j,
                l_bar * self.l_0,
                _fpsi(Z),
                v_bar * self.v_j,
                p_bar * self.f * self.Delta,
            )
            for (tag, t_bar, l_bar, Z, v_bar, p_bar) in bar_data
        )

        return data

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
        p_err = self.f * self.Delta * tol

        """
        technically speaking we should use Z_0 here, but that must be solved

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

        return (t_err, l_err, psi_err, v_err, p_err)

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
    M17 = compositions["M17 JAN-PD-26"]

    M17SHC = Propellant(M17, Geometry.SEVEN_PERF_ROSETTE, 0.1e-3, 0.5e-3, 3)

    # print(1 / M17SHC.rho_p / M17SHC.maxLF / 1)
    test = Gun(
        0.035, 1.0, M17SHC, 0.8, 0.001086470935115527, 30000000.0, 1.0, 1.1
    )
    print("\nnumerical: time")
    print(
        tabulate(
            test.integrate(5, 1e-5, dom="time"),
            headers=("tag", "t", "l", "phi", "v", "p"),
        )
    )
    print("\nnumerical: length")
    print(
        tabulate(
            test.integrate(5, 1e-5, dom="length"),
            headers=("tag", "t", "l", "phi", "v", "p"),
        )
    )

    print("\nErrors")
    print(
        tabulate(
            [test.getErr(1e-5)],
            headers=("t", "l", "phi", "v", "p"),
        )
    )

    # print(*test.integrate(10, 1e-5, dom="length"), sep="\n")
    # lbs/in^3 -> kg/m^3, multiply by 27680
    # in^3/lbs -> m^3/kg, divide by 27680
    # ft-lbs per lb ->J/kg multiply by 2.98907
