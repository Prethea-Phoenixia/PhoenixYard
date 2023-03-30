import enum

import csv

from math import pi, log
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

    SEVEN_PERF_CYLINDER = "7 Perf Cylinder", 1, 7, 0, 0.2956, 7
    SEVEN_PERF_ROSETTE = "7 Perf Rosette", 2, 8, 12 * 3**0.5 / pi, 0.1547, 7
    NINETEEN_PERF_ROSETTE = (
        "19 Perf Rosette",
        3,
        21,
        36 * 3**0.5 / pi,
        0.1547,
        19,
    )
    NINETEEN_PERF_CYLINDER = "19 Perf Cylinder", 1, 19, 0, 0.3559, 19
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

        self.maxLF = 1 - n * pi * self.d_0**2 / (4 * S_T)

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
            if self.psi_0 < 0:
                raise ValueError(
                    "Initial burnup fraction is solved to be negative."
                    + " Suggest reducing load fraction."
                )

    def integrate(self, steps=10, tol=1e-5, dom="time"):
        """
        Runs a full numerical solution for the gun in the specified domain sampled
        evenly at specified number of steps, using a scaled numerical tolerance as
        specified.
        """
        B = (
            self.S**2
            * self.e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_1**2)
            * (self.f * self.Delta) ** (2 * (1 - self.n))
        )

        def _fPsi(Z):
            if Z < 1.0:
                return self.chi * Z * (1 + self.labda * Z + self.mu * Z**2)
            elif Z < self.Z_b:
                return self.chi_s * Z * (1 + self.labda_s * Z)
            else:
                return 1.0

        def _fp_bar(Z, l_bar, v_bar):
            psi = _fPsi(Z)
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

        Z_0 = bisect(lambda z: _fPsi(z) - self.psi_0, 0, 1, tol)[0]
        """
        It is imperative to calculate this to a high degree of accuracy since
        if 0 is accepted as solution for Z0 the RKF45 integrator later will not
        be able to self-start the system.
        """
        l_g_bar = self.l_g / self.l_0
        if Z_0 == 0:
            raise ValueError(
                "Initial burnup solved to 0, impossible to initialize ODE."
                + " Suggest reducing tolerance."
            )
        """
        Subscript f indicate fracture condition
        ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile
        movement to charge fracture
        """
        try:
            t_bar_f, l_bar_f, v_bar_f = RKF45OverTuple(
                _ode_Z, (0, 0, 0), Z_0, 1, tol
            )
        except ValueError:
            raise ValueError(
                "Unable to integrate to propellant fracture within"
                + " reasonable cycles."
                + " Propellant size might be excessive."
                + " Suggest reducing web thickness or specifying faster"
                + " burning propellant."
            )

        """
        Subscript b indicated burnout condition
        ODE w.r.t Z is integrated from 1 to Z_b, from onset of projectile
        movement to charge burnout. THe idea here is to reuse some 
        calculation from before to reduce computational burden
        """
        t_bar_b, l_bar_b, v_bar_b = RKF45OverTuple(
            _ode_Z, (t_bar_f, l_bar_f, v_bar_f), 1, self.Z_b, tol
        )
        """
        Subscript e indicate exit condition
        integrate forth or backwards to the barrel exit condition, along
        length. Instead of Zb, Zi < Zb is used to allow for correct
        handling of burning in the reverse direction. 
        """
        t_bar_e, Z_e, v_bar_e = RKF45OverTuple(
            _ode_l, (t_bar_f, 1, v_bar_f), l_bar_f, l_g_bar, tol
        )
        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.
        """

        def f(t_bar):
            Z, l_bar, v_bar = RKF45OverTuple(_ode_t, (Z_0, 0, 0), 0, t_bar, tol)
            return _fp_bar(Z, l_bar, v_bar)

        t_bar_p_1, t_bar_p_2 = gss(f, 0, t_bar_e, tol=tol, findMin=False)
        t_bar_p = (t_bar_p_1 + t_bar_p_2) / 2
        if abs(t_bar_p - t_bar_e) < 2 * tol:
            t_bar_p, Z_p, l_bar_p, v_bar_p = t_bar_e, Z_e, l_g_bar, v_bar_e
        else:
            Z_p, l_bar_p, v_bar_p = RKF45OverTuple(
                _ode_t, (Z_0, 0, 0), 0, t_bar_p, tol
            )
        """
        populate data for output purposes
        """
        bar_data = []

        bar_data.append(
            (
                "SHOT START",
                0,
                0,
                Z_0,
                0,
                self.p_0,
            )
        )

        if dom == "time":
            t_bar_i, Z_i, l_bar_i, v_bar_i = 0, Z_0, 0, 0
            for i in range(1, steps):
                t_bar_i2 = t_bar_e / steps * i

                Z_i, l_bar_i, v_bar_i = RKF45OverTuple(
                    _ode_t,
                    (Z_i, l_bar_i, v_bar_i),
                    t_bar_i,
                    t_bar_i2,
                    tol,
                )

                bar_data.append(
                    (
                        "",
                        t_bar_i2,
                        l_bar_i,
                        Z_i,
                        v_bar_i,
                        _fp_bar(Z_i, l_bar_i, v_bar_i),
                    )
                )
                t_bar_i = t_bar_i2

        else:
            t_bar_i, Z_i, l_bar_i, v_bar_i = t_bar_f, 1, l_bar_f, v_bar_f
            for i in range(1, steps):
                l_bar_i2 = l_g_bar / steps * i

                t_bar_i, Z_i, v_bar_i = RKF45OverTuple(
                    _ode_l,
                    (t_bar_i, Z_i, v_bar_i),
                    l_bar_i,
                    l_bar_i2,
                    tol,
                )
                bar_data.append(
                    (
                        "",
                        t_bar_i,
                        l_bar_i2,
                        Z_i,
                        v_bar_i,
                        _fp_bar(Z_i, l_bar_i2, v_bar_i),
                    )
                )
                l_bar_i = l_bar_i2

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
        bar_data.sort(key=lambda x: x[1])  # sort by scaled time
        data = tuple(
            (
                tag,
                t_bar * self.l_0 / self.v_j,
                l_bar * self.l_0,
                _fPsi(Z),
                v_bar * self.v_j,
                p_bar * self.f * self.Delta,
            )
            for (tag, t_bar, l_bar, Z, v_bar, p_bar) in bar_data
        )

        return data

    def getEff(self, vg):
        te = (vg / self.v_j) ** 2
        be = te / self.phi
        return te, be

    def analyze(self, tol=1e-5):
        """run the psi-bar analytical solution on the defined weapon

        key assumptions:
            >burn rate scales linearlly with pressure
            >propellant volume fraction (ψ, psi) is a polynomial
             function of linear burn fraction (Z),
             (this involves discarding the 3rd order term entire-
             ly, which is not exactly rigorous, but is necessary
             for analytically solving the SOE.)
            >l_psi does not vary much for the entire burn duration
             (since with current propellant the increase in free
             chamber volume as a consequence of burning is out-
             weighed by the effect of covolume, this is approxima-
             tely true for moderate chamber loading fractions.)

        when u_1 refers to the linear burn rate model, it is substituted
        for u_0, deviating from the source's use.


        """
        # the equivalent of B used for linear burn rates.
        B_0 = (
            self.S**2
            * self.e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_0**2)
        )
        epsilon_0 = (1 + 4 * self.labda / self.chi * self.psi_0) ** 0.5
        Z_0 = (epsilon_0 - 1) / (2 * self.labda)
        K_1 = self.chi * epsilon_0
        B_1 = B_0 * self.theta * 0.5 - self.chi * self.labda
        v_k = self.S * self.e_1 / (self.u_0 * self.phi * self.m)
        gamma = B_1 * self.psi_0 / K_1**2

        # x = Z - Z_0
        def _v(x):
            return v_k * x

        def _psi(x):
            return self.psi_0 + K_1 * x + self.labda * self.chi * x**2

        def _l_psi_avg(x_i, x_j):
            psi_i, psi_j = _psi(x_i), _psi(x_j)
            psi_avg = (psi_i + psi_j) / 2
            return self.l_0 * (
                1
                - self.Delta / self.rho_p
                - self.Delta * (self.alpha - 1 / self.rho_p) * psi_avg
            )

        def _Z_x(
            x,
        ):  # distinct from Z which should be straightforwardly related to x
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

        def _l(x, l_psi_avg):
            if self.chi * self.labda < 0.5 * B_0 * self.theta:
                return l_psi_avg * (_Z_x(x) ** (-B_0 / B_1) - 1)
            elif self.chi * self.labda == 0.5 * B_0 * self.theta:
                x_prime = K_1 / self.psi_0 * x
                return (
                    B_0 * self.psi_0 / K_1**2 * (x_prime - log(1 + x_prime))
                )
            else:
                return l_psi_avg * (_Z_prime_x(x) ** (B_0 / B_1) - 1)

        def propagate(x, tol):
            """propagate l from 0 to x using global adaptive step control"""
            l = 0
            for i in range(10):  # 1024 steps limit is quite excessive
                step = x / 2**i
                x_i = 0
                l_i = 0
                for j in range(2**i):
                    x_j = step * (j + 1)
                    l_psi_avg = _l_psi_avg(x_i, x_j)
                    l_i += _l(x_j, l_psi_avg) - _l(x_i, l_psi_avg)
                    x_i = x_j
                if abs(l_i - l) > tol:  # change in iteration
                    l = l_i
                    continue
                else:
                    return l_i
            raise ValueError(
                "Unable to propagate system to required accuracy "
                + "in reasonable cycles."
            )

        def _p(x, l):
            l_psi = _l_psi_avg(x, x)
            return (
                self.f
                * self.omega
                * (self.psi_0 + K_1 * x + B_1 * x**2)
                / (self.S * (l + l_psi))
            )

        delta_1 = 1 / (self.alpha - 1 / self.rho_p)
        x_k = 1 - Z_0  # upper limit of x_m

        def _x_m(p_m):
            """
            x_m being a function of p_m, therefore must be solved iteratively
            """
            return K_1 / (
                B_0 * (1 + self.theta) / (1 + p_m / (self.f * delta_1))
                - 2 * self.chi * self.labda
            )

        """
        iteratively solve x_m in the range of (0,x_k)
        x_k signifies end of progressive burn
        
        Since it is known that this converge extremely vigorously,
        propagation will be done using reduced accuracy

        tolerance is specified against the unitless scaled value
        for each parameter.
        """
        p_m_i = 100e6  # initial guess, 100MPa
        for i in range(10):
            x_m_i = _x_m(p_m_i)

            if x_m_i > x_k:
                x_m_i = x_k
            elif x_m_i < 0:
                x_m_i = 0

            l_m_i = propagate(x_m_i, tol=self.l_0 * tol**0.5)  # see above

            p_m_j = _p(x_m_i, l_m_i)

            if abs(p_m_i - p_m_j) > self.f * self.Delta * tol:
                p_m_i = p_m_j
            else:
                break
        if i == 10:
            raise ValueError(
                "Unable to iteratively solve for peak pressure"
                + " within reasonable cycles"
            )

        p_m = p_m_j  # peak pressure has been found!
        x_m = x_m_i

    def __getattr__(self, attrName):
        try:
            return getattr(self.propellant, attrName)
        except:
            raise AttributeError("object has no '%s'" % attrName)


if __name__ == "__main__":
    """standard 7 hole cylinder has d_0=e_1, hole dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """

    compositions = GrainComp.readFile("data/propellants.csv")
    M17 = compositions["M17 JAN-PD-26"]

    M17SHC = Propellant(M17, Geometry.SEVEN_PERF_ROSETTE, 3e-3, 1.5e-3, 3)

    # print(1 / M17SHC.rho_p / M17SHC.maxLF / 1)
    test = Gun(
        0.2, 50.0, M17SHC, 40.0, 0.05194182390081051, 30000000.0, 12.5, 1.1
    )
    print(*test.integrate(10, 1e-5, dom="time"), sep="\n")
    # print(*test.integrate(10, 1e-5, dom="length"), sep="\n")
    test.analyze(1e-9)

    # lbs/in^3 -> kg/m^3, multiply by 27680
    # in^3/lbs -> m^3/kg, divide by 27680
    # ft-lbs per lb ->J/kg multiply by 2.98907
