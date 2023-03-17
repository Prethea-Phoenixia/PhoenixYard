import enum

import csv

from math import pi


def bisect(f, x_0, x_1, tol=1e-9, it=1000):
    """
    bisection method to numerically solve for zero for a function
    taking one variable
    two initial guesses must be of opposite sign.
    The root found is guaranteed to be within the range specified.

    tol: tolerance of f(x) if x is root
    it: maximum iteration
    """
    a = x_0
    b = x_1

    if abs(f(a)) < tol:
        return a, f(a)
    elif abs(f(b)) < tol:
        return b, f(b)
    elif f(a) * f(b) > 0:
        raise ValueError("Initial Guesses Must Be Of Opposite Sign")

    for _ in range(it):
        c = (a + b) / 2
        if abs(f(c)) < tol:
            return (c, f(c))
        if f(c) * f(a) > 0:
            a = c
        else:
            b = c

    raise ValueError("Maximum iteration exceeded at ({},{})".format(c, f(c)))


def RKF45OverTuple(dTupleFunc, iniValTuple, x_0, x_1, tol=1e-9):
    """
    Runge Kutta Fehlberg method, of the fourth and fifth order
    Even though this involves a lot more computation per cycle,
    in practice since the step size is adaptive with regard to
    the tolerance specified, significant amount of extraneous
    computation can be saved.
    """
    y_this = iniValTuple
    x = x_0
    beta = 0.9  # "safety" factor
    h = x_1 - x_0  # initial step size
    i = 0
    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        i += 1
        K1 = dTupleFunc(x, *y_this)
        K1 = tuple(k * h for k in K1)
        if any(isinstance(v, complex) for v in K1):
            h *= 0.5
            continue

        K2 = dTupleFunc(
            x + 0.2 * h, *(y + 0.2 * k1 for y, k1 in zip(y_this, K1))
        )
        K2 = tuple(k * h for k in K2)

        if any(isinstance(v, complex) for v in K2):
            h *= 0.5
            continue

        K3 = dTupleFunc(
            x + 0.375 * h,
            *(y + (3 * k1 + 9 * k2) / 32 for y, k1, k2 in zip(y_this, K1, K2))
        )
        K3 = tuple(k * h for k in K3)

        if any(isinstance(v, complex) for v in K3):
            h *= 0.5
            continue

        K4 = dTupleFunc(
            x + 12 / 13 * h,
            *(
                y + (1932 * k1 - 7200 * k2 + 7296 * k3) / 2197
                for y, k1, k2, k3 in zip(y_this, K1, K2, K3)
            )
        )
        K4 = tuple(k * h for k in K4)

        if any(isinstance(v, complex) for v in K4):
            h *= 0.5
            continue

        K5 = dTupleFunc(
            x + h,
            *(
                y + 439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4
                for y, k1, k2, k3, k4 in zip(y_this, K1, K2, K3, K4)
            )
        )
        K5 = tuple(k * h for k in K5)

        if any(isinstance(v, complex) for v in K5):
            h *= 0.5
            continue

        # 20520
        K6 = dTupleFunc(
            x + 0.5 * h,
            *(
                y
                + -8 / 27 * k1
                + 2 * k2
                - 3544 / 2565 * k3
                + 1859 / 4104 * k4
                - 11 / 40 * k5
                for y, k1, k2, k3, k4, k5 in zip(y_this, K1, K2, K3, K4, K5)
            )
        )
        K6 = tuple(k * h for k in K6)

        y_next = tuple(
            y + 25 / 216 * k1 + 1408 / 2565 * k3 + 2197 / 4104 * k4 - 1 / 5 * k5
            for y, k1, k3, k4, k5 in zip(y_this, K1, K3, K4, K5)
        )  # forth order estimation
        z_next = tuple(
            y
            + 16 / 135 * k1
            + 6656 / 12825 * k3
            + 28561 / 56430 * k4
            - 9 / 50 * k5
            + 2 / 55 * k6
            for y, k1, k3, k4, k5, k6 in zip(y_this, K1, K3, K4, K5, K6)
        )  # fifth order estimation

        epsilon = (
            sum(abs(z - y) ** 2 for z, y in zip(z_next, y_next)) ** 0.5
        )  # error estimation
        if epsilon >= tol:  # error is greater than acceptable
            h *= beta * (tol / (2 * epsilon)) ** 0.2

        else:
            y_this = y_next
            x += h
            if epsilon != 0:  # sometimes the error can be estimated to be 0
                h *= (
                    beta * (tol / (2 * epsilon)) ** 0.25
                )  # apply the new best estimate

        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x

    if abs(x - x_1) > tol:
        raise ValueError(
            "Premature Termination of Integration, after {} Cycles".format(i)
        )
    return y_this


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

        self.rhoDiv = rhoDiv  # rho divided by (e1+d0/2)
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

    r_dot = u1 * p^n

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
    ):
        self.name = name
        self.desc = desc
        self.f = propellantForce
        self.alpha = covolume
        self.rho_p = density
        self.theta = redAdbIndex
        self.u1 = detVel
        self.n = pressureExp

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
                )

                composition.append(newComp)

        return {i.name: i for i in composition}


class Propellant:
    # assumed to be multi-holed propellants
    def __init__(self, composition, propGeom, arcThick, holeDia, grainLDR):
        self.composition = composition
        self.geometry = propGeom

        self.e1 = arcThick / 2
        self.d0 = holeDia

        if propGeom == Geometry.SEVEN_PERF_CYLINDER:
            self.D0 = 8 * self.e1 + 3 * self.d0  # enforce geometry constraints
            b, a = self.D0, 0

        elif propGeom == Geometry.SEVEN_PERF_ROSETTE:
            self.D0 = 8 * self.e1 + 3 * self.d0
            b, a = self.d0 + 4 * self.e1, self.d0 + 2 * self.e1

        elif propGeom == Geometry.NINETEEN_PERF_ROSETTE:
            self.D0 = 12 * self.e1 + 5 * self.d0
            b, a = self.d0 + 4 * self.e1, self.d0 + 2 * self.e1

        elif propGeom == Geometry.NINETEEN_PERF_CYLINDER:
            self.D0 = 12 * self.e1 + 5 * self.d0
            b, a = self.D0, 0

        elif propGeom == Geometry.NINETEEN_PERF_HEXAGON:
            self.D0 = 12 * self.e1 + 5 * self.d0
            b, a = self.d0 + 2 * self.e1, self.d0 + 2 * self.e1

        elif propGeom == Geometry.NINETEEN_PERF_ROUNDED_HEXAGON:
            self.D0 = 12 * self.e1 + 5 * self.d0
            b, a = self.d0 + 2 * self.e1, self.d0 + 2 * self.e1

        else:
            raise ValueError(
                "unhandled propellant geometry {}".format(propGeom)
            )

        grainLength = self.D0 * grainLDR
        self.c = grainLength / 2

        A, B, C = propGeom.A, propGeom.B, propGeom.C
        self.rho = propGeom.rhoDiv * (self.e1 + self.d0 / 2)

        Pi1 = (A * b + B * self.d0) / (2 * self.c)
        Q1 = (C * a**2 + A * b**2 - B * self.d0**2) / (2 * self.c) ** 2

        S_T = Q1 * pi * self.c**2
        n = propGeom.nHole

        self.maxLF = 1 - n * pi * self.d0**2 / (4 * S_T)

        beta = self.e1 / self.c

        # first phase of burning rate parameters, Z from 0 to 1
        self.chi = (Q1 + 2 * Pi1) / Q1 * beta
        self.labda = (
            (n - 1 - 2 * Pi1) / (Q1 + 2 * Pi1) * beta
        )  # deliberate misspell to prevent issues with python lambda keyword
        self.mu = -(n - 1) / (Q1 + 2 * Pi1) * beta**2

        self.Z_b = (
            self.e1 + self.rho
        ) / self.e1  # second phase burning Z upper limit
        Z = 1
        psi_s = self.chi * Z * (1 + self.labda * Z + self.mu * Z**2)

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
    ):
        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.propellant = propellant
        self.omega = chargeMass
        self.V0 = chamberVolume
        self.p0 = squeezePressure
        self.l_g = lengthGun

        self.Delta = self.omega / self.V0
        self.phi = 1 + self.omega / (3 * self.m)  # extra work factor
        self.v_j = (
            2 * self.f * self.omega / (self.theta * self.phi * self.m)
        ) ** 0.5
        self.l0 = self.V0 / self.S

        self.B = (
            self.S**2
            * self.e1**2
            / (self.f * self.phi * self.omega * self.m * self.u1**2)
            * (self.f * self.Delta) ** (2 * (1 - self.n))
        )

        if self.p0 == 0:
            raise ValueError(
                "p0 = 0 will cause starting issue for ODE solver in this\n"
                + " implementation. If necessary, approximate result by\n"
                + " finding the limit of p0->0"
            )

        else:
            self.psi_0 = (1 / self.Delta - 1 / self.rho_p) / (
                self.f / self.p0 + self.alpha - 1 / self.rho_p
            )

    def fPsi(self, Z):
        if Z < 1.0:
            return self.chi * Z * (1 + self.labda * Z + self.mu * Z**2)
        elif Z < self.Z_b:
            return self.chi_s * Z / self.Z_b * (1 + self.labda_s * Z / self.Z_b)
        else:
            return 1.0

    def fp_bar(self, Z, l_bar, v_bar):
        psi = self.fPsi(Z)
        l_psi_bar = (
            1
            - self.Delta / self.rho_p
            - self.Delta * (self.alpha - 1 / self.rho_p) * psi
        )

        p_bar = (
            self.f * self.omega * psi
            - 0.5 * self.theta * self.phi * self.m * (v_bar * self.v_j) ** 2
        ) / (self.S * self.l0 * (l_bar + l_psi_bar) * self.f * self.Delta)
        return p_bar

    def ode_t(self, t_bar, Z, l_bar, v_bar):
        """return d/dt_bar for the inputs in that order, first dummy variable"""
        p_bar = self.fp_bar(Z, l_bar, v_bar)
        if Z < self.Z_b:
            dZ = (self.theta / (2 * self.B)) ** 0.5 * (
                p_bar
            ) ** self.n  # over dt
        else:
            dZ = 0

        dl_bar = v_bar  # over dt_bar
        dv_bar = self.theta * 0.5 * p_bar  # over dt_bar

        return (dZ, dl_bar, dv_bar)

    def ode_l(self, l_bar, Z, t_bar, v_bar):
        """return d/dl_bar for the inputs in that order, first dummy variable
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition."""

        p_bar = self.fp_bar(Z, l_bar, v_bar)

        if Z < self.Z_b:
            dZ = (
                (self.theta / (2 * self.B)) ** 0.5 * (p_bar) ** self.n / v_bar
            )  #
        else:
            dZ = 0

        dv_bar = self.theta * 0.5 * p_bar / v_bar
        dt_bar = 1 / v_bar
        return (dZ, dt_bar, dv_bar)

    def integrate(self, steps=100, tol=1e-9):
        """
        this step is calculated to the square of the specified accuracy since
        the result can be very close to 0 if the squeeze pressure is very low

        It is imperative to calculate this to a high degree of accuracy since
        if 0 is accepted as solution for Z0 the RKF45 integrator later will not
        be able to self-start the integration.
        """
        Z_0 = bisect(lambda z: self.fPsi(z) - self.psi_0, 0, 1, tol**2)[0]
        l_g_bar = self.l_g / self.l0

        print("finding burnout")
        """
        Subscript b indicate burnout condition
        Zi < Zb, Zj < Zb initially
        Integration interval is incremented until Zb>Zj
        at which point integration interval is cut by half

        Search is concluded when Zb-Zi is less than tolerance specified.
        Since Zi is only incremented when Zj < Zb, it is guaranteed that
        subscript i indicate a point just before burnout.

        sign of Zj-Zb is not guaranteed.

        Since the length based ODE suffers from starting issue (divide by 0)
        initially, time domain integration is used.
        """
        Z_i, l_bar_i, v_bar_i = Z_0, 0, 0
        t_bar_i = 0
        t_bar_j = 10
        delta_t_bar = t_bar_j - t_bar_i  # initialize the search interval
        burnoutFound = True
        while self.Z_b - Z_i >= tol:
            Z_j, l_bar_j, v_bar_j = RKF45OverTuple(
                self.ode_t, (Z_i, l_bar_i, v_bar_i), t_bar_i, t_bar_j, tol
            )
            if l_bar_j > l_g_bar:
                burnoutFound = False
                break

            if Z_j >= self.Z_b:  # burnout has already happened
                delta_t_bar /= 2  # shrink the interval by half
                t_bar_j = t_bar_i + delta_t_bar
            else:  # burnout has not happened, search further
                t_bar_i = t_bar_j  # increment the search interval
                t_bar_j += delta_t_bar
                Z_i, l_bar_i, v_bar_i = Z_j, l_bar_j, v_bar_j

        if burnoutFound:
            t_bar_b, l_bar_b, v_bar_b = t_bar_i, l_bar_i, v_bar_i

        print("finding exit condition")
        """
        Subscript e indicate exit condition
        integrate forth or backwards to the barrel exit condition, along
        length. Instead of Zb, Zi < Zb ~ tol is used to allow for correct
        handling of burning in the reverse direction. 
        """
        Z_e, t_bar_e, v_bar_e = RKF45OverTuple(
            self.ode_l,
            (Z_i, t_bar_i, v_bar_i),
            l_bar_i,
            l_g_bar,
            tol,
        )  # e for exit condition

        print("populating")
        t_bar_data = []

        t_bar_i, Z_i, l_bar_i, v_bar_i = 0, Z_0, 0, 0

        t_bar_data.append(
            (t_bar_i, l_bar_i, Z_i, v_bar_i, self.p0 / (self.f * self.Delta))
        )
        for i in range(steps):
            t_bar_i = (t_bar_e - 0) / steps * i
            t_bar_i2 = (t_bar_e - 0) / steps * (i + 1)
            Z_i, l_bar_i, v_bar_i = RKF45OverTuple(
                self.ode_t, (Z_i, l_bar_i, v_bar_i), t_bar_i, t_bar_i2, tol
            )

            t_bar_data.append(
                (
                    t_bar_i2,
                    l_bar_i,
                    Z_i,
                    v_bar_i,
                    self.fp_bar(Z_i, l_bar_i, v_bar_i),
                )
            )

        t_data = tuple(
            (
                t_bar * self.l0 / self.v_j,
                l_bar * self.l0,
                self.fPsi(Z),
                v_bar * self.v_j,
                p_bar * self.f * self.Delta,
            )
            for (t_bar, l_bar, Z, v_bar, p_bar) in t_bar_data
        )

        return t_data
        # ---------------

    def __getattr__(self, attrName):
        try:
            return getattr(self.propellant, attrName)
        except:
            raise AttributeError("object has no '%s'" % attrName)


if __name__ == "__main__":
    """standard 7 hole cylinder has d0=e1, hole dia = 0.5 * arc width
    D0 = 4*2e1+3*d0 = 11 * e1
    """

    compositions = GrainComp.readFile("propellants.csv")
    M17 = compositions["M17 JAN-PD-26"]

    M17SHC = Propellant(
        M17, Geometry.SEVEN_PERF_ROSETTE, 5.5e-3 * 2, 5.5e-3, 2.25
    )

    test = Gun(
        50e-3,
        1,
        M17SHC,
        1,
        1 / M17SHC.rho_p / M17SHC.maxLF / 0.999,
        10,
        3500e-3,
    )

    # print(1 / M17SHC.rho_p / M17SHC.maxLF / 1)
    # test = Gun(57e-3, 2.8, M17SHC, 1.16, 1.51e-3, 3e7, 3.624)

    print(*test.integrate(100, 1e-9), sep="\n")

    # lbs/in^3 -> kg/m^3, multiply by 27680
    # in^3/lbs -> m^3/kg, divide by 27680
    # ft-lbs per lb ->J/kg multiply by 2.98907
