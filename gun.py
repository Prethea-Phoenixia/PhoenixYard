import enum

import csv

from math import pi, log
from numerical import *


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
        self.V0 = chamberVolume
        self.p0 = squeezePressure
        self.l_g = lengthGun
        self.chi_k = chamberExpansion  # ration of l0 / l_chamber
        self.Delta = self.omega / self.V0
        self.l0 = self.V0 / self.S
        Labda = self.l_g / self.l0
        self.phi = 1 + self.omega / (3 * self.m) * (
            1 - (1 - 1 / self.chi_k) * 2.303 * log(Labda + 1) / Labda
        )  # extra work factor, chamberage effect averaged over entire length
        self.v_j = (
            2 * self.f * self.omega / (self.theta * self.phi * self.m)
        ) ** 0.5
        self.B = (
            self.S**2
            * self.e1**2
            / (self.f * self.phi * self.omega * self.m * self.u1**2)
            * (self.f * self.Delta) ** (2 * (1 - self.n))
        )

        if self.p0 == 0:
            raise ValueError(
                "p0 = 0 will cause starting issue for ODE solver in this"
                + " implementation. If necessary, approximate result by"
                + " finding the limit of p0->0"
            )

        else:
            self.psi_0 = (1 / self.Delta - 1 / self.rho_p) / (
                self.f / self.p0 + self.alpha - 1 / self.rho_p
            )

    def _fPsi(self, Z):
        if Z < 1.0:
            return self.chi * Z * (1 + self.labda * Z + self.mu * Z**2)
        elif Z < self.Z_b:
            return self.chi_s * Z * (1 + self.labda_s * Z)
        else:
            return 1.0

    def _fp_bar(self, Z, l_bar, v_bar):
        psi = self._fPsi(Z)
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

    def _ode_t(self, t_bar, Z, l_bar, v_bar):
        """return d/dt_bar for the inputs in that order, first dummy variable"""
        p_bar = self._fp_bar(Z, l_bar, v_bar)
        if Z < self.Z_b:
            dZ = (self.theta / (2 * self.B)) ** 0.5 * (
                p_bar
            ) ** self.n  # over dt
        else:
            dZ = 0

        dl_bar = v_bar  # over dt_bar
        dv_bar = self.theta * 0.5 * p_bar  # over dt_bar

        return (dZ, dl_bar, dv_bar)

    def _ode_l(self, l_bar, Z, t_bar, v_bar):
        """return d/dl_bar for the inputs in that order, first dummy variable
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition."""

        p_bar = self._fp_bar(Z, l_bar, v_bar)

        if Z < self.Z_b:
            dZ = (
                (self.theta / (2 * self.B)) ** 0.5 * (p_bar) ** self.n / v_bar
            )  #
        else:
            dZ = 0

        dv_bar = self.theta * 0.5 * p_bar / v_bar
        dt_bar = 1 / v_bar
        return (dZ, dt_bar, dv_bar)

    def integrate(self, steps=10, tol=1e-5, dom="time"):
        """
        this step is calculated to the square of the specified accuracy since
        the result can be very close to 0 if the squeeze pressure is very low

        It is imperative to calculate this to a high degree of accuracy since
        if 0 is accepted as solution for Z0 the RKF45 integrator later will not
        be able to self-start the integration.
        """
        Z_0 = bisect(lambda z: self._fPsi(z) - self.psi_0, 0, 1, tol)[0]
        l_g_bar = self.l_g / self.l0
        if Z_0 == 0:
            raise ValueError(
                "Initial burnup solved to 0, impossible to initialize ODE."
                + " Suggest reducing tolerance."
            )
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
                self._ode_t, (Z_i, l_bar_i, v_bar_i), t_bar_i, t_bar_j, tol
            )

            if Z_j >= self.Z_b:  # burnout has already happened
                delta_t_bar /= 2  # shrink the interval by half
                t_bar_j = t_bar_i + delta_t_bar
            else:
                Z_i, l_bar_i, t_bar_i, v_bar_i = Z_j, l_bar_j, t_bar_j, v_bar_j
                if l_bar_j > l_g_bar:  # barrel length has been exhausted
                    burnoutFound = False
                    break
                else:  # burnout has not happened, search further
                    t_bar_j += delta_t_bar  # increment search interval

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
            self._ode_l,
            (Z_i, t_bar_i, v_bar_i),
            l_bar_i,
            l_g_bar,
            tol,
        )  # e for exit condition

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.
        """

        print("finding peak pressure condition")

        def f(t_bar):
            Z, l_bar, v_bar = RKF45OverTuple(
                self._ode_t, (Z_0, 0, 0), 0, t_bar, tol
            )
            return self._fp_bar(Z, l_bar, v_bar)

        t_bar_p_1, t_bar_p_2 = gss(f, 0, t_bar_e, tol=tol, findMin=False)

        t_bar_p = (t_bar_p_1 + t_bar_p_2) / 2

        if abs(t_bar_p - t_bar_e) < 2 * tol:
            t_bar_p, Z_p, l_bar_p, v_bar_p = t_bar_e, Z_e, l_g_bar, v_bar_e
        else:
            Z_p, l_bar_p, v_bar_p = RKF45OverTuple(
                self._ode_t, (Z_0, 0, 0), 0, t_bar_p, tol
            )

        """
        populate data for output purposes
        """

        print("populating")
        bar_data = []

        bar_data.append(
            (
                "SHOT START",
                0,
                0,
                Z_0,
                0,
                self.p0 / (self.f * self.Delta),
            )
        )

        if dom == "time":
            t_bar_i, Z_i, l_bar_i, v_bar_i = 0, Z_0, 0, 0
            for i in range(steps):
                t_bar_i2 = t_bar_e / steps * (i + 1)

                if t_bar_i <= t_bar_p <= t_bar_i2:
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

                if burnoutFound:
                    if t_bar_i <= t_bar_b <= t_bar_i2:
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

                Z_i, l_bar_i, v_bar_i = RKF45OverTuple(
                    self._ode_t,
                    (Z_i, l_bar_i, v_bar_i),
                    t_bar_i,
                    t_bar_i2,
                    tol,
                )

                bar_data.append(
                    (
                        "SHOT EXIT" if i == steps - 1 else "",
                        t_bar_i2,
                        l_bar_i,
                        Z_i,
                        v_bar_i,
                        self._fp_bar(Z_i, l_bar_i, v_bar_i),
                    )
                )
                t_bar_i = t_bar_i2

        else:
            t_bar_i, Z_i, l_bar_i, v_bar_i = t_bar_p, Z_p, l_bar_p, v_bar_p

            for i in range(steps):
                l_bar_i2 = l_g_bar / steps * (i + 1)

                if l_bar_i <= l_bar_p <= l_bar_i2:
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

                if burnoutFound:
                    if l_bar_i <= l_bar_b <= l_bar_i2:
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

                Z_i, t_bar_i, v_bar_i = RKF45OverTuple(
                    self._ode_l,
                    (Z_i, t_bar_i, v_bar_i),
                    l_bar_i,
                    l_bar_i2,
                    tol,
                )
                bar_data.append(
                    (
                        "SHOT EXIT" if i == steps - 1 else "",
                        t_bar_i,
                        l_bar_i2,
                        Z_i,
                        v_bar_i,
                        self._fp_bar(Z_i, l_bar_i2, v_bar_i),
                    )
                )
                l_bar_i = l_bar_i2

        data = tuple(
            (
                tag,
                t_bar * self.l0 / self.v_j,
                l_bar * self.l0,
                self._fPsi(Z),
                v_bar * self.v_j,
                p_bar * self.f * self.Delta,
            )
            for (tag, t_bar, l_bar, Z, v_bar, p_bar) in bar_data
        )

        return data

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
        M17, Geometry.SEVEN_PERF_ROSETTE, 0.55e-3 * 2, 0.55e-3, 2.25
    )

    # print(1 / M17SHC.rho_p / M17SHC.maxLF / 1)
    test = Gun(57e-3, 2.8, M17SHC, 1.16, 1.51e-3, 3e7, 3.624, 1.2)
    print(*test.integrate(100, 1e-1, dom="time"), sep="\n")
    print(*test.integrate(100, 1e-9, dom="length"), sep="\n")

    # lbs/in^3 -> kg/m^3, multiply by 27680
    # in^3/lbs -> m^3/kg, divide by 27680
    # ft-lbs per lb ->J/kg multiply by 2.98907
