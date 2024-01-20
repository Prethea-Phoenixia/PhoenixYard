from highlow import Highlow
from pso import pso
from math import inf

# from gun import
from gun import (
    POINT_START,
    POINT_PEAK_AVG,
    POINT_PEAK_SHOT,
    POINT_FRACTURE,
    POINT_BURNOUT,
    POINT_EXIT,
)
from highlow import POINT_PEAK_HIGH, POINT_PEAK_BLEED


class psoHighlow:
    def __init__(
        self,
        caliber,
        propellant,
        shotMass,
        chargeMass,
        # chamberVolume,
        burstPressure,  # low pressure chamber starting pressure
        startPressure,  # shot start pressure
        dragCoefficient,
        chambrage,
        tol,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        **_,
    ):
        # cache the constants:
        self.caliber = caliber
        self.propellant = propellant
        self.shotMass = shotMass
        self.chargeMass = chargeMass
        # self.chamberVolume = chamberVolume
        self.burstPressure = burstPressure
        self.startPressure = startPressure

        self.chambrage = chambrage

        self.dragCoefficient = dragCoefficient

        self.tol = tol

        self.ambientRho = ambientRho
        self.ambientP = ambientP
        self.ambientGamma = ambientGamma

    def constrained(
        self,
        minWeb,
        maxWeb,
        minCV,
        maxCV,
        minEV,
        maxEV,
        minPortA,
        maxPortA,
        minLength,
        maxLength,
        control,  # targeted pressure
        designHighPressure,
        designLowPressure,
        designVelocity,
        iw=0.8,
        cog=0.5,
        soc=0.5,
        n=50,
    ):
        self.control = control

        self.designHighPressure = designHighPressure
        self.designLowPressure = designLowPressure
        self.designVelocity = designVelocity

        self.show = False

        constraints = [
            (minWeb, maxWeb),
            (minCV, maxCV),
            (minEV, maxEV),
            (minPortA, maxPortA),
            (minLength, maxLength),
        ]
        solution, dev = pso(
            self._f,
            constraints,
            # y_rel_tol=self.tol,
            # x_rel_tol=self.tol,
            y=0,
            y_abs_tol=self.tol,
            iw=iw,
            n=n,
            cog=cog,
            soc=soc,
        )
        self.show = True
        self._f(*solution)

        for (s_min, s_max), p in zip(constraints, solution):
            print(s_min, p, s_max)

        print(dev)

    def _f(self, grainSize, chamberVolume, expansionVolume, portArea, lengthGun):
        try:
            highlow = Highlow(
                caliber=self.caliber,
                shotMass=self.shotMass,
                propellant=self.propellant,
                grainSize=grainSize,
                chargeMass=self.chargeMass,
                # chamberVolume=self.chamberVolume,
                chamberVolume=chamberVolume,
                expansionVolume=expansionVolume,
                burstPressure=self.burstPressure,
                startPressure=self.startPressure,
                portArea=portArea,
                chambrage=self.chambrage,
                lengthGun=lengthGun,
                dragCoefficient=self.dragCoefficient,
            )

            data, _, _, _ = highlow.integrate(
                step=0,
                tol=self.tol,
                ambientRho=self.ambientRho,
                ambientP=self.ambientP,
                ambientGamma=self.ambientGamma,
                peaks=[self.control, POINT_PEAK_HIGH],
            )
        except Exception:
            return inf

        tags = [line[0] for line in data]

        i_high = tags.index(POINT_PEAK_HIGH)
        i_low = tags.index(self.control)
        i_vel = tags.index(POINT_EXIT)

        p_high = data[i_high][5]

        if self.control == POINT_PEAK_BLEED:
            p_low = data[i_low][6]
        elif self.control == POINT_PEAK_AVG:
            p_low = data[i_low][7]
        elif self.control == POINT_PEAK_SHOT:
            p_low = data[i_low][8]

        v_g = data[i_vel][4]

        if self.show:
            print(p_high, p_low, v_g)

        dev = (
            sum(
                v**2
                for v in (
                    (p_high - self.designHighPressure) / self.designHighPressure,
                    (p_low - self.designLowPressure) / self.designLowPressure,
                    (v_g - self.designVelocity) / self.designVelocity,
                )
            )
            ** 0.5
        )

        return dev


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
    lf = 0.5
    print("DELTA/rho:", lf)
    test = psoHighlow(
        caliber=0.082,
        propellant=M1C,
        shotMass=5,
        chargeMass=0.3,
        # chamberVolume=0.3 / M1C.rho_p / lf,
        burstPressure=10e6,
        startPressure=5e6,
        dragCoefficient=3e-3,
        chambrage=1.5,
        tol=1e-3,
    )

    result = test.constrained(
        minWeb=1e-6,
        maxWeb=1e-2,
        minCV=1e-4,
        maxCV=1e-2,
        minEV=1e-4,
        maxEV=1e-2,
        minPortA=1e-6,
        maxPortA=1e-2,
        minLength=0.01,
        maxLength=0.5,
        control=POINT_PEAK_AVG,  # targeted pressure
        designHighPressure=150e6,
        designLowPressure=50e6,
        designVelocity=75,
    )
