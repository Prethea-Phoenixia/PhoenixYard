from ballistics import ConstrainedHighlow
from ballistics import SimpleGeometry, MultPerfGeometry
from ballistics import POINT_PEAK_BLEED, POINT_PEAK_AVG
from ballistics import GrainComp, Propellant

from tabulate import tabulate
from math import pi

if __name__ == "__main__":

    compositions = GrainComp.readFile("ballistics/resource/propellants.csv")

    M17 = compositions["M17"]
    M1 = compositions["M1"]

    M1_19HEX = Propellant(M1, MultPerfGeometry.NINETEEN_PERF_HEXAGON, 1, 2.5)
    test = ConstrainedHighlow(
        caliber=0.05,
        propellant=M1_19HEX,
        shotMass=1.0,
        burstPressure=50000000.0,
        startPressure=30000000.0,
        dragCoefficient=0.03,
        chambrage=1.5,
        tol=0.001,
        designHighPressure=350000000.0,
        designLowPressure=70000000.0,
        designVelocity=1000.0,
        minWeb=1e-06,
        maxLength=100.0,
        maxEV=1,
        ambientRho=0,
        ambientP=0,
        ambientGamma=1,
        control=POINT_PEAK_AVG,
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
        loadFraction=0.49311307124633463,
        chargeMassRatio=0.5,
        portArea=0.0014726215563702157,
        lengthGun=None,
        knownBore=False,
        suppress=True,
    )

    # result = test.findMinV(
    #     chargeMassRatio=0.5 / 5,
    #     portArea=0.5 * pi * 0.082**2 * 0.25,
    # )

    print(result)
