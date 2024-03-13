from ballistics import Highlow
from ballistics import SimpleGeometry, MultPerfGeometry
from ballistics import POINT_PEAK_BLEED, POINT_PEAK_AVG, DOMAIN_TIME, DOMAIN_LENG
from ballistics import GrainComp, Propellant

from tabulate import tabulate
from math import pi

if __name__ == "__main__":
    """standard 7 port cylinder has d_0=e_1, port dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """

    compositions = GrainComp.readFile("ballistics/resource/propellants.csv")

    M17 = compositions["M17"]
    M1 = compositions["M1"]

    M17C = Propellant(M17, SimpleGeometry.CYLINDER, None, 25)
    M1C = Propellant(M1, SimpleGeometry.CYLINDER, None, 10)
    lf = 0.5
    print("DELTA/rho:", lf)
    test = Highlow(
        caliber=0.082,
        shotMass=5,
        propellant=M17C,
        grainSize=5e-3,
        chargeMass=0.5,
        chamberVolume=0.5 / M1C.rho_p / lf,
        expansionVolume=0.5 / M1C.rho_p / lf,
        startPressure=5e6,
        burstPressure=10e6,
        lengthGun=3.5,
        portArea=0.5 * (pi * 0.082**2 * 0.25),
        chambrage=1,
    )
    record = []

    print("\nnumerical: time")
    print(
        tabulate(
            test.integrate(10, 1e-3, dom=DOMAIN_TIME).getRawTableData(),
            headers=(
                "tag",
                "t",
                "l",
                "psi",
                "v",
                "p1",
                "pb",
                "p2",
                "ps",
                "T1",
                "T2",
                "eta",
            ),
        )
    )

    # input()

    print("\nnumerical: length")
    print(
        tabulate(
            test.integrate(9, 1e-3, dom=DOMAIN_LENG).getRawTableData(),
            headers=(
                "tag",
                "t",
                "l",
                "psi",
                "v",
                "p1",
                "pb",
                "p2",
                "ps",
                "T1",
                "T2",
                "eta",
            ),
        )
    )
    # print(test.getEff(942))
