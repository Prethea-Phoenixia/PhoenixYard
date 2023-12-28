from opt import Constrained
from prop import GrainComp, SimpleGeometry
from prop import Propellant
from gun import Gun
from gun import DOMAIN_TIME
from gun import SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV
from gun import (
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    POINT_BURNOUT,
)
import copy

from labellines import labelLine, labelLines

compositions = GrainComp.readFile("data/propellants.csv")
Py = copy.deepcopy(compositions["M10"])  # standin for pyroxylin
# beta does not affect the actual curve....
PyTu = Propellant(Py, SimpleGeometry.TUBE, 1, 100, 0.22)

caliber = 125e-3
tol = 1e-3
dragCoefficient = 3e-2
control = POINT_PEAK_BREECH
chamberVolume = 12.27e-3  # in liters
chargeMass = 10
loadFraction = chargeMass / chamberVolume / Py.rho_p

lengthGun = 5.2  # 52dm


def f(shotMass):
    target = Constrained(
        caliber=caliber,
        shotMass=shotMass,
        propellant=PyTu,
        startPressure=30e6,
        dragCoefficient=dragCoefficient,
        designPressure=457e6,
        designVelocity=1000,
        chambrage=1.25,  # almost exact, 4 * 12.27L/ pi (1.25dm)^2 * 8dm
    )
    chargeMassRatio = 10 / shotMass
    try:
        halfWeb, _ = target.solve(
            loadFraction,
            chargeMassRatio,
            tol,
            minWeb=1e-7,
            maxLength=1e3,
            sol=SOL_LAGRANGE,
            ambientRho=1.204,
            ambientP=101.325e3,
            ambientGamma=1.4,
            control=control,
            lengthGun=lengthGun,
            known_bore=True,
        )
        gun = Gun(
            caliber=caliber,
            shotMass=target.m,
            propellant=PyTu,
            grainSize=2 * halfWeb,
            chargeMass=chargeMass,
            chamberVolume=chamberVolume,
            startPressure=target.p_0,
            lengthGun=lengthGun,
            chambrage=target.chi_k,
            dragCoefficient=dragCoefficient,
        )
        datas = gun.integrate(
            step=0,
            tol=tol,
            dom=DOMAIN_TIME,
            sol=SOL_LAGRANGE,
            ambientRho=1.204,
            ambientP=101.325e3,
            ambientGamma=1.4,
        )[0]

        try:
            i = [line[0] for line in datas].index(POINT_BURNOUT)
            burnout = datas[i][2] / lengthGun
        except ValueError:
            burnout = 1

        velocity = datas[-1][4]
        eff = gun.getEff(velocity)

        web = halfWeb * 2e3
    except ValueError as e:
        web, burnout, velocity, eff = None, None, None, None

    return shotMass, web, burnout, velocity, eff


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import multiprocessing

    parameters = []
    for i in range(int(18 / 0.2) + 1):
        shotMass = 2 + 0.2 * i
        parameters.append(shotMass)

    with multiprocessing.Pool() as pool:
        results = pool.map(f, parameters)

    ms, es, bos, vs, es = zip(*[v for v in results if v[1] is not None])

    fig, axv = plt.subplots(layout="constrained")

    axv.plot(ms, vs)

    axe = axv.twinx()

    axe.plot(ms, es)
    axe.set_ylim(0, 1)

    """
    axp = axv.twinx()

    from math import tanh, exp

    for Y in (800e6, 1000e6, 1200e6):
        ps = []
        for m, v in zip(ms, vs):
            k = m / 5.67
            LDR = k * 14.39

            aLDRr = (LDR + 3.77 * (1 - tanh((LDR - 10) / 6.89))) / (
                14.39 + 3.77 * (1 - tanh((14.39 - 10) / 6.89))
            )

            relPen = aLDRr * exp(
                -25.9 * Y / (7850 * v**2) + 25.9 * Y / (7850 * 1800**2)
            )
            ps.append(relPen)

        axp.plot(ms, ps, color="red")
    axp.set_ylim(0, 2)
    # axe = axv.twinx()
    # axe.plot(ms, es, linestyle="dashed")
    """

    plt.show()
