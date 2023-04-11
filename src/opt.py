from gun import *


def guideMap(
    shotVel,
    peakPres,
    cal,
    shotMass,
    startPress,
    expRatio,
    composition,
    geometry,
    arcToPerf=2,
    lenToDia=2.5,
    tolerance=1e-5,
    maxiter=100,
):
    """
    given velocity and peak pressure as constraints
    """

    def f(wTom, lf, lg, a):
        """
        w (charge weight) and lf (load fraction) are supplied
        lg (length of gun) and a (arc thickness) should be solved
        to yield the given boudanry condition.

        Specifically,
        """
        w = shotMass * wTom
        prop = Propellant(composition, geometry, a, a / arcToPerf, lenToDia)
        cv = w / prop.rho_p / (prop.maxLF * lf)  # chamber volume
        gun = Gun(cal, shotMass, prop, w, cv, startPress, lg, expRatio)
        data = gun.integrate(steps=0, tol=tolerance, maxiter=maxiter)

        # tag, t, l, psi, v, p
        _, _, _, _, _, pb = data[
            tuple(d[0] for d in data).index("PEAK PRESSURE")
        ]

        _, _, _, _, ve, _ = data[tuple(d[0] for d in data).index("SHOT EXIT")]

        return (pb - peakPres) ** 2 + (ve - shotVel) ** 2

    print(f(1, 0.1, 1, 1e-3))


if __name__ == "__main__":
    compositions = GrainComp.readFile("data/propellants.csv")
    M17 = compositions["M17"]

    print(
        guideMap(
            1000, 300, 50e-3, 1, 30e6, 1.1, M17, Geometry.SEVEN_PERF_ROSETTE
        )
    )
