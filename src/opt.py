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

    def f(wTom, lf):
        """
        w (charge weight) and lf (load fraction) are supplied
        lg (length of gun) and a (arc thickness) should be solved
        to yield the given boudanry condition.

        Specifically,
        """
        w = shotMass * wTom

        # first, we need to solve for a that allows the peak pressure to exist
        def fa(a):
            prop = Propellant(composition, geometry, a, a / arcToPerf, lenToDia)
            cv = w / prop.rho_p / (prop.maxLF * lf)  # chamber volume
            gun = Gun(cal, shotMass, prop, w, cv, startPress, 1, expRatio)
            return gun.propagate(tol=tolerance, maxiter=maxiter)

        a = secant(lambda a: fa(a)[0] - peakPres, 0.1e-3, 2e-3)[0]
        _, vp, lp = fa(a)

        if vp > shotVel:
            raise ValueError(
                "Velocity at peak pressure exceeds specified gun vel"
                + " , contained solutions does not exist for this set of parameters"
            )

        def flg(lg):
            print(lg / lp)
            prop = Propellant(composition, geometry, a, a / arcToPerf, lenToDia)
            cv = w / prop.rho_p / (prop.maxLF * lf)  # chamber volume
            gun = Gun(cal, shotMass, prop, w, cv, startPress, 1, expRatio)
            return gun.propagate(tol=tolerance, maxiter=maxiter, l_g=lg)

        lg = secant(lambda lg: flg(lg) - shotVel, lp * 2, lp * 3, x_min=lp)[0]

        print(a, lg)

        """
        # tag, t, l, psi, v, p
        _, _, _, _, _, pb = data[
            tuple(d[0] for d in data).index("PEAK PRESSURE")
        ]
        _, _, _, _, ve, _ = data[tuple(d[0] for d in data).index("SHOT EXIT")]

        return (pb - peakPres) ** 2 + (ve - shotVel) ** 2
        """

    print(f(1, 0.5))


if __name__ == "__main__":
    compositions = GrainComp.readFile("data/propellants.csv")
    M17 = compositions["M17"]

    print(
        guideMap(
            1000, 250e6, 50e-3, 1, 30e6, 1.1, M17, Geometry.SEVEN_PERF_ROSETTE
        )
    )
