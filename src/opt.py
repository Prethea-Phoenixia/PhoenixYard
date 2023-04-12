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
        print(wTom, lf)
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

        a = bisect(lambda a: fa(a)[0] - peakPres, 1e-4, 1e-2, tol=tolerance)[
            0
        ]  # 0.1mm to 10mm

        pp, vp, lp = fa(a)
        if (pp - peakPres) / (peakPres) > tolerance:
            raise ValueError(
                "Impossible to achieve the specified design pressure"
            )

        if vp > shotVel:
            raise ValueError(
                "Velocity at peak pressure exceeds specified gun vel"
                + " , contained solutions does not exist for this set of parameters"
            )

        def flg(lg):
            prop = Propellant(composition, geometry, a, a / arcToPerf, lenToDia)
            cv = w / prop.rho_p / (prop.maxLF * lf)  # chamber volume
            gun = Gun(cal, shotMass, prop, w, cv, startPress, 1, expRatio)
            return gun.propagate(tol=tolerance, maxiter=maxiter, l_g=lg)

        lg = secant(lambda lg: flg(lg) - shotVel, lp * 2, lp * 3, x_min=lp)[0]
        print(lg)

        return (a, lg)

    aLst = []
    lgLst = []

    wToms = [round((i + 1) / 3, 2) for i in range(6)]
    lfs = [round((i + 1) / 3, 2) for i in range(3)]
    for wTom in wToms:
        aLine, lgLine = [], []
        for lf in lfs:
            try:
                a, lg = f(wTom, lf)
                aLine.append(a)
                lgLine.append(lg)
            except ValueError as e:
                aLine.append(0)
                lgLine.append(0)
        aLst.append(aLine)
        lgLst.append(lgLine)

    from tabulate import tabulate

    print(tabulate(aLst, headers=lfs))
    print(tabulate(lgLst, headers=lfs))


if __name__ == "__main__":
    compositions = GrainComp.readFile("data/propellants.csv")
    M17 = compositions["M17"]

    print(
        guideMap(
            1000, 250e6, 50e-3, 1, 30e6, 1.1, M17, Geometry.SEVEN_PERF_ROSETTE
        )
    )
