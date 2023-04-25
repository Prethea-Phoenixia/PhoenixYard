from gun import *
import multiprocessing

import time


def flatten(l):
    """
    https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
    """
    return [item for sublist in l for item in sublist]


def searchChargeMass(context, chargeMass):
    # for j in range(ystep):
    # uncompress the context object
    (
        caliber,
        shotMass,
        propComp,
        propGeom,
        chargeMassMin,
        chargeMassMax,
        grainPR,
        grainLDR,
        startPressure,
        chamberExpansion,
        designedPress,
        designedVel,
        lengthGunMinCal,
        lengthGunMaxCal,
        grainArcMin,
        grainArcMax,
        tol,  # numerical tolerance, unitless
        arcTol,  # tolerance for the arc, representing manufacturing limit
        xs,
    ) = context

    aSpaceCM = []
    lSpaceCM = []

    lengthGunMin = lengthGunMinCal * caliber
    lengthGunMax = lengthGunMaxCal * caliber

    for loadFraction in xs:

        def f_a(a):
            prop = Propellant(propComp, propGeom, a, grainPR, grainLDR)
            cv = chargeMass / (prop.rho_p * prop.maxLF * loadFraction)
            gun = Gun(
                caliber,
                shotMass,
                prop,
                chargeMass,
                cv,
                startPressure,
                1,  # this doesn't matter now.
                chamberExpansion,
            )
            try:
                l_b, v_b, p_p, p_max, p_min = gun.getBP(
                    abortLength=lengthGunMax,
                    abortVel=designedVel,
                    tol=tol,
                )

            except ValueError as e:
                """this approach is only possible since we have confidence
                in that the valid solutions are islands surrounded by invalid
                solutions. If instead the solution space was any other shape
                this would have been cause for a grievious error."""

                return inf, None, None, None, e

            return abs(p_p - designedPress), gun, p_max, p_min, None

        m = 1

        try:
            DeltaP_a, gun_a, pmax_a, pmin_a, err = f_a(grainArcMax)
            if gun_a is not None:  # and pmin_a < designedPress:
                maxValida = grainArcMax
            else:
                maxValida = grainArcMin

            DeltaP_a, gun_a, pmax_a, pmin_a, err = f_a(grainArcMin)
            if gun_a is not None:  # and pmax_a > designedPress:
                minValida = grainArcMin
            else:
                minValida = grainArcMax

            t0 = time.time()

            # FIXME: This segment of code is extremely long running

            while 2**-m > arcTol / (grainArcMax - grainArcMin):
                for n in range(0, 2**m):
                    if n % 2 != 0:
                        v = (
                            n * (grainArcMax - grainArcMin) / 2**m
                            + grainArcMin
                        )

                        if maxValida > v > minValida:
                            continue

                        (_, det, pmax_v, pmin_v, _) = f_a(v)

                        if det is not None:
                            if v < minValida:
                                minValida = v

                            if v > maxValida:
                                maxValida = v

                m += 1

            t1 = time.time()

            print(loadFraction, chargeMass)
            print("validation of arc range took:", t1 - t0)

            # print(minValida, maxValida)

            if minValida >= maxValida:
                """
                DNE a valid range of arc such that:
                    burnout is contained within max length
                    designed velocity is achieved after burnout

                """
                raise ValueError("!ARC")

            # print(minValida, maxValida)

            a_1, a_2 = gss(
                lambda x: f_a(x)[0],
                minValida,
                maxValida,
                tol=tol * grainArcMin,
                findMin=True,
            )  # minimize the deviation of actual developed pressure with desired

            a = 0.5 * (a_1 + a_2)

            t2 = time.time()

            print("minimizing pressure deviation took:", t2 - t1)

            DeltaP, gun, pmax, pmin, err = f_a(a)

            if gun is None:  # this **really** should not happen
                if isinstance(err, AbortedDueToLength):
                    raise ValueError("L<B.O.")
                elif isinstance(err, AbortedDueToVelocity):
                    raise ValueError("V<B.O.")
                else:
                    raise err

            if pmax < designedPress:
                """
                the required arc:
                    to develop sufficient pressure, and:
                    to contain the burnout within max length, and:
                    to contain the designed velocity within the same.
                is smaller than required
                """
                raise ValueError("A<MIN")

            elif pmin > designedPress:
                raise ValueError("A>MAX")

            lg = gun.findBL(designedVel, lengthGunMax, tol=tol)

            if lengthGunMin > lg:
                raise ValueError("L<MIN")
            elif lg > lengthGunMax:
                raise ValueError("L>MAX")
            else:
                aSpaceCM.append(a)
                lSpaceCM.append(lg / caliber)

        except ValueError as e:
            lSpaceCM.append(str(e))
            aSpaceCM.append(str(e))

    # aSpace.append(aSpaceCM)
    # lSpace.append(lSpaceCM)
    return (
        lSpaceCM,
        aSpaceCM,
    )


def constrained(
    caliber,
    shotMass,
    propComp,
    propGeom,
    chargeMassMin,
    chargeMassMax,
    grainPR,
    grainLDR,
    startPressure,
    chamberExpansion,
    designedPress,
    designedVel,
    lengthGunMinCal,
    lengthGunMaxCal,
    loadFractionMin,
    loadFractionMax,
    grainArcMin,
    grainArcMax,
    tol,
    arcTol,
    xstep,
    ystep,
):
    """does constrained design, finding the design space, with
    requisite grain size and
    length of gun to achieve the specified gun peak pressure and
    shot velocity, subjected to the constraint given.
    """

    if loadFractionMax >= 1 or loadFractionMin <= 0:
        raise ValueError(
            "The specified load fraction is impractical. (<=0.01/>=0.99)"
        )

    aSpace = []
    lSpace = []
    xs = [
        i * (loadFractionMax - loadFractionMin) / (xstep - 1) + loadFractionMin
        for i in range(xstep)
    ]
    ys = [
        j * (chargeMassMax - chargeMassMin) / (ystep - 1) + chargeMassMin
        for j in range(ystep)
    ]

    context = (
        caliber,
        shotMass,
        propComp,
        propGeom,
        chargeMassMin,
        chargeMassMax,
        grainPR,
        grainLDR,
        startPressure,
        chamberExpansion,
        designedPress,
        designedVel,
        lengthGunMinCal,
        lengthGunMaxCal,
        grainArcMin,
        grainArcMax,
        tol,
        arcTol,
        xs,
    )
    """
    # this is the serial version
    tabuleau = flatten(searchChargeMass(context, lf) for lf in ys)
    """

    with multiprocessing.Pool() as pool:
        tabuleau = flatten(
            pool.starmap(searchChargeMass, ((context, y) for y in ys))
        )

    from tabulate import tabulate

    print("designed vel: ", designedVel)
    print("designed Pmax: ", designedPress)

    print(
        tabulate(
            tabuleau,
            headers=(*(round(x, 2) for x in xs),),
        )
    )


if __name__ == "__main__":
    """standard 7 hole cylinder has d_0=e_1, hole dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """
    from tabulate import tabulate

    compositions = GrainComp.readFile("data/propellants.csv")
    M17 = compositions["M17"]

    M17SHC = Propellant(M17, Geometry.SEVEN_PERF_ROSETTE, 1.5e-3, 0, 2.5)

    test = constrained(
        0.05,
        1.0,
        M17,
        Geometry.SEVEN_PERF_ROSETTE,
        chargeMassMin=0.6,
        chargeMassMax=3,
        grainPR=0,
        grainLDR=2.5,
        startPressure=30e6,
        chamberExpansion=1.1,
        designedPress=300e6,
        designedVel=1500,
        lengthGunMinCal=0,
        lengthGunMaxCal=500,
        loadFractionMin=0.05,
        loadFractionMax=0.65,
        grainArcMin=0.1e-3,
        grainArcMax=3e-3,
        tol=1e-3,
        arcTol=1e-6,  # 1um
        xstep=2,
        ystep=2,
    )
