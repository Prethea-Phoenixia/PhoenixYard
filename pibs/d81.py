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

from labellines import labelLine, labelLines

compositions = GrainComp.readFile("data/propellants.csv")
M10 = compositions["M10"]  # standin for pyroxylin
# beta does not affect the actual curve....
M10Cy = Propellant(M10, SimpleGeometry.TUBE, 1, 100)

caliber = 125e-3
tol = 1e-3
dragCoefficient = 3e-2
D81 = Constrained(
    caliber=caliber,
    shotMass=5.67,
    propellant=M10Cy,
    startPressure=30e6,
    dragCoefficient=dragCoefficient,
    designPressure=392e6,
    designVelocity=1800,
    chambrage=1.25,  # almost exact, 4 * 12.27L/ pi (1.25dm)^2 * 8dm
)


def f(loadFraction, chargeMassRatio):
    chargeMass = D81.m * chargeMassRatio
    loadDensity = loadFraction * D81.propellant.rho_p
    try:
        halfWeb, lengthGun = D81.solve(
            loadFraction,
            chargeMassRatio,
            tol,
            minWeb=1e-4,
            maxLength=1e2,
            sol=SOL_LAGRANGE,
            ambientRho=1.204,
            ambientP=101.325e3,
            ambientGamma=1.4,
            control=POINT_PEAK_AVG,
        )

        chamberVolume = chargeMass / loadDensity

        gun = Gun(
            caliber=caliber,
            shotMass=D81.m,
            propellant=M10Cy,
            grainSize=2 * halfWeb,
            chargeMass=chargeMass,
            chamberVolume=chamberVolume,
            startPressure=D81.p_0,
            lengthGun=lengthGun,
            chambrage=D81.chi_k,
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

        tubeVolume = (lengthGun) * D81.S
        volume = chamberVolume + tubeVolume  # convert to liters

        halfWeb *= 2e3
        volume *= 1e3
    except ValueError:
        # print(e)
        halfWeb, lengthGun, volume, burnout = None, None, None, None

    return loadDensity, chargeMass, halfWeb, lengthGun, burnout, volume


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import multiprocessing

    parameters = []
    for i in range(int(3.1 / 0.05) + 1):
        chargeMassRatio = 0.9 + 0.05 * i
        for j in range(int(0.55 / 0.05) + 1):
            loadFraction = 0.2 + 0.05 * j
            parameters.append((loadFraction, chargeMassRatio))

    with multiprocessing.Pool() as pool:
        results = pool.starmap(f, parameters)

    # print(*results, sep="\n")

    xs, ys, es, ls, bos, vs = zip(*[v for v in results if v[2] is not None])

    fig, ax = plt.subplots(layout="constrained")

    for chamberVolume in [5, 10, 15, 20, 25, 30, 35]:
        ax.plot(
            (0, D81.propellant.rho_p),
            (0, chamberVolume * 1e-3 * D81.propellant.rho_p),
            c="grey",
            alpha=0.5,
            linestyle="dashed",
            label=f"{chamberVolume:.0f} L",
            linewidth=1,
        )

    labelLines(ax.get_lines(), xvals=(0, 600), fontsize=8, outline_width=1)

    def makeLevels(data, delta):
        return [
            (min(data) // delta) * delta + v * delta
            for v in range(int((max(data) - min(data)) // delta))
        ]

    gunLengthContours = ax.tricontour(
        xs,
        ys,
        ls,
        levels=[4 + 0.2 * v for v in range(10)],
        linestyles=["solid"] + ["dotted" for _ in range(4)],
        colors="black",
        alpha=0.5,
        linewidths=1,
    )

    clabels = ax.clabel(
        gunLengthContours,
        gunLengthContours.levels[:10:5],
        inline=True,
        fontsize=8,
        fmt=lambda x: f"{x:.1f} m",
        use_clabeltext=True,
    )

    gunLengthContours = ax.tricontour(
        xs,
        ys,
        ls,
        levels=[6, 7, 8, 9],
        linestyles="solid",
        colors="black",
        alpha=0.5,
        linewidths=1,
    )

    clabels = ax.clabel(
        gunLengthContours,
        gunLengthContours.levels,
        inline=True,
        fontsize=8,
        fmt=lambda x: f"{x:.1f} m",
        use_clabeltext=True,
    )

    print(max(vs))
    """
    gunVolumeContours = ax.tricontour(
        xs,
        ys,
        vs,
        levels=[70, 80, 90, 100],
        colors="blue",
        alpha=0.5,
        linestyles=["solid" for _ in range(5)],
        linewidths=1,
    )
    ax.clabel(
        gunVolumeContours,
        gunVolumeContours.levels,
        inline=True,
        fontsize=8,
        fmt=lambda x: f"{x:.0f} L",
    )
    """

    propArchContours = ax.tricontour(
        xs,
        ys,
        es,
        levels=[0.2 * v for v in range(31)],
        colors="red",
        linestyles=["solid"] + ["dotted" for _ in range(4)],
        alpha=0.5,
        linewidths=1,
    )
    clabels = ax.clabel(
        propArchContours,
        propArchContours.levels[::5],
        inline=True,
        fontsize=8,
        fmt=lambda x: f"{x:.0f} mm",
        use_clabeltext=True,
    )

    gunBurnoutContours = ax.tricontour(
        xs,
        ys,
        bos,
        levels=[0.2 * v for v in range(5)],
        colors="purple",
        alpha=0.5,
        linewidths=1,
    )

    clabels = ax.clabel(
        gunBurnoutContours,
        gunBurnoutContours.levels,
        inline=True,
        fontsize=8,
        fmt=lambda x: f"{x:.0%}",
        use_clabeltext=True,
    )
    """
    for label in clabels:
        label.set_va("top")
        # label.set_rotation(0)
    """
    ax.set_xlim(0, 0.8 * D81.propellant.rho_p)
    ax.set_ylim(0, 4 * D81.m)

    plt.show()
