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
PyTu = Propellant(Py, SimpleGeometry.TUBE, 1, 100, 0.225)

caliber = 125e-3
tol = 1e-3
dragCoefficient = 2e-2
control = POINT_PEAK_AVG

target = Constrained(
    caliber=caliber,
    shotMass=5.67,
    propellant=PyTu,
    startPressure=30e6,
    dragCoefficient=dragCoefficient,
    designPressure=392e6,
    designVelocity=1800,
    chambrage=1.25,  # almost exact, 4 * 12.27L/ pi (1.25dm)^2 * 8dm
)


def f(loadFraction, chargeMassRatio):
    chargeMass = target.m * chargeMassRatio
    loadDensity = loadFraction * target.propellant.rho_p
    try:
        halfWeb, lengthGun = target.solve(
            loadFraction,
            chargeMassRatio,
            tol,
            minWeb=1e-6,
            maxLength=1e3,
            sol=SOL_LAGRANGE,
            ambientRho=1.204,
            ambientP=101.325e3,
            ambientGamma=1.4,
            control=control,
        )

        chamberVolume = chargeMass / loadDensity

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

        tubeVolume = (lengthGun) * target.S
        volume = chamberVolume + tubeVolume  # convert to liters

        halfWeb *= 2e3
        volume *= 1e3
        lengthGun *= 10
    except ValueError:
        # print(e)
        halfWeb, lengthGun, volume, burnout = None, None, None, None

    return loadDensity, chargeMass, halfWeb, lengthGun, burnout, volume


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import multiprocessing

    parameters = []
    for i in range(int(3.1 / 0.1) + 1):
        chargeMassRatio = 0.9 + 0.1 * i
        for j in range(int(0.8 / 0.025) + 1):
            loadFraction = 0.1 + 0.025 * j
            parameters.append((loadFraction, chargeMassRatio))

    with multiprocessing.Pool() as pool:
        results = pool.starmap(f, parameters)

    # print(*results, sep="\n")

    xs, ys, es, ls, bos, vs = zip(*[v for v in results if v[2] is not None])

    fig, ax = plt.subplots(layout="constrained", figsize=(11.69, 8.27))

    for chamberVolume in [5, 10, 12.27, 15, 20, 25, 30, 35]:
        ax.plot(
            (0, target.propellant.rho_p),
            (0, chamberVolume * 1e-3 * target.propellant.rho_p),
            c="black",
            alpha=0.5,
            linestyle="dotted",
            label=(
                f"{chamberVolume:.0f} L"
                if chamberVolume != 12.27
                else f"{chamberVolume:.2f} L"
            ),
            linewidth=1,
        )

    labelLines(ax.get_lines(), xvals=(0, 600), fontsize=10, outline_width=1)

    def makeLevels(data, delta):
        return [
            (min(data) // delta) * delta + v * delta
            for v in range(int((max(data) - min(data)) // delta))
        ]

    gunLengthContours = ax.tricontour(
        xs,
        ys,
        ls,
        levels=[48 + 4 * v for v in range(10)],
        linestyles="solid",
        colors="black",
        alpha=1,
        linewidths=1,
    )

    clabels = ax.clabel(
        gunLengthContours,
        gunLengthContours.levels,
        inline=True,
        fontsize=10,
        use_clabeltext=True,
    )
    """
    gunVolumeContours = ax.tricontour(
        xs,
        ys,
        vs,
        levels=[70, 76.1, 80, 90, 100, 110],
        colors=["black", "red", "black", "black", "black", "black", "black"],
        alpha=1,
        linestyles="solid",
        linewidths=1,
    )
    clabels = ax.clabel(
        gunVolumeContours,
        gunVolumeContours.levels,
        inline=False,
        fontsize=8,
        fmt=lambda x: f"{x:.0f} L",
        use_clabeltext=True,
    )"""

    propArchContours = ax.tricontour(
        xs,
        ys,
        es,
        levels=[0, 0.5, 1, 1.5, 2, 2.5, 3],
        colors="black",
        linestyles="dashed",
        alpha=1,
        linewidths=1,
    )
    clabels = ax.clabel(
        propArchContours,
        propArchContours.levels,
        inline=True,
        fontsize=10,
        use_clabeltext=True,
    )
    gunBurnoutContours = ax.tricontour(
        xs,
        ys,
        bos,
        levels=[0.8],
        colors="black",
        linestyles="dashdot",
        alpha=1,
        # linewidths=1,
    )

    clabels = ax.clabel(
        gunBurnoutContours,
        gunBurnoutContours.levels,
        inline=False,
        fontsize=10,
        fmt=lambda x: f"{x:.0%}",
        use_clabeltext=True,
    )
    for label in clabels:
        label.set_va("bottom")

    ax.tricontourf(
        xs,
        ys,
        bos,
        levels=[0.0, 0.8],
        colors="black",
        alpha=0.1,
    )

    ax.scatter(10 / 12.27 * 1e3, 10, s=100, marker="*")
    ax.set_xlim(0, 0.8 * target.propellant.rho_p)
    ax.set_ylim(0, 4 * target.m)

    plt.show()
