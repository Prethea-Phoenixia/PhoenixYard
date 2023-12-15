from opt import Constrained
from prop import GrainComp, SimpleGeometry
from prop import Propellant
from gun import SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV
from gun import POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT

from math import atan2, cos, sin
from labellines import labelLine, labelLines

compositions = GrainComp.readFile("data/propellants.csv")
M10 = compositions["M10"]  # standin for pyroxylin
# beta does not affect the actual curve....
M10Cy = Propellant(M10, SimpleGeometry.TUBE, 1, 100)

tol = 1e-5
D81 = Constrained(
    caliber=125e-3,
    shotMass=5.67,
    propellant=M10Cy,
    startPressure=30e6,
    dragCoefficient=3e-2,
    designPressure=392e6,
    designVelocity=1800,
    chambrage=1.25,  # almost exact, 4 * 12.27L/ pi (1.25dm)^2 * 8dm
)


def f(loadFraction, chargeMassRatio):
    charge = D81.m * chargeMassRatio
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

        chamberVolume = charge / D81.propellant.rho_p
        tubeVolume = (lengthGun) * D81.S
        volume = chamberVolume + tubeVolume  # convert to liters

        halfWeb *= 2e3
        volume *= 1e3
    except ValueError:
        halfWeb, lengthGun, volume = None, None, None

    return loadFraction, charge, halfWeb, lengthGun, volume


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import multiprocessing

    parameters = []
    for i in range(int(3.1 / 0.05) + 1):
        chargeMassRatio = 0.9 + 0.05 * i
        for j in range(int(0.5 / 0.05) + 1):
            loadFraction = 0.2 + 0.05 * j
            parameters.append((loadFraction, chargeMassRatio))

    with multiprocessing.Pool() as pool:
        results = pool.starmap(f, parameters)

    xs, ys, es, ls, vs = zip(*[v for v in results if v[2] is not None])

    fig, ax = plt.subplots(layout="constrained")

    for chamberVolume in [5, 10, 15, 20, 25, 30, 35]:
        ax.plot(
            (0, 1),
            (0, chamberVolume * 1e-3 * D81.propellant.rho_p),
            c="grey",
            alpha=0.5,
            linestyle="dashed",
            label=f"{chamberVolume:.0f} L",
            linewidth=1,
        )

    labelLines(ax.get_lines(), xvals=(0, 0.2), fontsize=8, outline_width=1)

    def makeLevels(data, delta):
        return [
            (min(data) // delta) * delta + v * delta
            for v in range(int((max(data) - min(data)) // delta))
        ]

    gunLengthContours = ax.tricontour(
        xs,
        ys,
        ls,
        levels=[4 + 0.2 * v for v in range(10)] + [6, 7, 8, 9],
        linestyles=["solid"]
        + ["dashed" for _ in range(4)]
        + ["solid"]
        + ["dashed" for _ in range(4)]
        + ["solid" for _ in range(5)],
        colors="black",
        alpha=0.5,
        linewidths=1,
    )

    clabels = ax.clabel(
        gunLengthContours,
        gunLengthContours.levels,
        inline=False,
        fontsize=8,
        fmt=lambda x: f"{x:.1f} m",
        inline_spacing=0,
    )

    for label in clabels:
        label.set_va("top")

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
        propArchContours.levels,
        inline=False,
        fontsize=8,
        fmt=lambda x: f"{x:.1f} mm",
        inline_spacing=0,
    )

    for label in clabels:
        label.set_va("top")

    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, 4 * D81.m)

    plt.show()
