from opt import Constrained
from prop import GrainComp, SimpleGeometry
from prop import Propellant
from gun import SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV
from gun import POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT

compositions = GrainComp.readFile("data/propellants.csv")
M10 = compositions["M10"]  # standin for pyroxylin
# beta does not affect the actual curve....
M10Cy = Propellant(M10, SimpleGeometry.TUBE, 1, 100)

tol = 1e-3
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
        chamberVolume = D81.m * chargeMassRatio / D81.propellant.rho_p
        tubeVolume = (lengthGun) * D81.S
        volume = chamberVolume + tubeVolume  # convert to liters

        halfWeb *= 2e3
        volume *= 1e3
    except ValueError:
        halfWeb, lengthGun, volume = None, None, None

    return loadFraction, chargeMassRatio, halfWeb, lengthGun, volume


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import multiprocessing

    parameters = []
    for i in range(63):
        chargeMassRatio = 0.9 + 0.05 * i
        for j in range(13):
            loadFraction = 0.2 + 0.05 * j
            parameters.append((loadFraction, chargeMassRatio))

    with multiprocessing.Pool() as pool:
        results = pool.starmap(f, parameters)

    xs, ys, es, ls, vs = zip(*[v for v in results if v[2] is not None])

    fig, ax = plt.subplots()

    def makeLevels(data, delta):
        return [
            (min(data) // delta) * delta + v * delta
            for v in range(int((max(data) - min(data)) // delta))
        ]

    print(max(ls))

    gunLengthContours = ax.tricontour(
        xs,
        ys,
        ls,
        levels=[5, 6, 7, 8, 9],
        colors="black",
        alpha=0.5,
        linewidths=1,
    )
    ax.clabel(
        gunLengthContours,
        gunLengthContours.levels,
        inline=True,
        fontsize=8,
        fmt=lambda x: f"{x:.1f} m",
    )

    print(max(vs))

    gunVolumeContours = ax.tricontour(
        xs,
        ys,
        vs,
        levels=[60, 70, 80, 90, 100],
        colors="blue",
        alpha=0.5,
        linestyles="dotted",
        linewidths=1,
    )
    ax.clabel(
        gunVolumeContours,
        gunVolumeContours.levels,
        inline=True,
        fontsize=8,
        fmt=lambda x: f"{x:.0f} L",
    )

    propArchContours = ax.tricontour(
        xs,
        ys,
        es,
        levels=[0.4 + 0.2 * v for v in range(9)],
        colors="red",
        alpha=0.5,
        linewidths=1,
    )
    ax.clabel(
        propArchContours,
        propArchContours.levels,
        inline=True,
        fontsize=8,
        fmt=lambda x: f"{x:.1f} mm",
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 4)
    plt.show()
