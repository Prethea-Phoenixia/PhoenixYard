equilibriumKT = [
    [800.0, 0.2478, 0, 0, 0, 0, 0, 0, 31.25, 2.91e-3],
    [1000.0, 0.7286, 0, 0, 0, 0, 0, 0, 3.72e-2, 5.64e-4],
    [1200.0, 1.435, 0, 0, 0, 0, 0, 0, 3.99e-4, 1.85e-4],
    [1400.0, 2.270, 0, 0, 0, 0, 0, 0, 1.55e-5, 8.2e-5],
    [1500.0, 2.704, 0, 0, 0, 0, 0, 0, 4.2e-6, 5.9e-5],
    [1600.0, 3.132, 6e-5, 0, 0, 0, 0, 0, 1.35e-6, 4.5e-5],
    [1700.0, 3.555, 2e-4, 0, 0, 0, 0, 0, 4.9e-7, 3.5e-5],
    [1800.0, 3.975, 4e-4, 2e-5, 0, 0, 0, 0, 2.0e-7, 2.8e-5],
    [1900.0, 4.385, 9e-4, 6e-5, 0, 0, 0, 0, 9.2e-8, 2.3e-5],
    [2000.0, 4.782, 0.0016, 2e-4, 0, 0, 0, 0, 4.5e-8, 1.9e-5],
    [2100.0, 5.161, 0.0031, 4e-4, 2e-5, 0, 0, 0, 2.4e-8, 1.6e-5],
    [2200.0, 5.520, 0.0057, 8e-4, 5e-5, 0, 0, 0, 1.3e-8, 1.4e-5],
    [2300, 5.852, 0.0097, 0.0016, 1e-4, 1e-5, 1e-5, 0, 0, 1.2e-5],
    [2400.0, 6.155, 0.0159, 0.0029, 2e-4, 3e-5, 3e-5, 1e-5, 0, 1.1e-5],
    [2500.0, 6.433, 0.0251, 0.0052, 4e-4, 7e-5, 9e-5, 3e-5, 0, 1.0e-5],
    [2600.0, 6.694, 0.0383, 0.0090, 7e-4, 2e-4, 2e-4, 1e-4, 0, 0.9e-5],
    [2700.0, 6.939, 0.0566, 0.0146, 0.0012, 3e-4, 5e-4, 2e-4, 0, 0.8e-5],
    [2800.0, 7.167, 0.0814, 0.0231, 0.0020, 5e-4, 0.0012, 5e-4, 0, 0.8e-5],
    [2900.0, 7.379, 0.1143, 0.0355, 0.0034, 8e-4, 0.0026, 0.0010, 0, 0.7e-5],
    [3000.0, 7.574, 0.1574, 0.0529, 0.0055, 0.0014, 0.0053, 0.0020, 0, 0.6e-5],
    [3100.0, 7.753, 0.2125, 0.0768, 0.0086, 0.0022, 0.0103, 0.0038, 0, 0],
    [3200.0, 7.917, 0.2813, 0.1089, 0.0129, 0.0035, 0.0190, 0.0071, 0, 0],
    [3300.0, 8.068, 0.3658, 0.1513, 0.0188, 0.0053, 0.0339, 0.0126, 0, 0],
    [3400.0, 8.205, 0.4682, 0.2064, 0.0271, 0.0078, 0.0586, 0.0216, 0, 0],
    [3500.0, 8.330, 0.5910, 0.2760, 0.0387, 0.0113, 0.0978, 0.0358, 0, 0],
    [3600.0, 8.443, 0.7367, 0.3626, 0.0539, 0.0161, 0.1591, 0.0580, 0, 0],
    [3700.0, 8.544, 0.9079, 0.4693, 0.0734, 0.0225, 0.2516, 0.0917, 0, 0],
    [3800.0, 8.633, 1.1070, 0.6001, 0.0982, 0.0309, 0.3886, 0.1412, 0, 0],
    [3900.0, 8.710, 1.3360, 0.7589, 0.1300, 0.0417, 0.5878, 0.2132, 0, 0],
    [4000.0, 8.775, 1.5970, 0.9495, 0.1693, 0.0554, 0.8711, 0.3157, 0, 0],
]
"""
Kp1 = pO / pO2^1/2
Kp2 = pH / pH21/2
Kp3 = pOH / (pO2^1/2 pH2^1/2)
Kp4 = pH2O / (pH2 pO2^1/2)
Kp5 = pN / pN2^1/2
Kp6 = pNO / (pN2 pO2^1/2)
Kp7 = pC(c,diamond) / pC(c,graphite)
Kp8 = pC (g) / pC
Kp9 = pCO / (pC pO2^1/2)
Kp10 = pCO2 / (pC pO2)
Kp11 = pCH4 / (pC (pH2)^2)
Kp12 = pC2H2 / ((pC)^2 pH2)
Kp13 = pCl / pCl21/2
Kp14 = pHCl / (pH2^1/2 pCl2^1/2)
Kp15 = pNH / (pN2^1/2 pH2^1/2)
Kp16 = pNH3 / (pN2^1/2 pH2^3/2)
Kp17 = pCH2 / (pC pH2)
Kp18 = pNO2 / (pNO pO2^1/2)
Kp19 = pNO / ( pN2^1/2 pO2^1/2)
"""

from nasa7 import Specie, Reaction

Specie.read("data/nasa7.dat")

CO2 = Specie.get("CO2")
H2O = Specie.get("H2O")
CO = Specie.get("CO")
H2 = Specie.get("H2")
H = Specie.get("H")
N2 = Specie.get("N2")
OH = Specie.get("OH")
NO = Specie.get("NO")
N = Specie.get("N")
O = Specie.get("O")
O2 = Specie.get("O2")
CH4 = Specie.get("CH4")
NH3 = Specie.get("NH3")

k0 = Reaction("water-gas", LHS={CO2: 1, H2: 1}, RHS={CO: 1, H2O: 1})
k1 = Reaction("hydrogen decomposition", LHS={H2: 0.5}, RHS={H: 1})
k2 = Reaction("water-hydroxyl", LHS={H2O: 1}, RHS={H2: 0.5, OH: 1})
k3 = Reaction("water-nitroxide", LHS={H2O: 1, N2: 0.5}, RHS={NO: 1, H2: 1})
k4 = Reaction("nitrogen decomposiiton", LHS={N2: 0.5}, RHS={N: 1})
k5 = Reaction("water-oxygen radical", LHS={H2O: 1}, RHS={O: 1, H2: 1})
k6 = Reaction("water-decomposition", LHS={H2O: 1}, RHS={O2: 0.5, H2: 1})
k7 = Reaction(
    "methane synthesis w/ carbon dioxide",
    LHS={CO: 2, H2: 2},
    RHS={CH4: 1, CO2: 1},
)
k8 = Reaction("ammonia synthesis", LHS={N2: 0.5, H2: 1.5}, RHS={NH3: 1})

if __name__ == "__main__":
    import numpy as np
    from matplotlib import pyplot as plt

    KTTransposed = list(zip(*[line for line in equilibriumKT if line[8] != 0]))
    fit7 = np.polynomial.polynomial.Polynomial.fit(
        np.power(KTTransposed[0], -1), np.log10(KTTransposed[8]), 5
    )  # this is an inverse fit!

    fig, axs = plt.subplots(3, 3)

    # these are the 1500-5000K curves
    kp1_h = np.polynomial.polynomial.Polynomial(
        [0, 1.8636e-10, -2.3475e-06, 1.0500e-02, -1.6377e01][::-1]
    )  # 1500-5000K

    kp2_h = np.polynomial.polynomial.Polynomial(
        [0, 1.6505e-10, -2.0819e-06, 9.4000e-03, -1.4586e01][::-1]
    )  # 1500-5000K

    kp3_h = np.polynomial.polynomial.Polynomial(
        [0, 3.0374e-11, -3.8141e-07, 1.7000e-03, -2.4360e00][::-1]
    )  # 1500-5000K

    kp4_l = np.polynomial.polynomial.Polynomial(
        [5.1208e-11, -2.2402e-07, 3.6861e-04, -2.8230e-01, 9.6794e01][::-1]
    )  # 1500-5000K
    kp4_h = np.polynomial.polynomial.Polynomial(
        [0, -1.8559e-10, 2.3375e-06, -1.0500e-02, 1.6715e01][::-1]
    )  # 1500-5000K

    kp5_h = np.polynomial.polynomial.Polynomial(
        [0, 2.6841e-10, -3.3790e-06, 1.5200e-02, -2.5117e01][::-1]
    )  # 1500-5000K

    kp6_h = np.polynomial.polynomial.Polynomial(
        [0, 6.6881e-11, -8.4371e-07, 3.8000e-03, -6.4496e00][::-1]
    )  # 1500-5000K

    kp9_l = np.polynomial.polynomial.Polynomial(
        [2.3424e-11, -1.0237e-07, 1.6815e-04, -1.2840e-01, 4.9814e01][::-1]
    )  # 298 K - 1500K
    kp9_h = np.polynomial.polynomial.Polynomial(
        [0, -1.3698e-10, 1.4770e-06, -5.9000e-03, 1.4414e01][::-1]
    )  # 1500K - 4000K

    kp10_l = np.polynomial.polynomial.Polynomial(
        [8.3822e-11, -3.6648e-07, 6.0237e-04, -4.6010e-01, 1.6147e02][::-1]
    )  # 298K - 1500K
    kp10_h = np.polynomial.polynomial.Polynomial(
        [0, -4.8969e-10, 5.2746e-06, -2.0600e-02, 3.4372e01][::-1]
    )  # 1500K - 4000K

    kp11 = np.polynomial.polynomial.Polynomial(
        [1.5091e-11, -6.6349e-08, 1.1032e-04, -8.6400e-02, 2.6451e01][::-1]
    )  # 298-1500K

    kp16_l = np.polynomial.polynomial.Polynomial(
        [1.0847e-12, -7.0150e-09, 1.7592e-05, -2.1400e-02, 6.4824e00][::-1]
    )  # 600K- 2000K
    kp16_h = np.polynomial.polynomial.Polynomial(
        [
            0,
            -4.4574e-11,
            5.2571e-07,
            -2.2000e-03,
            -1.9461e00 - 0.1145479999999992,
        ][::-1]
    )  # 2000K - 4000K

    for i in range(9):
        ax = axs[i // 3][i % 3]
        x, y = list(
            zip(*[(v[0], v[i + 1]) for v in equilibriumKT if v[i + 1] != 0])
        )  # transpose the data
        ax.plot(x, y, ".")

        f1, f2 = None, None
        mint, transt, maxt = 800, 1500, 4000
        invFit1, invFit2 = False, False
        if i == 0:
            # k0 = kp9 * kp4/ kp10
            f1 = np.polynomial.polynomial.Polynomial(
                kp9_l.coef + kp4_l.coef - kp10_l.coef
            )
            f2 = np.polynomial.polynomial.Polynomial(
                kp9_h.coef + kp4_h.coef - kp10_h.coef
            )
        elif i == 1:
            f2 = kp2_h
        elif i == 2:
            # k2 = kp3/kp4
            f2 = np.polynomial.polynomial.Polynomial(kp3_h.coef - kp4_h.coef)
        elif i == 3:
            # k3 = kp6/kp4
            f2 = np.polynomial.polynomial.Polynomial(kp6_h.coef - kp4_h.coef)
        elif i == 4:
            f2 = kp5_h
        elif i == 5:
            # k5 = kp1/kp4:
            f2 = np.polynomial.polynomial.Polynomial(kp1_h.coef - kp4_h.coef)
        elif i == 6:
            f2 = np.polynomial.polynomial.Polynomial(-2 * kp4_h.coef)

        elif i == 8:
            mint = 600
            f1 = kp16_l
            transt = 2000
            f2 = kp16_h
        """
        elif i == 7:
            # k7 = kp10 * kp11 / kp9^2
            f1 = np.polynomial.polynomial.Polynomial(
                kp10_l.coef + kp11.coef - 2 * kp9_l.coef
            )
            f2 = fit7
            invFit2 = True
        """

        if f1 is not None:
            t1 = np.linspace(mint, transt, 100)
            ax.plot(t1, 10 ** f1(1 / t1) if invFit1 else 10 ** f1(t1))
        if f2 is not None:
            t2 = np.linspace(transt, maxt, 100)
            ax.plot(t2, 10 ** f2(1 / t2) if invFit2 else 10 ** f2(t2))

        ax.set_yscale("log")

    for i, k in enumerate((k0, k1, k2, k3, k4, k5, k6, k7, k8)):
        t = np.linspace(800, 4000, 100)
        ax = axs[i // 3][i % 3]
        ax.plot(t, [k(v) for v in t], c="grey", ls="dashed")

    plt.show()
    input()
    ###########################################################################

    newK = []

    logK_l = [
        np.polynomial.polynomial.Polynomial(kp9_l.coef + kp4_l.coef - kp10_l.coef),
        None,
        None,
        None,
        None,
        None,
        None,
        np.polynomial.polynomial.Polynomial(kp10_l.coef + kp11.coef - 2 * kp9_l.coef),
        kp16_l,
    ]

    isInvFit_l = [False] * 9
    isInvFit_h = [False] * 7 + [True, False]

    logK_h = [
        np.polynomial.polynomial.Polynomial(kp9_h.coef + kp4_h.coef - kp10_h.coef),
        kp2_h,
        np.polynomial.polynomial.Polynomial(kp3_h.coef - kp4_h.coef),
        np.polynomial.polynomial.Polynomial(kp6_h.coef - kp4_h.coef),
        kp5_h,
        np.polynomial.polynomial.Polynomial(kp1_h.coef - kp4_h.coef),
        np.polynomial.polynomial.Polynomial(-2 * kp4_h.coef),
        fit7,
        kp16_h,
    ]

    transT = [1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 2000]
    f = None
    isInv = None
    for line in equilibriumKT:
        T, Ks = line[0], line[1:]
        if 4000 >= T > 1500:
            newLine = [T]
            for i in range(9):
                if Ks[i] == 0:
                    if T >= transT[i]:
                        f = logK_h[i]
                        isInv = isInvFit_h[i]
                    else:
                        f = logK_l[i]
                        isInv = isInvFit_l[i]
                    if f is not None:
                        newLine.append(
                            float(
                                "{:.4g}".format(10 ** f(1 / T) if isInv else 10 ** f(T))
                            )
                        )
                    else:
                        newLine.append(None)
                else:
                    newLine.append(Ks[i])

            newK.append(newLine)

    from tabulate import tabulate

    print("oldK")
    print(
        tabulate(
            equilibriumKT,
            headers=["T"] + ["K{:}".format(i) for i in range(9)],
        ),
    )
    print("newK")
    print(
        tabulate(
            newK,
            headers=["T"] + ["K{:}".format(i) for i in range(9)],
        ),
    )
    print(newK)

    # fmt:off
    newK = [
            [1600.0, 3.132, 6e-05, 1.964e-06, 4.033e-08, 4.489e-09, 1.091e-10, 5.277e-11, 1.35e-06, 4.5e-05],
            [1700.0, 3.555, 0.0002, 6.2e-06, 1.556e-07, 1.89e-08, 7.866e-10, 3.827e-10, 4.9e-07, 3.5e-05],
            [1800.0, 3.975, 0.0004, 2e-05, 5.504e-07, 7.251e-08, 4.987e-09, 2.442e-09, 2e-07, 2.8e-05],
            [1900.0, 4.385, 0.0009, 6e-05, 1.79e-06, 2.546e-07, 2.795e-08, 1.378e-08, 9.2e-08, 2.3e-05],
            [2000.0, 4.782, 0.0016, 0.0002, 5.374e-06, 8.209e-07, 1.392e-07, 6.909e-08, 4.5e-08, 1.9e-05],
            [2100.0, 5.161, 0.0031, 0.0004, 2e-05, 2.44e-06, 6.191e-07, 3.096e-07, 2.4e-08, 1.6e-05],
            [2200.0, 5.52, 0.0057, 0.0008, 5e-05, 6.709e-06, 2.472e-06, 1.245e-06, 1.3e-08, 1.4e-05],
            [2300, 5.852, 0.0097, 0.0016, 0.0001, 1e-05, 1e-05, 4.523e-06, 7.725e-09, 1.2e-05],
            [2400.0, 6.155, 0.0159, 0.0029, 0.0002, 3e-05, 3e-05, 1e-05, 4.745e-09, 1.1e-05],
            [2500.0, 6.433, 0.0251, 0.0052, 0.0004, 7e-05, 9e-05, 3e-05, 3.032e-09, 1e-05],
            [2600.0, 6.694, 0.0383, 0.009, 0.0007, 0.0002, 0.0002, 0.0001, 2.007e-09, 9e-06],
            [2700.0, 6.939, 0.0566, 0.0146, 0.0012, 0.0003, 0.0005, 0.0002, 1.369e-09, 8e-06],
            [2800.0, 7.167, 0.0814, 0.0231, 0.002, 0.0005, 0.0012, 0.0005, 9.608e-10, 8e-06],
            [2900.0, 7.379, 0.1143, 0.0355, 0.0034, 0.0008, 0.0026, 0.001, 6.909e-10, 7.825e-06],
            [3000.0, 7.574, 0.1574, 0.0529, 0.0055, 0.0014, 0.0053, 0.002, 5.079e-10, 7.366e-06],
            [3100.0, 7.753, 0.2125, 0.0768, 0.0086, 0.0022, 0.0103, 0.0038, 3.809e-10, 6.975e-06],
            [3200.0, 7.917, 0.2813, 0.1089, 0.0129, 0.0035, 0.019, 0.0071, 2.909e-10, 6.638e-06],
            [3300.0, 8.068, 0.3658, 0.1513, 0.0188, 0.0053, 0.0339, 0.0126, 2.258e-10, 6.346e-06],
            [3400.0, 8.205, 0.4682, 0.2064, 0.0271, 0.0078, 0.0586, 0.0216, 1.779e-10, 6.09e-06],
            [3500.0, 8.33, 0.591, 0.276, 0.0387, 0.0113, 0.0978, 0.0358, 1.421e-10, 5.864e-06],
            [3600.0, 8.443, 0.7367, 0.3626, 0.0539, 0.0161, 0.1591, 0.058, 1.149e-10, 5.661e-06],
            [3700.0, 8.544, 0.9079, 0.4693, 0.0734, 0.0225, 0.2516, 0.0917, 9.398e-11, 5.477e-06],
            [3800.0, 8.633, 1.107, 0.6001, 0.0982, 0.0309, 0.3886, 0.1412, 7.769e-11, 5.306e-06],
            [3900.0, 8.71, 1.336, 0.7589, 0.13, 0.0417, 0.5878, 0.2132, 6.485e-11, 5.144e-06],
            [4000.0, 8.775, 1.597, 0.9495, 0.1693, 0.0554, 0.8711, 0.3157, 5.462e-11, 4.989e-06]
    ]
