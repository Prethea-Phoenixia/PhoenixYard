from math import log10, exp
from num import quadratic

"""
TABLE 2.07 from Hunt
T in Kelvin, -DeltaB in cc/(gm.mol), -DeltaC/2 in (cc/(gm.mol))^2
"""

negDeltaTable = [
    [1600, 34.2, 490],
    [1700, 33.8, 460],
    [1800, 33.5, 435],
    [1900, 33.2, 410],
    [2000, 33.0, 390],
    [2100, 32.8, 370],
    [2200, 32.6, 355],
    [2300, 32.5, 340],
    [2400, 32.4, 325],
    [2500, 32.2, 310],
    [2600, 32.1, 300],
    [2700, 32.0, 290],
    [2800, 31.9, 280],
    [2900, 31.8, 270],
    [3000, 31.7, 260],
    [3100, 31.6, 255],
    [3200, 31.5, 245],
    [3300, 31.4, 235],
    [3400, 31.4, 230],
    [3500, 31.3, 225],
    [3600, 31.2, 215],
    [3700, 31.1, 210],
    [3800, 31.1, 205],
    [3900, 31.0, 200],
    [4000, 30.9, 195],
]
"""
TABLE 2.04 from Hunt
T in Kelvin, K0(T)..... K6(T)
"""
equilibriumKT = [
    [1000, 0.7185, 3e-11, 1e-14, 8e-21, 2e-10, 3e-9, 5e-16],
    [1200, 1.406, 7e-9, 1e-11, 2e-16, 3e-8, 2e-7, 7e-13],
    [1400, 2.212, 3e-7, 1e-9, 2e-13, 1e-6, 5e-6, 1e-10],
    [1600, 3.043, 6e-6, 5e-8, 5e-11, 2e-5, 6e-5, 5e-9],
    [1700, 3.438, 2e-5, 2e-7, 5e-10, 6e-5, 0.0002, 2e-8],
    [1800, 3.832, 6e-5, 8e-7, 3e-9, 0.0001, 0.0004, 1e-7],
    [1900, 4.206, 0.0002, 3e-6, 2e-8, 0.0003, 0.0008, 4e-7],
    [2000, 4.574, 0.0004, 8e-6, 9e-8, 0.0007, 0.0017, 1e-6],
    [2100, 4.909, 0.0009, 2e-5, 4e-7, 0.0015, 0.0033, 3e-6],
    [2200, 5.235, 0.0018, 5e-5, 1e-6, 0.0029, 0.0060, 8e-6],
    [2300, 5.533, 0.0035, 0.0001, 5e-6, 0.0053, 0.0103, 2e-5],
    [2400, 5.821, 0.0063, 0.0002, 1e-5, 0.0091, 0.0168, 4e-5],
    [2500, 6.081, 0.0109, 0.0004, 5e-5, 0.0152, 0.0265, 9e-5],
    [2600, 6.331, 0.0182, 0.0009, 0.0001, 0.0244, 0.0403, 0.0002],
    [2700, 6.552, 0.0292, 0.0016, 0.0002, 0.0378, 0.0596, 0.0004],
    [2800, 6.764, 0.0453, 0.0026, 0.0005, 0.0567, 0.0856, 0.0006],
    [2900, 6.949, 0.0680, 0.0043, 0.0011, 0.0828, 0.1201, 0.0010],
    [3000, 7.127, 0.0995, 0.0070, 0.0023, 0.1179, 0.1648, 0.0017],
    [3100, 7.281, 0.1419, 0.0108, 0.0044, 0.1639, 0.2217, 0.0027],
    [3200, 7.428, 0.1981, 0.0163, 0.0081, 0.2233, 0.2926, 0.0042],
    [3300, 7.552, 0.2709, 0.0240, 0.0144, 0.2989, 0.3799, 0.0063],
    [3400, 7.670, 0.3638, 0.0345, 0.0248, 0.3931, 0.4858, 0.0092],
    [3500, 7.764, 0.4801, 0.0486, 0.0415, 0.5094, 0.6130, 0.0133],
    [3600, 7.854, 0.6236, 0.0671, 0.0674, 0.6503, 0.7631, 0.0188],
    [3700, 7.923, 0.799, 0.0911, 0.1065, 0.820, 0.940, 0.0261],
    [3800, 7.989, 1.009, 0.1216, 0.1642, 1.021, 1.144, 0.0357],
    [3900, 8.037, 1.260, 0.1599, 0.2478, 1.257, 1.380, 0.0480],
    [4000, 8.082, 1.556, 0.2073, 0.3663, 1.531, 1.648, 0.0637],
]
""" TABLE 2.06 from Hunt
            B                    C
T in K, H2, N2/CO, CO2, H2O, H2, N2/CO, CO2, H2O
"""
BCTable = [
    [1600, 16.4, 32.1, 45.7, -4.2, 20, 210, 1385, 220],
    [1700, 16.3, 32.3, 47.3, -2.5, 20, 200, 1305, 210],
    [1800, 16.2, 32.4, 48.7, -1.1, 20, 190, 1235, 195],
    [1900, 16.1, 32.6, 49.9, 0.2, 20, 180, 1170, 185],
    [2000, 16.0, 32.6, 50.9, 1.2, 15, 170, 1110, 175],
    [2100, 15.9, 32.7, 51.8, 2.2, 15, 160, 1055, 170],
    [2200, 15.8, 32.7, 52.6, 3.0, 15, 155, 1010, 160],
    [2300, 15.7, 32.8, 53.2, 3.7, 15, 150, 965, 155],
    [2400, 15.6, 32.8, 53.8, 4.4, 15, 140, 925, 145],
    [2500, 15.6, 32.8, 54.4, 5.0, 15, 135, 885, 140],
    [2600, 15.5, 32.7, 54.8, 5.5, 15, 130, 855, 135],
    [2700, 15.4, 32.7, 55.3, 6.0, 15, 125, 825, 130],
    [2800, 15.3, 32.7, 55.6, 6.4, 10, 120, 795, 125],
    [2900, 15.3, 32.6, 56.0, 6.8, 10, 120, 765, 120],
    [3000, 15.2, 32.6, 56.2, 7.1, 10, 115, 740, 120],
    [3100, 15.1, 32.6, 56.5, 7.5, 10, 110, 720, 115],
    [3200, 15.0, 32.5, 56.7, 7.7, 10, 105, 695, 110],
    [3300, 15.0, 32.4, 56.9, 8.0, 10, 105, 675, 105],
    [3400, 14.9, 32.4, 57.1, 8.3, 10, 100, 650, 105],
    [3500, 14.8, 32.3, 57.3, 8.5, 10, 95, 635, 100],
    [3600, 14.8, 32.3, 57.4, 8.7, 10, 95, 615, 100],
    [3700, 14.7, 32.2, 57.5, 8.9, 10, 90, 600, 95],
    [3800, 14.7, 32.2, 57.6, 9.1, 10, 90, 585, 95],
    [3900, 14.6, 32.1, 57.7, 9.3, 10, 85, 570, 90],
    [4000, 14.5, 32.0, 57.8, 9.4, 10, 85, 555, 90],
]

"""
Table 2.05 from Hunt
first line is found in Table 3.2 from Corner.
Add 1.987 for constant pressure values.
"Mean molecular heats over the temperature range of 300 deg K to T deg K"
"at zero desnity" (from Corner)
in calories per gram-mole per degree at constant volume
T (K), CO2, H2O, CO, H2, N2, OH, NO, O2
for monoatomic gas, take 2.980
"""
MMHTable = [
    [800, 8.896, 6.599, 5.244, 5.019, 5.179, 5.092, 5.432, 5.570],
    [1000, 9.409, 6.823, 5.403, 5.055, 5.326, 5.136, 5.592, 5.751],
    [1200, 9.824, 7.107, 5.553, 5.115, 5.468, 5.212, 5.744, 5.907],
    [1400, 10.165, 7.388, 5.684, 5.189, 5.597, 5.298, 5.869, 6.037],
    [1600, 10.449, 7.661, 5.799, 5.272, 5.712, 5.390, 5.975, 6.147],
    [1700, 10.577, 7.797, 5.852, 5.318, 5.766, 5.439, 6.024, 6.198],
    [1800, 10.690, 7.918, 5.899, 5.359, 5.814, 5.482, 6.067, 6.244],
    [1900, 10.798, 8.044, 5.945, 5.405, 5.861, 5.529, 6.110, 6.290],
    [2000, 10.896, 8.157, 5.986, 5.447, 5.904, 5.572, 6.149, 6.331],
    [2100, 10.990, 8.273, 6.012, 5.492, 5.945, 5.617, 6.186, 6.373],
    [2200, 11.075, 8.378, 6.036, 5.533, 5.983, 5.657, 6.220, 6.412],
    [2300, 11.157, 8.484, 6.086, 5.577, 6.020, 5.700, 6.252, 6.452],
    [2400, 11.233, 8.581, 6.131, 5.617, 6.053, 5.739, 6.282, 6.488],
    [2500, 11.306, 8.678, 6.162, 5.659, 6.086, 5.780, 6.310, 6.525],
    [2600, 11.373, 8.768, 6.191, 5.697, 6.116, 5.818, 6.336, 6.559],
    [2700, 11.438, 8.857, 6.218, 5.736, 6.146, 5.856, 6.361, 6.594],
    [2800, 11.498, 8.940, 6.244, 5.773, 6.173, 5.891, 6.384, 6.626],
    [2900, 11.556, 9.022, 6.269, 5.811, 6.199, 5.927, 6.406, 6.659],
    [3000, 11.611, 9.099, 6.293, 5.846, 6.224, 5.960, 6.427, 6.690],
    [3100, 11.664, 9.174, 6.315, 5.882, 6.248, 5.993, 6.448, 6.722],
    [3200, 11.714, 9.245, 6.336, 5.915, 6.270, 6.024, 6.467, 6.752],
    [3300, 11.762, 9.315, 6.357, 5.948, 6.292, 6.055, 6.486, 6.782],
    [3400, 11.808, 9.380, 6.376, 5.980, 6.312, 6.085, 6.504, 6.811],
    [3500, 11.853, 9.444, 6.395, 6.012, 6.332, 6.115, 6.521, 6.840],
    [3600, 11.895, 9.505, 6.412, 6.042, 6.350, 6.143, 6.538, 6.868],
    [3700, 11.936, 9.565, 6.429, 6.073, 6.368, 6.171, 6.554, 6.895],
    [3800, 11.975, 9.621, 6.446, 6.102, 6.385, 6.198, 6.569, 6.921],
    [3900, 12.013, 9.677, 6.462, 6.131, 6.402, 6.225, 6.585, 6.947],
    [4000, 12.050, 9.730, 6.477, 6.158, 6.418, 6.251, 6.600, 6.972],
]

"""
Table 2.08 from Hunt
E1 /100
H2, N2 & CO, CO2, H2O

for E2 * 1e^-4:
H2:     3
N2, CO: 34,
CO2:    220
H2O:    35
"""

E1Table = [
    [1600, 45, -110, -870, -935],
    [1700, 50, -95, -840, -895],
    [1800, 55, -80, -810, -855],
    [1900, 65, -70, -785, -825],
    [2000, 70, -55, -755, -795],
    [2100, 75, -40, -730, -770],
    [2200, 80, -25, -700, -745],
    [2300, 90, -15, -675, -725],
    [2400, 95, -0, -645, -705],
    [2500, 100, 15, -620, -685],
    [2600, 105, 25, -595, -665],
    [2700, 110, 40, -565, -650],
    [2800, 120, 55, -540, -635],
    [2900, 125, 65, -515, -620],
    [3000, 130, 80, -490, -605],
    [3100, 135, 95, -465, -590],
    [3200, 140, 105, -440, -580],
    [3300, 145, 120, -415, -565],
    [3400, 150, 130, -390, -555],
    [3500, 160, 145, -365, -540],
    [3600, 165, 155, -340, -530],
    [3700, 170, 170, -315, -520],
    [3800, 175, 185, -290, -510],
    [3900, 180, 195, -265, -500],
    [4000, 185, 210, -240, -490],
]


"""
(X): nbr. molecule per unit mass of propellant
{Y}: nbr. atom per unit mass of propellant

unless otherwise stated, mol/gram is assumed for
propellant work derived form Hunt.

Eq. 2.04
(CO) * (H2O) / [(CO2) * (H2)] = K0

Eq. 2.06
(N2)                    = 0.5 * {N}
(CO) + (CO2)            = {C}
(H2) + (H2O)            = 0.5 * {H}
(CO) + 2(CO2) + (H2O)   = {O}

From Eq.2.06 and Eq.2.04
Major                              Minor
(N2)  = 0.5  {N}                   | - 0.5 (N) - 0.5 (NO)
(CO)  = {C} - (CO2)                |
(H2O) = {O} - {C} - (CO2)          | - (OH) - (NO) - 2 [ (O2) + .5 (O) ]
(H2)  = 0.5 {H} + {C} - {O} + (CO2)| - 0.5 (H) + 2 (O2) + (o) + (NO) + 0.5 (OH)

n = {C} + 0.5 * {H} + 0.5 * {N}    | + (OH) + (H) + (NO) + (O2) + (O) + (N)
   CO/CO2     H2O/H2        N2

dissociaiton considered ,in descending order of significance:

2 H2O <-> 2 OH + H2             (OH)  K1 (H)
H2 <-> 2 H                      (H)   K5 (H)
2 H2O + N2 <-> 2 H2 + 2 NO      (NO)  K2 (H)
2 H2O <-> 2 H2 + O2             (O2)  K3 (H)
O2 <-> 2 O                      (O)   K4 (H)
N2 <-> 2 N                      (N)   K6 (H)
"""


def readHuntTable(T):
    for i in range(len(negDeltaTable) - 1):
        Tlow, negDeltaBlow, neghalfDeltaClow = negDeltaTable[i]
        Thigh, negDeltaBhigh, neghalfDeltaChigh = negDeltaTable[i + 1]
        if Tlow <= T <= Thigh:
            k = (T - Tlow) / (Thigh - Tlow)
            negDeltaB = negDeltaBlow * (1 - k) + negDeltaBhigh * k
            neghalfDeltaC = neghalfDeltaClow * (1 - k) + neghalfDeltaChigh * k
            break

    for i in range(len(equilibriumKT) - 1):
        Tlow, *Klow = equilibriumKT[i]
        Thigh, *Khigh = equilibriumKT[i + 1]
        if Tlow <= T <= Thigh:
            k = (1 / T - 1 / Tlow) / (1 / Thigh - 1 / Tlow)
            K = [
                10 ** (log10(Ki) * (1 - k) + log10(Kj) * k)
                for Ki, Kj in zip(Klow, Khigh)
            ]
            break

    for i in range(len(BCTable) - 1):
        Tlow, *BClow = BCTable[i]
        Thigh, *BChigh = BCTable[i + 1]
        if Tlow <= T <= Thigh:
            k = (T - Tlow) / (Thigh - Tlow)
            BC = [vi * (1 - k) + vj * k for vi, vj in zip(BClow, BChigh)]
            break

    """
    Find the internal energy of the gaseous products.
    E = MMH * (T-300K) + E1 * n/V + E2 * (n/V)**2

    the E1 and E2 corrections are only done for major products.
    """
    for i in range(len(MMHTable) - 1):
        Tlow, *MMHlow = MMHTable[i]
        Thigh, *MMHhigh = MMHTable[i + 1]
        if Tlow <= T <= Thigh:
            if Tlow >= 2000:
                k = (1 / T - 1 / Tlow) / (1 / Thigh - 1 / Tlow)
            else:
                k = (T - Tlow) / (Thigh - Tlow)
            MMH = [vi * (1 - k) + vj * k for vi, vj in zip(MMHlow, MMHhigh)]
            break

    for i in range(len(E1Table) - 1):
        Tlow, *E1low = E1Table[i]
        Thigh, *E1high = E1Table[i + 1]
        if Tlow <= T <= Thigh:
            k = (T - Tlow) / (Thigh - Tlow)
            E1 = [(vi * (1 - k) + vj * k) * 1e2 for vi, vj in zip(E1low, E1high)]
            break

    return negDeltaB, neghalfDeltaC, K, BC, MMH, E1


def balance(T, Ci, Hi, Oi, Ni, V=1 / 0.1, tol=1e-5):
    """
    Ci, Hi, Oi, Ni are in mol/g.
    Consequently,
    CO2j, OHj, Hj, NOj, O2j, Oj, Nj are also in mol/g

    T: temperature in Kelvin
    V: volume per unit mass of gas, in cc/g
    1/V is also commonly known as Delta, or load density, if the propellant
    solid covolume can be ignored.
    in literature, usually 1/V = 0.2 g/cc = 200 kg/m^3 is adopted for IB
    purposes. However, very often covolume values are given for a specific
    load condition.

        > a value of 0.089(0.084) g/cc is adopted for propellant referenced
        from ADA043778. When this value is adopted, prediciton is good to
        within 1% of tabulated.

        >It is possible that for modern software, even if the user sets a 1/V
        of 0.2 g/cc, the program will apply its own averaging routine resulting
        in different result than is calculated here, given the apparent accuracy
        of all other parameter estimation performed here, and the closeness the
        composition can be solved to published data.

    Covolume is is suitably invariant over typical internal ballistic
    tempearture, range, however it varies not-insignificantly by load density
    in test chambers, and thus must be chosen with reference to gas density
    achieved in actual gun system.

    n: number of gram-molecules per unit mass of gas, in mol/g
    """
    if T > 4000 or T < 1000:
        raise ValueError("T Not supported")

    negDeltaB, neghalfDeltaC, K, BC, MMH, E1 = readHuntTable(T)

    BH2, BN2CO, BCO2, BH2O, CH2, CN2CO, CCO2, CH2O = BC

    R = 82.06  # in cc atm /(K mol)
    sqrtVdivRT = (V / (R * T)) ** 0.5

    OHj, Hj, NOj, O2j, Oj, Nj = 0, 0, 0, 0, 0, 0

    G = Ci
    H = Oi - Ci
    I = 0.5 * Hi + Ci - Oi

    oldCO2j = -1
    CO2j = 0.5 * (
        min(G, H) - I
    )  # initiate with a value that would not result in a negative fraction

    while True:
        n = Ci + 0.5 * Hi + 0.5 * Ni + OHj + Hj + NOj + O2j + Oj + Nj  # mol/g

        N2j = 0.5 * Ni - 0.5 * Nj - 0.5 * NOj
        G = Ci
        COj = G - CO2j
        H = Oi - Ci - OHj - NOj - 2 * O2j - Oj
        H2Oj = H - CO2j
        I = 0.5 * Hi + Ci - Oi - 0.5 * Hj + 2 * O2j + Oj + NOj + 0.5 * OHj
        H2j = I + CO2j
        K0 = K[0] * exp(n / V * negDeltaB + (n / V) ** 2 * neghalfDeltaC)
        oldCO2j = CO2j
        CO2j = quadratic((1 - K0), -(G + H + K0 * I), G * H)[1]
        CO2j = max(CO2j, 0)

        if abs(oldCO2j - CO2j) / CO2j < tol:
            break

        OHj = H2Oj / (H2j) ** 0.5 * sqrtVdivRT * K[1] * exp(-20 * n / V)
        Hj = H2j**0.5 * sqrtVdivRT * K[5]
        NOj = H2Oj * N2j**0.5 / H2j * sqrtVdivRT * K[2] * exp(-20 * n / V)
        O2j = (H2Oj / H2j) ** 2 * sqrtVdivRT**2 * K[3]
        Oj = O2j**0.5 * sqrtVdivRT * K[4]
        Nj = N2j**0.5 * sqrtVdivRT * K[6]

    n = N2j + COj + CO2j + H2j + H2Oj + OHj + Hj + NOj + O2j + Oj + Nj  # mol/g

    speciesList = [
        ("CO", COj * 28.01, COj),
        ("CO2", CO2j * 44.01, CO2j),
        ("H2O", H2Oj * 18.02, H2Oj),
        ("NO", NOj * 30.01, NOj),
        ("OH", OHj * 17.01, OHj),
        ("N2", N2j * (14.01 * 2), N2j),
        ("N", Nj * 14.01, Nj),
        ("O2", O2j * (16.0 * 2), O2j),
        ("O", Oj * 16.0, Oj),
        ("H2", H2j * 1.008 * 2, H2j),
        ("H", Hj * 1.008, Hj),
    ]

    speciesList.sort(key=lambda x: x[1], reverse=True)

    B = BCO2 * CO2j + BN2CO * (N2j + COj) + BH2O * H2Oj + BH2 * H2j
    C = CCO2 * CO2j + CN2CO * (N2j + COj) + CH2O * H2Oj + CH2 * H2j

    b = (B * V**2 + n * C * V) / (V**2 + B * V + n * C)
    p = n * R * T / (V - b) / 9.869

    DeltaT = T - 300
    HCO2, HH2O, HCO, HH2, HN2, HOH, HNO, HO2 = (v * DeltaT for v in MMH)
    HH, HO, HN = (2.980 * DeltaT for _ in range(3))  # monoatomic

    E1H2, E1N2, E1CO2, E1H2O = E1
    E1CO = E1N2

    """
    2nd part of Hunt table 2.08
    """
    HH2 += E1H2 * (n / V) + 3e4 * (n / V) ** 2
    HN2 += E1N2 * (n / V) + 34e4 * (n / V) ** 2
    HCO2 += E1CO2 * (n / V) + 220e4 * (n / V) ** 2
    HH2O += E1H2O * (n / V) + 35e4 * (n / V) ** 2
    HCO += E1CO * (n / V) + 34e4 * (n / V) ** 2

    E = (  # the convention in propellant work is kinda weird
        HCO2 * CO2j
        + HH2O * H2Oj
        + HCO * COj
        + HH2 * H2j
        + HN2 * N2j
        + HOH * OHj
        + HNO * NOj
        + HO2 * O2j
        + HH * Hj
        + HO * Oj
        + HN * Nj
    )
    """add the heat of formation of the products
    according to their proportions, Hunt table 2.02, constant volume"""
    Hf = (
        COj * 26.70e3
        + CO2j * 94.02e3
        + H2Oj * 57.51e3  # this is in the gaseous form!
        + NOj * -21.50e3
        + OHj * -5.95e3
        + Nj * -84.15e3
        + Oj * -58.85e3
        + Hj * -51.53e3
    )
    f = n * T * 8.314  # 8.314 j/ mol K force constant is calculated

    """ b: covolume
        p: pressure implied by load density
        f: propellant force
        E: internal energy of products.
    """
    return E - Hf, E, n, speciesList, b, p, f


if __name__ == "__main__":
    print(balance(3100, Ci=0.02232, Hi=0.03010, Ni=0.01046, Oi=0.03469, V=1 / 0.2))
