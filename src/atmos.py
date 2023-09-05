from math import cos, pi, exp

"""
USA Standard Atmospere 1976, from left to right:
Geopotential Altitude Above Sea Level (meter)
Temperature (C)
Acceleration due Gravity (g)
Absolute Pressure (10^4 Pa)
Density (kg/m^3)
Dynamic viscosity (10^5 Ns / m^2)

These data are used for checking the validity of the
ICAO atmosphere calculations
"""
USSA1976 = [
    [-1000, 21.50, 9.810, 11.39, 1.347, 1.821],
    [0, 15.00, 9.807, 10.13, 1.225, 1.789],
    [1000, 8.50, 9.804, 8.988, 1.112, 1.758],
    [2000, 2.00, 9.801, 7.950, 1.007, 1.726],
    [3000, -4.49, 9.797, 7.012, 0.9093, 1.694],
    [4000, -10.98, 9.794, 6.166, 0.8194, 1.661],
    [5000, -17.47, 9.791, 5.405, 0.7364, 1.628],
    [6000, -23.96, 9.788, 4.722, 0.6601, 1.595],
    [7000, -30.45, 9.785, 4.111, 0.5900, 1.561],
    [8000, -36.94, 9.782, 3.565, 0.5258, 1.527],
    [9000, -43.42, 9.779, 3.080, 0.4671, 1.493],
    [10000, -49.90, 9.776, 2.650, 0.4135, 1.458],
    [15000, -56.50, 9.761, 1.211, 0.1948, 1.422],
    [20000, -56.50, 9.745, 0.5529, 0.08891, 1.422],
    [25000, -51.60, 9.730, 0.2549, 0.04008, 1.448],
    [30000, -46.64, 9.715, 0.1197, 0.01841, 1.475],
    [40000, -22.80, 9.684, 0.0287, 0.003996, 1.601],
    [50000, -2.5, 9.654, 0.007978, 0.001027, 1.704],
    [60000, -26.13, 9.624, 0.002196, 0.0003097, 1.584],
    [70000, -53.57, 9.594, 0.00052, 0.00008283, 1.438],
    [80000, -74.51, 9.564, 0.00011, 0.00001846, 1.321],
]


"""
ICAO adopted vertical temperature gradients.
"""
ICAOTbbeta = [
    [-5.00, -6.50],
    [0.00, -6.50],
    [11.00, 0.00],
    [20.00, +1.00],
    [32.00, +2.80],
    [47.00, 0.00],
    [51.00, -2.80],
    [71.00, -2.00],
    [84.8520, 0],
    [89.7157, 0],
]

R_e = 6356766  # nominal earth's radius


def atmosphere(h, T0=288.15, P0=101.325e3):
    """
    ICAO standard atmosphere
        h   : Geometrical altitude
     #   lat : Latitiude (in degrees)
        p0  : pressure ASL
    """

    # Standard acceleration due to gravtiy, it conforms with latitude phi =
    # 42deg 32min 33sec using Lambert's equation of the acceleration due to
    # gravity as a function of latitude

    # 45 deg 32 min 33 sec is approximated by 45.5425 deg

    # phi = lat * pi / 180

    # gphi = 9.80616 * (1 - 2.6373e-3 * cos(2 * phi) + 5.9e-6 * cos(2 * phi) ** 2)
    # g0 = gphi

    g0 = 9.80665  # nominal earth's gravitational acceleartion
    r = R_e

    g = g0 * (r / (r + h)) ** 2  # approximate local gravitational acceleration
    # <0.001% error at 60km altitude
    H = (
        r * h / (r + h) * 1e-3
    )  # H: approximate geopotential height in kilometer

    if H < ICAOTbbeta[0][0] or H > ICAOTbbeta[-1][0]:
        raise ValueError(
            "Geopotential altitude {:.3f} km is outside of range {:} to {:} km".format(
                H, ICAOTbbeta[0][0], ICAOTbbeta[-1][0]
            )
        )

    R = 287.05287  # J/(K · kg), R / M

    if H >= 0:
        i, j = 1, 2
        Pb = P0
        T_b = T0
        H_b, beta = ICAOTbbeta[i]  # bottom
        H_c, beta_c = ICAOTbbeta[j]  # ceiling

        while H_c < H:
            T_c = T_b + beta * (H_c - H_b)
            if beta != 0:
                Pb = Pb * (1 + beta / T_b * (H_c - H_b)) ** (
                    -g0 / (beta * 1e-3 * R)
                )
            else:
                Pb = Pb * exp(-g0 * (H_c - H_b) * 1e3 / (R * T_c))

            i += 1
            j += 1

            T_b = T_c
            H_b = H_c
            beta = beta_c

            H_b, beta = ICAOTbbeta[i]  # bottom
            H_c, beta_c = ICAOTbbeta[j]  # ceiling

        T = T_b + beta * (H - H_b)
        if beta != 0:
            P = Pb * (1 + beta / T_b * (H - H_b)) ** (-g0 / (beta * 1e-3 * R))
        else:
            P = Pb * exp(-g0 * (H - H_b) * 1e3 / (R * T))

    else:
        i, j = 0, 1
        Pc = P0
        T_c = T0
        H_b, beta = ICAOTbbeta[i]  # bottom
        H_c, beta_c = ICAOTbbeta[j]  # ceiling

        T = T_c + beta * (H - H_c)
        if beta != 0:
            P = Pc * (1 + beta / T_c * (H - H_c)) ** (-g0 / (beta * 1e-3 * R))
        else:
            P = Pc * exp(-g0 * (H - H_c) * 1e3 / (R * T))

    # derived terms:
    kappa = 1.4
    a = (kappa * R * T) ** 0.5

    rho = P / (R * T)

    """
    Return values:
        T : Temperature in K
        P : Pressure in Pa
        a : Speed of sound in m/s
      rho : Density in kg m^-3
        g : Acceleration due gravity.
    """
    return T, P, a, rho, g


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    import numpy as np

    ax = plt.subplot()

    x = np.linspace(-4.9, 91, 1000)
    ax.plot(x, [atmosphere(h * 1e3)[0] - 273.15 for h in x])
    ax.scatter(
        [line[0] * 1e-3 for line in USSA1976],
        [line[1] for line in USSA1976],
        s=4,
    )

    ax.plot(x, [atmosphere(h * 1e3)[1] * 1e-4 for h in x])
    ax.scatter(
        [line[0] * 1e-3 for line in USSA1976],
        [line[3] for line in USSA1976],
        s=4,
    )
    plt.show()
