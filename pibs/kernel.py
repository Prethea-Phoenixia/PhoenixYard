"""
枣核弹计算相关代码
"""
from math import log, asin, acos
from num import intg, dekker


def Cx(d, ln, M_inf, rho_inf, v_inf, mu_inf):
    """
    d   : diameter of the round
    ln  : length of the head-section
    fn  : length to diameter ratio of the head-section

    Cxw : wave drag parameter
    Cxf : friction drag parameter
    """

    Re_inf = rho_inf * v_inf * ln / mu_inf

    fn = ln / d
    gamma = 1.4
    Cp0 = (2 / (gamma * M_inf**2)) * (
        (0.5 * (gamma + 1) * M_inf**2) ** (gamma / (gamma - 1))
        * ((gamma + 1) / (2 * gamma * M_inf**2 - (gamma - 1))) ** (1 / (gamma - 1))
        - 1
    )

    # T_w__T_inf = 1 + 0.9 * 0.5 * (gamma - 1) * M_inf**2
    # A = (0.5 * (gamma - 1) * M_inf**2 / T_w__T_inf) ** 0.5
    # B = (1 + 2 * (gamma - 1) * M_inf**2) / T_w__T_inf - 1
    # C1 = (2 * A**2 - B) / (B**2 + 4 * A**2) ** 0.5  # alpha
    # C2 = B / (B**2 + 4 * A**2) ** 0.5  # beta

    # print(Re_inf)

    # def f_Cf_inf(Cf_inf):
    #     return (
    #         0.242 / (A * (Cf_inf * T_w__T_inf) ** 0.5) * (asin(C1) + asin(C2))
    #         - log(Re_inf * Cf_inf, 10)
    #         + 0.5 * (1 + 2 * 0.76) * log(T_w__T_inf, 10)
    #     )

    # Cf_inf, _ = dekker(f_Cf_inf, 1e-6, 10)'

    Cf_inf = 0.032 / Re_inf**0.145 * (1 + 0.18 * M_inf**2) ** -0.5

    n = 1 / 16 * (fn + 9) if fn >= 3 else 3 / 4

    print(n)

    def f_y(x):
        return 0.5 * d * (x / ln) ** n

    def f_dy(x):
        return 0.5 * d * n / ln**n * x ** (n - 1)

    print(Cp0, Cf_inf)

    def f_dCx(x):
        y = f_y(x)
        dy = f_dy(x)
        return (8 / d**2) * (
            Cp0 * y * dy**3 / (1 + dy**2) + Cf_inf * y * (1 + dy**2) ** 0.5
        )

    Cd, _ = intg(f_dCx, 0, ln, tol=1e-6)
    print(Cd)
    return Cd


if __name__ == "__main__":
    from atmos import atmosphere

    T, P, a, rho, g, mu = atmosphere(0)

    M = 5

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    fns = []
    Cxs = []
    for fn in range(2, 20):
        fns.append(fn)
        Cxs.append(
            Cx(152e-3, ln=152e-3 * fn, M_inf=M, rho_inf=rho, v_inf=a * M, mu_inf=mu)
        )

    ax.plot(fns, Cxs)
    plt.show()
