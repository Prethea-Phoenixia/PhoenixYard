from math import cos, sin, pi, atan, tan, acos

"""
#Mars
R = 3390e3
M = 6.39e23
"""

R = 1740e3
M = 7.347e22
G = 6.6743e-11


def f_v_bo(r_b, theta=None, h_bo=0):
    """
    given:
        ballistic range r_b (not including range covered under thrust),
        angle at burnt-out theta, w.r.t. vertical
        height at burnt-out h_bo,

    calculate the necessary velocity
    at burnt out v_bo
    """
    phi = r_b / R
    R_0 = R + h_bo
    GM = G * M

    if theta is None:
        theta = (phi + pi) / 4
        # this is the minimum energy solution

    v_0 = (
        (GM * R * (1 - cos(phi)))
        / (R_0**2 * sin(theta) ** 2 - R * R_0 * sin(theta - phi) * sin(theta))
    ) ** 0.5

    return theta, v_0


def f_r_b(v_bo, gamma, h_bo=0):
    """ """
    R_0 = R + h_bo
    g_0 = G * M / R_0**2
    V_0_bar_sq = v_bo**2 / (g_0 * R_0)
    K = (1 - V_0_bar_sq * (2 - V_0_bar_sq) * cos(gamma) ** 2) ** 0.5

    theta_0 = acos((1 - V_0_bar_sq * cos(gamma) ** 2) / K)
    r_b = 2 * R_0 * theta_0

    ttt = (
        (2 * v_bo)
        / (g_0 * cos(gamma))
        / (2 - V_0_bar_sq)
        * (
            K * sin(theta_0) / (1 - K * cos(theta_0))
            + 2
            / (1 - K**2) ** 0.5
            * atan((1 - K**2) ** 0.5 * tan(0.5 * theta_0) / (1 - K))
        )
    )
    return r_b, ttt


if __name__ == "__main__":
    from tabulate import tabulate

    data = []
    for i in range(1, 36):
        phi = i * 5 / 180 * pi
        r_b = R * phi
        theta, v_bo = f_v_bo(r_b, h_bo=0e3)
        r_b_prime, ttt = f_r_b(v_bo, pi / 2 - theta, 0e3)
        data.append(
            (
                phi * 180 / pi,
                theta * 180 / pi,
                r_b / 1e3,
                v_bo,
            )
        )

    print(
        tabulate(
            data,
            headers=(
                "Angle",
                "Launch Ang.",
                "Range (km)",
                "Vel (m/s)",
            ),
        ),
    )
