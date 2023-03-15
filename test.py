def RK4OverTuple(dTupleFunc, iniValTuple, x_min, x_max, step=25):
    curValues = iniValTuple
    x = x_min

    h = (x_max - x_min) / step
    while x < x_max:
        K1 = dTupleFunc(x, *curValues)
        K2 = dTupleFunc(x + h * 0.5, *(i + j * h * 0.5 for i, j in zip(curValues, K1)))
        K3 = dTupleFunc(x + h * 0.5, *(i + j * h * 0.5 for i, j in zip(curValues, K2)))
        K4 = dTupleFunc(x + h, *(i + j * h for i, j in zip(curValues, K3)))
        curValues = tuple(
            c + (i + 2 * j + 2 * k + l) * h / 6
            for (c, i, j, k, l) in zip(curValues, K1, K2, K3, K4)
        )
        x += h

    return curValues


def bisect(f, x_0, x_1, tol=1e-6, it=1000):
    """bisection method to numerically solve for zero
    two initial guesses must be of opposite sign.
    The root found is guaranteed to be within the range specified.

    tol: tolerance of f(x) if x is root
    it: maximum iteration
    """
    a = x_0
    b = x_1

    if abs(f(a)) < tol:
        return a, f(a)
    elif abs(f(b)) < tol:
        return b, f(b)
    elif f(a) * f(b) > 0:
        raise ValueError("Initial Guesses Must Be Of Opposite Sign")

    for _ in range(it):
        c = (a + b) / 2

        if abs(f(c)) < tol:
            return (c, f(c))

        if f(c) * f(a) > 0:
            a = c
        else:
            b = c

    raise ValueError("Maximum iteration exceeded at ({},{})".format(c, f(c)))


if __name__ == "__main__":
    rho_p = 1600
    omega = 1.16
    f = 950e3
    theta = 0.25
    alpha = 0.001
    u1 = 5.127e-8
    n = 0.83
    e1 = 0.55e-3
    d0 = 0.55e-3
    chi = 0.75
    labda = 0.12
    mu = 0
    chi_s = 1.696
    labda_s = -0.4104
    m = 2.8

    S = 0.00266
    V0 = 0.00151
    l_g = 3.624
    phi = 1.168
    p0 = 3e7

    rho = 0.2956 * (e1 + d0 / 2)

    l0 = V0 / S
    Delta = omega / V0
    v_j = (2 * f * omega / (theta * phi * m)) ** 0.5
    Z_b = (e1 + rho) / e1

    def fPsi(Z):
        if Z < 1:
            return chi * Z * (1 + labda * Z + mu * Z**2)
        elif Z < Z_b:
            return chi_s * Z / Z_b * (1 + labda_s * Z / Z_b)
        else:
            return 1

    B = (
        S**2
        * e1**2
        / (f * phi * omega * m * u1**2)
        * (f * Delta) ** (2 * (1 - n))
    )

    # set up the ODE with regard to dt
    def ode(
        t, Z, l, v
    ):  # return d/dt for the inputs in that order, first dummy variable
        psi = fPsi(Z)
        l_psi = l0 * (1 - Delta / rho_p - Delta * (alpha - 1 / rho_p) * psi)
        p = (f * omega * psi - theta / 2 * phi * m * v**2) / (S * (l + l_psi))
        if Z < Z_b:
            dZ_dt = u1 / e1 * p**n
        else:
            dZ_dt = 0

        dl_dt = v
        dv_dt = S * p / (phi * m)

        return (dZ_dt, dl_dt, dv_dt)

    # initialize variable
    psi_i = (1 / Delta - 1 / rho_p) / (f / p0 + alpha - 1 / rho_p)
    Z_i = bisect(lambda z: fPsi(z) - psi_i, 0, 1)[0]
    l_i = 0
    v_i = 0

    i = 0
    while l_i < l_g:
        print(
            "{:^10.3g} s: {:^10.3g} m: {:^10.3g} m/s @ {:^10.3g}".format(
                i / 10000, l_i, v_i, Z_i
            )
        )
        Z_i, l_i, v_i = RK4OverTuple(ode, (Z_i, l_i, v_i), i / 10000, (i + 1) / 10000)
        i += 1
