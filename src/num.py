import math
import cmath


def sign(x):
    if x > 0:
        return 1
    elif x == 0:
        return 0
    elif x < 0:
        return -1


def intg(f, l, u, tol=1e-3):
    """
    To apply the quadrature procedure, first the problem is transformed on interval to:

    u              1                        given:
    ∫ f(x) dx -> a ∫ f(ax+b) dx             a = (u - l) / 2
    l             -1                        b = (u + l) / 2

    another transformation on the variable of integration eliminates the need to sample
    at either end points, which makes it possible to evaluate improper integrals if the
    asymptotes is at either end point.

    1                                        1
    ∫ f(u) du -- let u = 1.5v-0.5v**3 -> 1.5 ∫ f(1.5v-0.5v^3)*(1-v^2) dv
    -1                                      -1

    as the weight (1-v^2) is exactly 0 on both end points. We then sample evenly along
    v, take quadrature using the mid-point rule and doubling the number of nodes taken
    for each pass. This helps with suppressing harmonics if the function integrated is
    periodic. In addition, all of the previously calcualted quadratures can be reused
    in the next round, after dividing by half. This is especially important when funct-
    ion calls are expensive. Specifically, for pass k (k>=1 & integer) we consider 2^k-1
    points (besides the end points):

    v(i) = -1 + 2^(1-k) * i as i increments from 1 to 2^k-1 (inclusive).

                                   2^k-1
    let I(k) =  2^(1-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2)
                                   i=1

                                     2^k+1
    then I(k+1) = 2^(-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2) for every odd i + I(k)/2
                                     i=1

    as a rough approximation, the error is simply taken to be the change in estimated value
    between two successive evaluations:

    ΔI(k) = I(k) - I(k-1)

    if the quadrature procedure results in a converging result, then the error should dec-
    rease faster than the increment in the result, speaking in absolute terms. Although
    this is no-way guaranteed, it is convenient to take the increment as an upper bound
    on error. Therefore we check for three consecutive increments smaller than the specified
    tolerance before submitting the result as a good enough estimate for the integral.
    """

    tol = abs(tol)

    a = (u - l) / 2
    b = (u + l) / 2

    k = 1
    I = 0
    c = 0

    while c < 3:
        dI = 0
        for i in range(1, 2**k, 2):
            v = -1 + 2 ** (1 - k) * i
            u = 1.5 * v - 0.5 * v**3
            dI += f(a * u + b) * (1 - v**2)

        dI *= 1.5 * a * 2 ** (1 - k)
        d = abs(I * 0.5 - dI)
        I = I * 0.5 + dI
        k += 1

        if d < tol:
            c += 1
        else:
            c = 0

    return I, d


invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2


def gss(f, a, b, tol=1e-9, findMin=True):
    """Golden-section search. improved from the example
    given on wikipedia. Reuse half the evaluatins.

    Given a function f with a single local extremum in
    the interval [a,b], gss returns a subset interval
    [c,d] that contains the extremum with d-c <= tol.

    Example:
    >>> f = lambda x: (x-2)**2
    >>> a = 1
    >>> b = 5
    >>> tol = 1e-5
    >>> (c,d) = gss(f, a, b, tol)
    >>> print(c, d)
    1.9999959837979107 2.0000050911830893
    """

    (a, b) = (min(a, b), max(a, b))
    h = b - a
    if h <= tol:
        return (a, b)

    # Required steps to achieve tolerance
    n = int(math.ceil(math.log(tol / h) / math.log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    for k in range(n - 1):
        if (yc < yd and findMin) or (
            yc > yd and not findMin
        ):  # yc > yd to find the maximum
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            yc = f(c)
        else:
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            yd = f(d)

    if yc < yd:
        return (a, d)
    else:
        return (c, b)


def RKF45OverTuple(
    dTupleFunc, iniValTuple, x_0, x_1, tol, termAbv=(None, None, None)
):
    """
    Runge Kutta Fehlberg method, of the fourth and fifth order
    Even though this involves a lot more computation per cycle,
    in practice since the step size is adaptive with regard to
    the tolerance specified, significant amount of extraneous
    computation can be saved.

    In addition we specify a premature terminating condition
    where if this condition is hit, the current iteration of
    integration is returned.
    """
    y_this = iniValTuple
    x = x_0
    beta = 0.9  # "safety" factor
    h = 0.1 * (x_1 - x_0)  # initial step size
    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        if (x + h) == x:
            break  # catch the error using the final lines
        try:
            K1 = dTupleFunc(x, *y_this)
            K1 = tuple(k * h for k in K1)

            K2 = dTupleFunc(
                x + 0.25 * h, *(y + 0.25 * k1 for y, k1 in zip(y_this, K1))
            )
            K2 = tuple(k * h for k in K2)

            K3 = dTupleFunc(
                x + 0.375 * h,
                *(
                    y + (3 * k1 + 9 * k2) / 32
                    for y, k1, k2 in zip(y_this, K1, K2)
                )
            )
            K3 = tuple(k * h for k in K3)

            K4 = dTupleFunc(
                x + 12 / 13 * h,
                *(
                    y + (1932 * k1 - 7200 * k2 + 7296 * k3) / 2197
                    for y, k1, k2, k3 in zip(y_this, K1, K2, K3)
                )
            )
            K4 = tuple(k * h for k in K4)

            K5 = dTupleFunc(
                x + h,
                *(
                    y
                    + 439 / 216 * k1
                    - 8 * k2
                    + 3680 / 513 * k3
                    - 845 / 4104 * k4
                    for y, k1, k2, k3, k4 in zip(y_this, K1, K2, K3, K4)
                )
            )
            K5 = tuple(k * h for k in K5)

            # 20520
            K6 = dTupleFunc(
                x + 0.5 * h,
                *(
                    y
                    + -8 / 27 * k1
                    + 2 * k2
                    - 3544 / 2565 * k3
                    + 1859 / 4104 * k4
                    - 11 / 40 * k5
                    for y, k1, k2, k3, k4, k5 in zip(y_this, K1, K2, K3, K4, K5)
                )
            )
            K6 = tuple(k * h for k in K6)

            if any(isinstance(i, complex) for i in K1 + K2 + K3 + K4 + K5 + K6):
                raise TypeError

        except (
            TypeError,
            ZeroDivisionError,
        ):  # complex value has been encountered during calculation
            # or that through unfortuante chance we got a divide by zero
            h *= beta
            continue

        y_next = tuple(
            y + 25 / 216 * k1 + 1408 / 2565 * k3 + 2197 / 4104 * k4 - 0.2 * k5
            for y, k1, k3, k4, k5 in zip(y_this, K1, K3, K4, K5)
        )  # forth order estimation
        z_next = tuple(
            y
            + 16 / 135 * k1
            + 6656 / 12825 * k3
            + 28561 / 56430 * k4
            - 9 / 50 * k5
            + 2 / 55 * k6
            for y, k1, k3, k4, k5, k6 in zip(y_this, K1, K3, K4, K5, K6)
        )  # fifth order estimation

        epsilon = sum(
            abs(z - y) for z, y in zip(z_next, y_next)
        )  # error estimation

        if epsilon >= tol:  # error is greater than acceptable
            h *= beta * (tol / epsilon) ** 0.2
        else:  # error is acceptable
            y_this = y_next
            x += h
            if any(
                cv > pv if pv is not None else False
                for cv, pv in zip(y_this, termAbv)
            ):  # premature terminating cond. is met
                return y_this
            if epsilon != 0:  # sometimes the error can be estimated to be 0
                h *= (
                    beta * (tol / epsilon) ** 0.25
                )  # apply the new best estimate

        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x

    if abs(x - x_1) > tol:
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size,"
            + " x at {}, h at {}.".format(x, h)
        )

    return y_this


def cubic(a, b, c, d):
    """
    returns the 3 roots of
    ax^3 + bx^2 + cx + d = 0
    assuming **real** coefficients.
    """
    if any(isinstance(i, complex) for i in (a, b, c, d)):
        raise ValueError("coefficients must be real")

    Delta = (
        18 * a * b * c * d
        - 4 * b**3 * d
        + b**2 * c**2
        - 4 * a * c**3
        - 27 * a**2 * d**2
    )
    """
    Δ>0: distinct real roots.
    Δ=0: repeating real roots.
    Δ<0: one real and 2 imaginary roots.
    """
    Delta_0 = b**2 - 3 * a * c
    Delta_1 = 2 * b**3 - 9 * a * b * c + 27 * a**2 * d

    C_1 = (0.5 * (Delta_1 + (Delta_1**2 - 4 * Delta_0**3) ** 0.5)) ** (
        1 / 3
    )
    C_2 = (0.5 * (Delta_1 - (Delta_1**2 - 4 * Delta_0**3) ** 0.5)) ** (
        1 / 3
    )

    xs = []
    if any(C != 0 for C in (C_1, C_2)):
        C = C_1 if C_1 != 0 else C_2
        epsilons = (
            1,
            complex(-0.5, 3**0.5 / 2),
            complex(-0.5, -(3**0.5) / 2),
        )
        for epsilon in epsilons:
            x = -1 / (3 * a) * (b + C * epsilon + Delta_0 / (C * epsilon))
            xs.append(x)
    else:
        for _ in range(3):
            xs.append(-b / (3 * a))

    if Delta >= 0:
        xs = list(z.real for z in xs)
    else:
        # one real and 2 imaginary roots.
        xs = list(
            z.real if abs(z.imag) == min(abs(z.imag) for z in xs) else z
            for z in xs
        )
    # put the first real solution at first.
    xs.sort(key=lambda z: 1 if isinstance(z, complex) else 0)
    return tuple(xs)


def bisect(f, x_0, x_1, tol, it=100):
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
    elif sign(f(a)) * sign(f(b)) > 0:
        raise ValueError("Initial guesses must be of opposite sign")

    for _ in range(it):
        c = (a + b) / 2

        if abs(f(c)) < tol:
            return (c, f(c))

        if sign(f(c)) == sign(f(a)):
            a = c
        else:
            b = c

    raise ValueError("Maximum iteration exceeded at ({},{})".format(c, f(c)))


def secant(f, x_0, x_1, x_min=None, x_max=None, tol=1e-6, it=100):
    """secant method that solves f(x) = 0 subjected to x in [x_min,x_max]"""
    if x_min is not None:
        if x_0 < x_min:
            x_0 = x_min
        if x_1 < x_min:
            x_1 = x_min
    if x_max is not None:
        if x_0 > x_max:
            x_0 = x_max
        if x_1 > x_max:
            x_1 = x_max

    fx_0 = f(x_0)
    fx_1 = f(x_1)

    if x_0 == x_1 or fx_0 == fx_1:
        raise ValueError("Initial guess must evaluate to different values")

    for _ in range(it):
        x_2 = x_1 - fx_1 * (x_1 - x_0) / (fx_1 - fx_0)
        if x_min is not None and x_2 < x_min:
            x_2 = x_min
        if x_max is not None and x_2 > x_max:
            x_2 = x_max
        x_0, x_1, fx_0, fx_1 = x_1, x_2, fx_1, f(x_2)
        if abs(fx_1) < tol or (fx_0 == fx_1):
            return x_1, fx_1

    raise ValueError("Maximum iteration exceeded at ({},{})".format(x_1, fx_1))


if __name__ == "__main__":
    from random import uniform

    for _ in range(10):
        print(
            cubic(
                uniform(-1, 1), uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)
            )
        )
