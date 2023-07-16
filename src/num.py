import math


def sign(x):
    if x > 0:
        return 1
    elif x == 0:
        return 0
    elif x < 0:
        return -1


invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2


def gss(f, a, b, tol=1e-9, findMin=True):
    """Golden-section search. improved from the example
    given on wikipedia. Reuse half the evaluatins.

    Given a function f with a single local extremum in
    the interval [a,b], gss returns a subset interval
    [c,d] that contains the extremum with d-c <= relTol.

    a----c--d----b
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
        if (yc < yd and findMin) or (yc > yd and not findMin):
            # a---c---d
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            yc = f(c)
        else:
            # c--d---b
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            yd = f(d)

    if (yc < yd and findMin) or (yc > yd and not findMin):
        return (a, d)
    else:
        return (c, b)


def GSS(f, a, b, yRelTol, xTol=0, findMin=True):
    """
    conduct Gold Section Search using both the relative deviance of
    the functional value, and the absolute deviance of x as dual-control.
    """

    i = 0

    (a, b) = (min(a, b), max(a, b))
    ya = f(a)
    yb = f(b)
    h = b - a
    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    if xTol != 0:
        n = int(math.ceil(math.log(xTol / h) / math.log(invphi)))
    else:
        n = math.inf

    while i < n:
        if (yc < yd and findMin) or (yc > yd and not findMin):
            # a---c---d
            b = d
            d = c
            yb = yd
            yd = yc

            h = invphi * h
            c = a + invphi2 * h
            yc = f(c)

            if c == a:
                break

            if abs(ya - yd) / min(ya, yd) < yRelTol:
                break

        else:
            # c--d---b
            a = c
            c = d
            ya = yc
            yc = yd
            h = invphi * h
            d = a + invphi * h
            yd = f(d)

            if d == a:
                break

            if abs(yc - yb) / min(yc, yb) < yRelTol:
                break

        i += 1

    if (yc < yd and findMin) or (yc > yd and not findMin):
        return (a, d)
    else:
        return (c, b)


a2 = 2 / 27
a3 = 1 / 9
a4 = 1 / 6
a5 = 5 / 12
a6 = 1 / 2
a7 = 5 / 6
a8 = 1 / 6
a9 = 2 / 3
a10 = 1 / 3
a11 = 1
a12 = 0
a13 = 1


b21 = 2 / 27

b31 = 1 / 36
b32 = 1 / 12

b41 = 1 / 24
b43 = 1 / 8

b51 = 5 / 12
b53 = -25 / 16
b54 = 25 / 16

b61 = 1 / 20
b64 = 1 / 4
b65 = 1 / 5

b71 = -25 / 108
b74 = 125 / 108
b75 = -65 / 27
b76 = 125 / 54

b81 = 31 / 300
b85 = 61 / 225
b86 = -2 / 9
b87 = 13 / 900

b91 = 2
b94 = -53 / 6
b95 = 704 / 45
b96 = -107 / 9
b97 = 67 / 90
b98 = 3

b101 = -91 / 108
b104 = 23 / 108
b105 = -976 / 135
b106 = 311 / 54
b107 = -19 / 60
b108 = 17 / 6
b109 = -1 / 12

b111 = 2383 / 4100
b114 = -341 / 164
b115 = 4496 / 1025
b116 = -301 / 82
b117 = 2133 / 4100
b118 = 45 / 82
b119 = 45 / 164
b1110 = 18 / 41

b121 = 3 / 205
b126 = -6 / 41
b127 = -3 / 205
b128 = -3 / 41
b129 = 3 / 41
b1210 = 6 / 41

b131 = -1777 / 4100
b134 = -341 / 164
b135 = 4496 / 1025
b136 = -289 / 82
b137 = 2193 / 4100
b138 = 51 / 82
b139 = 33 / 164
b1310 = 12 / 41
b1312 = 1

c1 = 41 / 840
c6 = 34 / 105
c7 = 9 / 35
c8 = 9 / 35
c9 = 9 / 280
c10 = 9 / 280
c11 = 41 / 840

c_hat = -41 / 840


def RKF78(
    dFunc,
    iniVal,
    x_0,
    x_1,
    relTol,
    absTol,
    minTol=1e-16,
    abortFunc=None,
    record=None,
    show=False,
):
    """
    use Runge Kutta Fehlberg of 7(8)th power to solve the System of Equation

    Arguments:
        dFunc   : d/dx|x=x(y1, y2, y3....) = dFunc(x, y1, y2, y3...)
        iniVal  : initial values for (y1, y2, y3...)
        x_0, x_1: integration
        relTol  : relative tolerance, per component
        absTol  : absolute tolerance, per component

        abortFunc
                : optional, accepts following arguments:
            x   : current value of integrand
            ys  : current value of the SoE
            dys : estimated first derivative
                and terminates the integrator on a boolean value of True



        minTol  : optional, minimum magnitude of error

        record  : optional, if supplied will record all committed datapoints



    Returns:
        (y1, y2, y3...)|x = x_1, (e1, e2, e3....)
        where e1, e2, e3...
        is the estimated maximum deviation (in absolute) for that individual
        component
    """
    y_this = iniVal
    x = x_0

    beta = 0.84  # "safety" factor
    """
    When beta<1: 
        to reject the initial choice of h at the ith step and repeat the calculations using beta h
    , and
    When q≥1
        to accept the computed value at the ith step using the step size h, but change the step size 
        to beta h for the (i + 1)st step.
    """
    h = x_1 - x_0  # initial step size

    Rm = [0 for _ in iniVal]

    allK = [[] for _ in range(13)]

    if h == 0:
        return x, y_this, Rm

    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        if (x + h) == x:
            break  # catch the error using the final lines
        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x

        try:
            # fmt: off
            K1 = [k * h for k in dFunc(x, *y_this)]
            allK[0] = K1

            # p
            K2 = [
                k * h
                for k in dFunc(
                    x + a2 * h,
                    *[y + b21 * k1 for y, k1 in zip(y_this, *allK[:1])]
                )
            ]
            allK[1] = K2

            K3 = [
                k * h
                for k in dFunc(
                    x + a3 * h,
                    *[
                        y + b31 * k1 + b32 * k2
                        for y, k1, k2 in zip(y_this, *allK[:2])
                    ]
                )
            ]
            allK[2] = K3

            K4 = [
                k * h
                for k in dFunc(
                    x + a4 * h,
                    *[
                        y + b41 * k1 + b43 * k3
                        for y, k1, k2, k3 in zip(y_this, *allK[:3])
                    ]
                )
            ]
            allK[3] = K4

            K5 = [
                k * h
                for k in dFunc(
                    x + a5 * h,
                    *[
                        y + b51 * k1 + b53 * k3 + b54 * k4
                        for y, k1, k2, k3, k4 in zip(y_this, *allK[:4])
                    ]
                )
            ]
            allK[4] = K5

            K6 = [
                k * h
                for k in dFunc(
                    x + a6 * h,
                    *[
                        y + b61 * k1 + b64 * k4 + b65 * k5
                        for y, k1, k2, k3, k4, k5 in zip(y_this, *allK[:5])
                    ]
                )
            ]
            allK[5] = K6

            K7 = [
                k * h
                for k in dFunc(
                    x + a7 * h,
                    *[
                        y + b71 * k1 + b74 * k4 + b75 * k5 + b76 * k6
                        for y, k1, k2, k3, k4, k5, k6 in zip(y_this, *allK[:6])
                    ]
                )
            ]
            allK[6] = K7

            K8 = [
                k * h
                for k in dFunc(
                    x + a8 * h,
                    *[
                        y + b81 * k1 + b85 * k5 + b86 * k6 + b87 * k7
                        for y, k1, k2, k3, k4, k5, k6, k7 in zip(
                            y_this, *allK[:7]
                        )
                    ]
                )
            ]
            allK[7] = K8

            K9 = [
                k * h
                for k in dFunc(
                    x + a9 * h,
                    *[
                        y + b91 * k1 + b94 * k4 + b95 * k5 + b96 * k6
                        + b97 * k7 + b98 * k8
                        for y, k1, k2, k3, k4, k5, k6, k7, k8 in 
                        zip(y_this, *allK[:8])
                    ]
                )
            ]
            allK[8] = K9

            K10 = [
                k * h
                for k in dFunc(
                    x + a10 * h,
                    *[
                        y + b101 * k1 + b104 * k4 + b105 * k5 + b106 * k6
                        + b107 * k7 + b108 * k8 + b109 * k9
                        for y, k1, k2, k3, k4, k5, k6, k7, k8, k9
                        in zip(y_this, *allK[:9])
                    ]
                )
            ]
            allK[9] = K10

            K11 = [
                k * h
                for k in dFunc(
                    x + a11 * h,
                    *[
                        y + b111 * k1 + b114 * k4 + b115 * k5 + b116 * k6
                        + b117 * k7 + b118 * k8 + b119 * k9 + b1110 * k10
                        for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 
                        in zip(y_this, *allK[:10])
                    ]
                )
            ]
            allK[10] = K11

            K12 = [
                k * h
                for k in dFunc(
                    x + a12 * h,
                    *[
                        y + b121 * k1 + b126 * k6 + b127 * k7 + b128 * k8
                        + b129 * k9 + b1210 * k10
                        for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11 
                        in zip(y_this, *allK[:11])
                    ]
                )
            ]
            allK[11] = K12

            K13 = [
                k * h
                for k in dFunc(
                    x + a13 * h,
                    *[
                        y + b131 * k1 + b134 * k4 + b135 * k5 + b136 * k6
                        + b137 * k7 + b138 * k8 + b139 * k9 + b1310 * k10
                        + b1312 * k12
                        for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12
                        in zip(y_this, *allK[:12])
                    ]
                )
            ]
            allK[12] = K13

            # fmt: on

            if any(isinstance(i, complex) for k in allK for i in k):
                raise TypeError

        except (
            TypeError,
            ZeroDivisionError,
            OverflowError,
        ):  # complex value has been encountered during calculation
            # or that through unfortuante chance we got a divide by zero
            # or that a step is too large that some operation overflowed

            h *= beta
            continue

        y_next = [
            y
            + c1 * k1
            + c6 * k6
            + c7 * k7
            + c8 * k8
            + c9 * k9
            + c10 * k10
            + c11 * k11
            for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, _, _ in zip(
                y_this, *allK
            )
        ]

        te = [
            c_hat * (k1 + k11 - k12 - k13)
            for y, k1, _, _, _, _, _, _, _, _, _, k11, k12, k13 in zip(
                y_this, *allK
            )
        ]  # local truncation error, or difference per step

        """
        Extrapolating global error from local truncation error.
        Using the entire range is considered more conservative (results in larger error)
        than scaling using the remaining.
        """
        Rs = [abs(e) * (x_1 - x_0) / h for e in te]
        """
        Construct a relative error specification, comparing the global extrapolated
        error to the smaller of current and next values.
        """
        # print(*Rs)

        if absTol is None:
            R = max(
                abs(r) / (minTol + relTol * min(abs(y1), abs(y2)))
                for r, y1, y2 in zip(Rs, y_this, y_next)
            )
        else:
            R = max(
                abs(r)
                / (
                    minTol
                    + max(
                        (relTol * min(abs(y1), abs(y2))),
                        absTol,
                    )
                )
                for r, y1, y2 in zip(Rs, y_this, y_next)
            )

        delta = 1

        if R >= 1:  # error is greater than acceptable
            delta = beta * abs(1 / R) ** (1 / 8)

        else:  # error is acceptable
            y_last = y_this
            y_this = y_next
            ox = x
            x += h
            Rm = [max(Rmi, Rsi) for Rmi, Rsi in zip(Rm, Rs)]

            # print(x, *y_this)

            if show:
                print(x, y_this)
            if abortFunc is not None and abortFunc(
                x=x, ys=y_this, o_x=ox, o_ys=y_last
            ):  # premature terminating cond. is met
                return x, y_this, Rm

            if record is not None:
                record.append([x, [*y_this]])

            if R != 0:  # sometimes the error can be estimated to be 0
                delta = beta * abs(1 / R) ** (1 / 7)

            else:
                """
                if this continues to be true, we are integrating a polynomial,
                in which case the error should be independent of the step size
                Therefore we aggressively increase the step size to seek forward.
                """
                delta = 2

        h *= min(max(delta, 0.125), 2)
        """
        The step size cannot be allowed to jump too much as that would in theory, invalidate
        the assumption made to allow us to extrapolate a global error given local.
        """

    if abs((x - x_1) / (x_1 - x_0)) > relTol:
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size,"
            + " x at {}, h at {}.".format(x, h)
        )

    return x, y_this, Rm


def cubic(a, b, c, d):
    """
    returns the 3 roots of
    ax^3 + bx^2 + cx + d = 0
    assuming **real** coefficients.
    """
    if any(isinstance(i, complex) for i in (a, b, c, d)):
        raise ValueError("coefficients must be real")
    if a == 0:
        return quadratic(b, c, d)
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


def quadratic(a, b, c):
    """
    solve the quadratic equation
    defined by:
    y = a*x**2 + b * x + c
    """

    Delta = b**2 - 4 * a * c

    x_1 = 0.5 * (-b - Delta**0.5) / a
    x_2 = 0.5 * (-b + Delta**0.5) / a

    if Delta > 0:
        return min(x_1, x_2), max(x_1, x_2)
    else:
        return x_1, x_2


def secant(f, x_0, x_1, x_min=None, x_max=None, tol=1e-6, it=1000):
    """secant method that solves f(x) = 0 subjected to x in [x_min,x_max]"""

    fx_0 = f(x_0)
    fx_1 = f(x_1)

    for i in range(it):
        x_2 = x_1 - fx_1 * (x_1 - x_0) / (fx_1 - fx_0)
        if x_min is not None and x_2 < x_min:
            x_2 = 0.9 * x_min + 0.1 * x_1
        if x_max is not None and x_2 > x_max:
            x_2 = 0.9 * x_max + 0.1 * x_1

        fx_2 = f(x_2)

        if any((x_2 == x_1, x_2 == x_0)):
            raise ValueError(
                "Secant search stalled at ({},{})".format(x_2, fx_2)
            )

        if abs(fx_2) < tol:
            return x_2, fx_2
        else:
            x_0, x_1, fx_0, fx_1 = x_1, x_2, fx_1, fx_2

    raise ValueError(
        "Maximum iteration exceeded at ({},{})->({},{})".format(
            x_1, fx_1, x_2, fx_2
        )
    )


def bisect(f, x_0, x_1, tol=1e-4):
    """bisection method to numerically solve for zero
    two initial guesses must be of opposite sign.
    The root found is guaranteed to be within the range specified.

    """
    a, b = min(x_0, x_1), max(x_0, x_1)
    fa = f(a)
    fb = f(b)

    n = math.ceil(math.log((b - a) / tol, 2))

    if sign(fa) * sign(fb) > 0:
        raise ValueError("Initial Guesses Must Be Of Opposite Sign")

    for i in range(n):
        c = 0.5 * (a + b)
        fc = f(c)
        # print("a", a, "b", b)
        # print("fa", fa, "fb", fb)

        if sign(f(c)) == sign(f(a)):
            a = c
            fa = fc
        else:
            b = c
            fb = fc

    return a, b


def matMul(A, B):
    dimA = len(A), len(A[0])
    if any(len(row) != dimA[1] for row in A):
        raise ValueError("Matrix A is not consistent")
    dimB = len(B), len(B[0])
    if any(len(row) != dimB[1] for row in B):
        raise ValueError("Matrix B is not consistent")
    if dimA[1] != dimB[0]:
        raise ValueError("Dimension mistmatch for matrix A and B")

    R = [[0 for _ in range(dimB[1])] for _ in range(dimA[0])]
    BT = [*zip(*B)]

    i = 0
    for rowA in A:
        j = 0
        for columnB in BT:
            R[i][j] = sum(a * b for a, b in zip(rowA, columnB))
            j += 1
        i += 1
    return R


def solveMat(A, B):
    """
    Solve the linear system defined by Ax = B,
    where A is given in nested lists with the inner list representing the
    row entries, and B given in a flattened list representing the only column
    in the result vectory. A flattened list respresenting the x vector is
    returned.

    Specifically, we use Gauss-Jordanian elimination to calculate A^-1,
    and left multiply it such that A^-1*A*x = A^-1*B.

    """
    dim = len(A)

    if dim != len(B):
        raise ValueError("Dimension mismatch between A,x and B")

    if any(len(row) != dim for row in A):
        raise ValueError("Matrix A is not square")

    I = [[1 if i == j else 0 for i in range(dim)] for j in range(dim)]

    def swapRow(i, j):
        rowI = A[i], I[i]
        rowJ = A[j], I[j]

        A[i], I[i] = rowJ
        A[j], I[j] = rowI

    """
    print(
        *[
            " ".join("{:^8.4g}".format(v) for v in a)
            + "|"
            + " ".join("{:^8.4g}".format(v) for v in i)
            for a, i in zip(A, I)
        ],
        sep="\n"
    )
    print("")
    """
    h = 0  # pivot row
    k = 0  # pivot column

    while h < dim and k < dim:
        # choose the largest possible absolute value as partial pivot
        imax = max(
            ((A[i][k], i) for i in range(h, dim)),
            key=lambda x: x[0],
        )[1]

        if A[imax][k] == 0:
            # no pivot in this column
            k += 1
        else:
            swapRow(h, imax)
            for i in range(h + 1, dim):
                f = A[i][k] / A[h][k]
                A[i][k] = 0  # fill the lower part of pivot column
                # do for all remaining elements in current row
                for j in range(k + 1, dim):
                    A[i][j] -= A[h][j] * f

                for j in range(0, dim):
                    # apply the same operation to the identity matrix.
                    I[i][j] -= I[h][j] * f

            h += 1
            k += 1

    for i in range(dim - 1, -1, -1):
        if A[i][i] != 0:
            for j in range(i):
                f = A[j][i] / A[i][i]
                A[j][i] = 0
                for k in range(0, dim):
                    I[j][k] -= I[i][k] * f

    # convert the leading entries to 1
    for i in range(dim):
        if A[i][i] != 0:
            f = 1 / A[i][i]
            for j in range(i, dim):
                A[i][j] *= f
            for j in range(0, dim):
                I[i][j] *= f

    """

    print(
        *[
            " ".join("{:^8.4g}".format(v) for v in a)
            + "|"
            + " ".join("{:^8.4g}".format(v) for v in i)
            for a, i in zip(A, I)
        ],
        sep="\n"
    )
    """

    # now the matrix I is converted into A^-1
    Ix = matMul(I, [[b] for b in B])
    result = [i[0] for i in Ix]

    return result


def intg(f, l, u, tol=1e-3):
    """
    Integration, a.la the HP-34C. For more info see:
    "Handheld Calculator Evaluates Integrals", William M.Kahan
    Hewlett Packard Journal, August 1980 Volume 31, number 8.

    f: function, single variable.
    l: lower limit
    u: upper limit of integration
    tol: tolerance, see below

    To apply the quadrature procedure, first the problem is transformed on
    interval to:

    u              1                        given:
    ∫ f(x) dx -> a ∫ f(ax+b) dx             a = (u - l) / 2
    l             -1                        b = (u + l) / 2

    another transformation on the variable of integration eliminates the need
    to sample at either end points, which makes it possible to evaluate improper
    integrals if the asymptotes is at either end point.

    1                                        1
    ∫ f(u) du -- let u = 1.5v-0.5v**3 -> 1.5 ∫ f(1.5v-0.5v^3)*(1-v^2) dv
    -1                                      -1

    as the weight (1-v^2) is exactly 0 on both end points. We then sample
    evenly along v, take quadrature using the mid-point rule and doubling
    the number of nodes taken for each pass. This helps with suppressing
    harmonics if the function integrated is periodic. In addition, all of
    the previously calcualted quadratures can be reused in the next round,
    after dividing by half. This is especially important when function calls
     are expensive. Specifically, for pass k (k>=1 & integer) we consider 2^k-1
    points (besides the end points):

    v(i) = -1 + 2^(1-k) * i as i increments from 1 to 2^k-1 (inclusive).

                                   2^k-1
    let I(k) =  2^(1-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2)
                                   i=1

                                     2^k+1
    then I(k+1) = 2^(-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2) for every odd i + I(k)/2
                                     i=1

    as a rough approximation, the error is simply taken to be the change in
    estimated value between two successive evaluations:

    ΔI(k) = I(k) - I(k-1)

    if the quadrature procedure results in a converging result, then the error
    should decrease faster than the increment in the result, speaking in
    absolute terms. Although this is no-way guaranteed, it is convenient to
    take the increment as an upper bound on error. Therefore we check for three
    consecutive increments smaller than the specified tolerance before
    submitting the result as a good enough estimate for the integral.

    Since specifying a relative error does not work well for values extremely
    close to 0. Instead we define a error as abs(x_true - x_ref) /
    (tolerance + abs(x_ref)) x 100%
    """

    a = (u - l) / 2
    b = (u + l) / 2

    tol = abs(tol)  # ensure positive

    k = 1  # iteration counter
    I = 0  # integral counter
    c = 0  # trend counter, No. of iterations with reducing delta.

    while c < 3:
        dI = 0  # change to integral
        for i in range(1, 2**k, 2):
            v = -1 + 2 ** (1 - k) * i
            u = 1.5 * v - 0.5 * v**3
            dI += f(a * u + b) * (1 - v**2)

        dI *= 1.5 * a * 2 ** (1 - k)
        I1 = I * 0.5 + dI
        d = abs(I1 - I)  # delta, change per iteration
        I = I1
        k += 1

        if d < tol * (abs(I) + tol):
            c += 1
        else:
            c = 0

    return I, d


if __name__ == "__main__":
    print(cubic(1, 1, 2, 3))

    def df1(x, y):
        return (7 * y**2 * x**3,)

    _, v, e = RKF78(df1, (3,), 2, 0, relTol=1e-3, absTol=1e-3, minTol=1e-14)

    print(v)
    print(e)

    print(e[0] / v[0])
    print("expected value")
    print(-1 / (7 / 4 * 0**4 - 85 / 3))

    A = [[2, 1, -1], [-3, -1, 2], [-2, 1, 2]]
    print(solveMat(A, [8, -11, -3]))
