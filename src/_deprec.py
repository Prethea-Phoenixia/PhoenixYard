"""
Depreceated code that i spent way too much time on and can't be parted
with. 
"""



print("\nanalytical")
print(tabulate(test.analyze(1e-3), headers=("tag", "t", "l", "phi", "v", "p")))






 def propagate(self, l_g=None, tol=1e-5, maxiter=100):
        """
        this is a stripped down version used in numerical optimization
        where the length of gun is considerd as undefined, and instead
        the general burning characteristic is sought after.

        returns peak pressure and velocity & shot travel if not supplied
        with a gun length.
        """

        N = 1
        Delta_Z = self.Z_b - self.Z_0
        Z_i = self.Z_0

        Z_j = Z_i + Delta_Z / N
        t_bar_i, l_bar_i, v_bar_i = 0, 0, 0
        p_bar_i = self.p_0 / (self.f * self.Delta)

        """
        Same as before, although here the terminating condition is
        set as a decline in pressure
        """
        while Z_i < self.Z_b:  # terminates if burnout is achieved
            try:
                t_bar_j, l_bar_j, v_bar_j = RKF45OverTuple(
                    self._ode_Z,
                    (t_bar_i, l_bar_i, v_bar_i),
                    Z_i,
                    Z_j,
                    tol=tol,
                    imax=maxiter * N**0.25,
                )
                p_bar_j = self._fp_bar(Z_j, l_bar_j, v_bar_j)

                if p_bar_i > p_bar_j:
                    if p_bar_i - p_bar_j > tol or l_bar_i == 0:
                        N *= 2
                        Z_j = Z_i + Delta_Z / N
                    else:
                        break  # l_bar_i is solved to within a tol of l_bar_g
                else:
                    t_bar_i, l_bar_i, v_bar_i, p_bar_i = (
                        t_bar_j,
                        l_bar_j,
                        v_bar_j,
                        p_bar_j,
                    )
                    Z_i = Z_j
                    """
                    this way the group of values denoted by _i is always updated
                    as a group.
                    """
                    Z_j += Delta_Z / N
                    if Z_j > self.Z_b:
                        Z_j = self.Z_b
            except ValueError:
                N *= 2
                Z_j = Z_i + Delta_Z / N

            if N > 0.1 / tol:
                raise ValueError("Excessive division", N)

        def f(t_bar):
            Z, l_bar, v_bar = RKF45OverTuple(
                self._ode_t,
                (self.Z_0, 0, 0),
                0,
                t_bar,
                tol=tol,
                imax=maxiter,
            )
            return self._fp_bar(Z, l_bar, v_bar)

        """

        # tolerance is specified a bit differently for gold section search
        t_bar_p_1, t_bar_p_2 = gss(
            f, t_bar_i, t_bar_i + tol, tol=tol, findMin=False
        )
        t_bar_p = (t_bar_p_1 + t_bar_p_2) / 2

        Z_p, l_bar_p, v_bar_p = RKF45OverTuple(
            self._ode_t, (self.Z_0, 0, 0), 0, t_bar_p, tol=tol, imax=maxiter
        )
        """

        if l_g is not None:
            t_bar_e, Z_bar_e, v_bar_e = RKF45OverTuple(
                self._ode_l,
                (t_bar_i, Z_i, v_bar_i),
                l_bar_i,
                l_g,
                tol=tol,
            )
            return v_bar_e * self.v_j
        else:
            return (
                p_bar_i * (self.f * self.Delta),
                v_bar_i * self.v_j,
                l_bar_i * self.l_0,
            )




    def analyze(
        self,
    ):
        """
        Run the Mayer-Hart analytical model, developed and used by the
        ballistics research library in the period of 1945-1970s

        Following assumptions:

        > no shot starting pressure
        > covolume of propellant is exactly the same as the initial volume
            of the propellant
        > propellant burns steadily
        > burn rate is proportional to pressure

        with the S.O.E:
                    ψ = z
                   dz = (p / I_k) dt
                   dv = S p / (φ m) dt
                   dl = v dt
        S p (l + l_1) = f ω ψ - 0.5 θ φ m v^2
        Where:
        l_1 = l_0 (1 - αΔ)
        I_k = e / u

        e: grain thickness
        u: linear burn rate coefficient

        psi = z is taken from 0 to 1. This is exactly correct for infinitely
        long cylinder shaped grain, and closely approximates the behaviour of
        strip shaped grains. In practice these shapes are all (slightly) degressive
        burning, i.e. the psi-z curve curves towards higher psi between (0,1)

        Although for multi-perf grained propellant, the exponential burn
        rate phenomnon is important for arriving at the correct result,
        since they tends to exhibit strong progressive or degressive tendencies
        which drives production of peak pressure,

        in practice simpler shapes that are close to constant burning is well
        modeled by assuming linear burn rate.

        In the following context, phase I refers to shot start to burnout
        subscript k indicates condition at end of phase I.
        phase II refers to the adiabatic expansion phase.
        """

        I_k = self.e_1 / self.u_0
        # equations for phase I
        l_1 = self.l_0 * (1 - self.alpha * self.Delta)
        B = (self.S * I_k) ** 2 / (self.f * self.omega * self.phi * self.m)
        v_k = self.S * I_k / (self.phi * self.m)  # velocity at end of burn
        p_1 = self.f * self.omega / (self.S * l_1)

        def _v1_psi(psi):
            return v_k * psi

        def _l1_psi(psi):
            return l_1 * (
                (1 - 0.5 * B * self.theta * psi) ** (-2 / self.theta) - 1
            )

        def _p1_psi(psi):
            return (
                p_1
                * (psi - 0.5 * B * self.theta * psi**2)
                * (1 - 0.5 * B * self.theta * psi) ** (2 / self.theta)
            )

        def _v1_l(l):
            y = l / l_1
            return (
                2
                * self.f
                * self.omega
                / (self.theta * self.S * I_k)
                * (1 - (1 + y) ** (-0.5 * self.theta))
            )

        def _p1_l(l):
            y = l / l_1
            return (
                p_1
                * 2
                / (B * self.theta)
                * (1 - (1 + y) ** (-0.5 * self.theta))
                * (1 + y) ** (-0.5 * self.theta - 1)
            )

        psi_m = 1 / (1 + B * self.theta)

        l_m = _l1_psi(psi_m)
        v_m = _v1_psi(psi_m)
        p_m = _p1_psi(psi_m)

        """
            p_m = _p1_l(l_m)
            v_m = _v1_l(l_m)
            """
        # interface values between phase I and phase II

        l_k = _l1_psi(1)
        p_k = _p1_psi(1)

        # equations for phase II

        def _v2_l(l):
            return (
                self.v_j
                * (
                    1
                    - ((l_1 + l_k) / (l_1 + l)) ** self.theta
                    * (1 - 0.5 * B * self.theta)
                )
                ** 0.5
            )

        def _p2_l(l):
            return (
                self.f
                * self.omega
                / (self.S * (l + l_1))
                * (1 - (_v2_l(l) / self.v_j) ** 2)
            )

        if self.l_g > l_k:
            p_e = _p2_l(self.l_g)
            v_e = _v2_l(self.l_g)
            psi_e = 1.0
        else:
            p_e = _p1_l(self.l_g)
            v_e = _v1_l(self.l_g)
            psi_e = v_e / v_k

        data = []

        # t, l, phi, v, p
        data.append(("SHOT START", 0, 0, 0, 0, 0))
        data.append(("PEAK PRESSURE", 0, l_m, psi_m, v_m, p_m))
        data.append(("BURNOUT", 0, l_k, 1.0, v_k, p_k))
        data.append(("SHOT EXIT", 0, self.l_g, psi_e, v_e, p_e))
        data.sort(key=lambda x: x[1])
        return data



 """
            Solving for peak pressure
        """

        def _x_m(p_m):
            """
            x_m being a function of p_m, therefore must be solved iteratively
            apparently the equation "- self.chi * self.labda" is correct
            """
            return K_1 / (
                B_0 * (1 + self.theta) / (1 + p_m / (self.f * delta_1))
                - chi_prime * labda_prime
            )

        """
            iteratively solve x_m in the range of (0,x_k)
            x_k signifies end of progressive burn
        
            tolerance is specified against the unitless scaled value
            for each parameter.
            """
        p_m_i = self.p_0 * 2  # initial guess, 100MPa

        for i in range(it):  # 8192
            x_m_i = _x_m(p_m_i)
            if x_m_i > x_k:
                x_m_i = x_k
            elif x_m_i < 0:
                x_m_i = 0

            l_m_i, t_m_i = propagate(x_m_i, it=it)  # see above
            p_m_j = _p(x_m_i, l_m_i)

            if abs(p_m_i - p_m_j) > self.f * self.Delta * tol:
                p_m_i = p_m_j
            else:
                break

        # value reflecting peak pressure point
        p_m = p_m_j
        x_m = x_m_i
        l_m = l_m_i
        t_m = t_m_i
        v_m = _v(x_m)
        psi_m = _psi(x_m)



def RKF45(dFunc, iniVal, x_0, x_1, tol, termAbv=None):
    """
    Runge Kutta Fehlberg method, of the fourth and fifth order
    Even though this involves a lot more computation per cycle,
    in practice since the step size is adaptive with regard to
    the tolerance specified, significant amount of extraneous
    computation can be saved.

    In addition we specify a premature terminating condition
    where if this condition is hit, the current iteration of
    integration is returned.

    For a RKF-45 method, the local truncation error is O(h^5)
    while the global accumulated error is O(h^4)
    """
    # i = 0
    if termAbv is None:
        termAbv = tuple(None for _ in iniVal)
    y_this = iniVal
    x = x_0
    beta = 0.9  # "safety" factor
    h = tol * (x_1 - x_0)  # initial step size
    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        if (x + h) == x:
            break  # catch the error using the final lines

        try:
            K1 = tuple(k * h for k in dFunc(x, *y_this))
            K2 = tuple(
                k * h
                for k in dFunc(
                    x + 0.25 * h, *(y + 0.25 * k1 for y, k1 in zip(y_this, K1))
                )
            )
            K3 = tuple(
                k * h
                for k in dFunc(
                    x + 0.375 * h,
                    *(
                        y + (3 * k1 + 9 * k2) / 32
                        for y, k1, k2 in zip(y_this, K1, K2)
                    )
                )
            )
            K4 = tuple(
                k * h
                for k in dFunc(
                    x + 12 / 13 * h,
                    *(
                        y + (1932 * k1 - 7200 * k2 + 7296 * k3) / 2197
                        for y, k1, k2, k3 in zip(y_this, K1, K2, K3)
                    )
                )
            )
            K5 = tuple(
                k * h
                for k in dFunc(
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
            )
            K6 = tuple(
                k * h
                for k in dFunc(
                    x + 0.5 * h,
                    *(
                        y
                        + -8 / 27 * k1
                        + 2 * k2
                        - 3544 / 2565 * k3
                        + 1859 / 4104 * k4
                        - 11 / 40 * k5
                        for y, k1, k2, k3, k4, k5 in zip(
                            y_this, K1, K2, K3, K4, K5
                        )
                    )
                )
            )

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

        R = (
            sum(abs(z - y) for z, y in zip(z_next, y_next)) / h
        )  # error estimation

        delta = 1
        if R >= tol:  # error is greater than acceptable
            delta = beta * abs(tol / R) ** 0.2

        else:  # error is acceptable
            y_this = y_next
            x += h
            if any(
                cv > pv if pv is not None else False
                for cv, pv in zip(y_this, termAbv)
            ):  # premature terminating cond. is met
                return y_this
            if R != 0:  # sometimes the error can be estimated to be 0
                delta = (
                    beta * abs(tol / R) ** 0.25
                )  # apply the new best estimate
            else:
                """
                if this continues to be true, we are integrating a polynomial,
                in which case the error should be independent of the step size
                Therefore we aggressively increase the step size to seek forward.
                """
                delta = 2

        h *= min(max(delta, 0.3), 2)  # to ensure that this does not jump
        # print(*y_this)
        # print(h, tol)
        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x
    if abs(x - x_1) > tol:
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size,"
            + " x at {}, h at {}.".format(x, h)
        )

    return y_this




def RKF23(dFunc, iniVal, x_0, x_1, tol, termAbv=None):
    """
    An extension of the thinking of Fehlberg down to a lower
    order Runge Kutta integrator.
    dFunc is interpreted as:
        d/dx|_x(y1,y2,y3...)  = dFUnc(x,y1,y2,y3....)
    """
    if termAbv is None:
        termAbv = tuple(None for _ in iniVal)
    y_this = iniVal
    x = x_0
    beta = 0.9  # "safety" factor
    h = tol * (x_1 - x_0)  # initial step size
    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        # print(x, y_this, h)
        if (x + h) == x:
            break  # catch the error using the final lines

        try:
            K1 = tuple(k * h for k in dFunc(x, *y_this))
            K2 = tuple(
                k * h
                for k in dFunc(x + h, *(y + k1 for y, k1 in zip(y_this, K1)))
            )
            K3 = tuple(
                k * h
                for k in dFunc(
                    x + 0.5 * h,
                    *(y + 0.25 * (k1 + k2) for y, k1, k2 in zip(y_this, K1, K2))
                )
            )

            if any(isinstance(i, complex) for i in K1 + K2 + K3):
                raise TypeError

        except (
            TypeError,
            ZeroDivisionError,
        ):  # complex value has been encountered during calculation
            # or that through unfortuante chance we got a divide by zero
            h *= beta
            continue

        y_next = tuple(
            y + 0.5 * (k1 + k2) for y, k1, k2 in zip(y_this, K1, K2)
        )  # 2nd order estimation
        z_next = tuple(
            y + (k1 + k2 + 4 * k3) / 6
            for y, k1, k2, k3 in zip(y_this, K1, K2, K3)
        )  # 3rd order estimation

        R = (
            max(abs((z - y) / y) for z, y in zip(z_next, y_next)) / h
        )  # error estimation

        delta = 1
        if R >= tol:  # error is greater than acceptable
            delta = beta * (tol / R) ** (1 / 3)
        else:  # error is acceptable
            y_this = y_next
            x += h
            if any(
                cv > pv if pv is not None else False
                for cv, pv in zip(y_this, termAbv)
            ):  # premature terminating cond. is met
                return y_this
            if R != 0:  # sometimes the error can be estimated to be 0
                delta = beta * (tol / R) ** 0.5  # apply the new best estimate
            else:
                """
                if this continues to be true, we are integrating a polynomial,
                in which case the error should be independent of the step size
                Therefore we aggressively increase the step size to seek forward.
                """
                delta = 2

        h *= delta
        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x

    if abs(x - x_1) > tol:
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size,"
            + " x at {}, h at {}.".format(x, h)
        )

    return y_this


def RK45(dFunc, iniVal, x_0, x_1, steps):
    y_this = iniVal
    x = x_0
    h = (x_1 - x_0) / steps
    for i in range(steps):
        K1 = dFunc(x, *y_this)
        K1 = tuple(k * h for k in K1)

        K2 = dFunc(x + 0.5 * h, *(y + 0.5 * k1 for y, k1 in zip(y_this, K1)))
        K2 = tuple(k * h for k in K2)

        K3 = dFunc(
            x + 0.5 * h, *(y + 0.5 * k2 for y, k1, k2 in zip(y_this, K1, K2))
        )
        K3 = tuple(k * h for k in K3)

        K4 = dFunc(
            x + h, *(y + k3 for y, k1, k2, k3 in zip(y_this, K1, K2, K3))
        )
        K4 = tuple(k * h for k in K4)

        y_next = tuple(
            y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            for y, k1, k2, k3, k4 in zip(y_this, K1, K2, K3, K4)
        )

        x += h
        y_this = y_next

    return y_this

def findExtBnd(f, a, b, tol=1e-9, findMin=True):
    """
    Although theoretically this saves a few cycles comapred to straight
    gold section search, it is also conditionally worse, and due to the
    extra calculation involved, it is actually slower than straight
    gss. Therefore it is only here as a testament to the effrot put into
    making the code faster.
    """
    n = int(math.ceil(math.log(tol / (b - a)) / math.log(invphi)))
    (p, q) = (min(a, b), max(a, b))
    if q - p <= tol:
        return (p, q)

    yp = f(p)
    yq = f(q)
    r = 0.5 * (a + b)
    yr = f(r)

    i = 0

    # p---r---q#
    while (q - r) > tol or (r - p) > tol:
        if r == p or r == q:
            gss = True
        else:
            alpha = (yq - yp) / (q - p)
            beta = (yr - yp - alpha * (r - p)) / ((r - p) * (r - q))
            if (beta > 0 and findMin) or (beta < 0 and not findMin):
                x = 0.5 * (a + b - alpha / beta)
                yx = f(x)
                if p < x < q and (
                    (yx < yr and findMin) or (yx > yr and not findMin)
                ):
                    if x < r:
                        q = r
                        yq = yr
                        r = x
                        yr = yx
                    else:
                        p = r
                        yp = yr
                        r = x
                        yr = yx
                    gss = False
                else:
                    gss = True
            else:
                gss = True
        if gss:
            a = p
            b = q
            h = b - a

            c = a + invphi2 * h
            d = a + invphi * h
            yc = f(c)
            yd = f(d)

            if (yc < yd and findMin) or (yc > yd and not findMin):
                # a---c---d     b
                # p---r---q
                q = d
                r = c
                yr = yc
            else:
                # a     c--d---b
                #       p--r---q
                p = c
                r = d
                yr = yd

        i += 1

    if (yp < yq and findMin) or (yp > yq and not findMin):
        return (p, r)
    else:
        return (r, q)



def RKF45(dFunc, iniVal, x_0, x_1, tol, absTol=1e-14, termAbv=None):
    """
    use Runge Kutta Fehlberg of 4(5)th power to solve the System of Equation

    Arguments:
        dFunc   : d/dx|x=x(y1, y2, y3....) = dFunc(x, y1, y2, y3...)
        iniVal  : initial values for (y1, y2, y3...)
        x_0, x_1: integration
        tol     : relative tolerance, per component
        absTol  : absolute tolerance, per component
        termAbv : premature termination condition

    Returns:
        (y1, y2, y3...)|x = x_1, (e1, e2, e3....)
        where e1, e2, e3...
        is the estimated maximum deviation (in absolute) for that individual
        component
    """
    if termAbv is None:
        termAbv = tuple(None for _ in iniVal)
    y_this = iniVal
    x = x_0
    beta = 0.84  # "safety" factor
    h = x_1 - x_0  # initial step size
    Rm = tuple(0 for _ in iniVal)
    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        if (x + h) == x:
            break  # catch the error using the final lines
        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x

        try:
            K1 = tuple(k * h for k in dFunc(x, *y_this))
            K2 = tuple(
                k * h
                for k in dFunc(
                    x + 0.25 * h, *(y + 0.25 * k1 for y, k1 in zip(y_this, K1))
                )
            )
            K3 = tuple(
                k * h
                for k in dFunc(
                    x + 0.375 * h,
                    *(
                        y + (3 * k1 + 9 * k2) / 32
                        for y, k1, k2 in zip(y_this, K1, K2)
                    )
                )
            )
            K4 = tuple(
                k * h
                for k in dFunc(
                    x + 12 / 13 * h,
                    *(
                        y + (1932 * k1 - 7200 * k2 + 7296 * k3) / 2197
                        for y, k1, k2, k3 in zip(y_this, K1, K2, K3)
                    )
                )
            )
            K5 = tuple(
                k * h
                for k in dFunc(
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
            )
            K6 = tuple(
                k * h
                for k in dFunc(
                    x + 0.5 * h,
                    *(
                        y
                        + -8 / 27 * k1
                        + 2 * k2
                        - 3544 / 2565 * k3
                        + 1859 / 4104 * k4
                        - 11 / 40 * k5
                        for y, k1, k2, k3, k4, k5 in zip(
                            y_this, K1, K2, K3, K4, K5
                        )
                    )
                )
            )

            if any(isinstance(i, complex) for i in K1 + K2 + K3 + K4 + K5 + K6):
                raise TypeError

        except (
            TypeError,
            ZeroDivisionError,
            OverflowError,
        ) as e:  # complex value has been encountered during calculation
            # or that through unfortuante chance we got a divide by zero
            # print(e)
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
        )

        te = tuple(z - y for z, y in zip(z_next, y_next))

        Rs = tuple(abs(e / h) for e in te)  # projected global error
        # print("Rs:", Rs)
        R = max(
            r / (absTol + tol * (abs(y) + abs(k1)))
            for r, y, k1 in zip(Rs, y_this, K1)
        )  # relative error as compared to specification
        # print("R :", R)
        delta = 1
        if R > 1:  # error is greater than acceptable
            delta = beta * abs(1 / R) ** 0.2

        else:  # error is acceptable
            y_this = y_next
            x += h
            Rm = tuple(max(Rmi, Rsi) for Rmi, Rsi in zip(Rm, Rs))
            if any(
                cv > pv if pv is not None else False
                for cv, pv in zip(y_this, termAbv)
            ):  # premature terminating cond. is met
                return y_this, Rm
            if R != 0:  # sometimes the error can be estimated to be 0
                delta = beta * abs(1 / R) ** 0.25  # apply the new best estimate
            else:
                """
                if this continues to be true, we are integrating a polynomial,
                in which case the error should be independent of the step size
                Therefore we aggressively increase the step size to seek forward.
                """
                delta = 4
        # print("Delta", delta)
        h *= min(max(delta, 0.25), 4)  # to ensure that this does not jump

    if abs((x - x_1) / (x_1 - x_0)) > tol:
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size,"
            + " x at {}, h at {}.".format(x, h)
        )

    return y_this, Rm

    from tabulate import tabulate

    print("verifying RKF45 is of 4th order")
    tols = [1e-5, 1e-6, 1e-7]
    devs = []
    for i in range(0, 6):
        dev_i = [i]
        for tol in tols:
            val, e = RKF45(
                lambda x, y: (x**i,), (0,), 4, 5, tol=tol, absTol=1e-16
            )

            act = (5 ** (i + 1) - (4) ** (i + 1)) / (i + 1)
            dev_i.append((act - val[0]) / act)

        devs.append(dev_i)

    print(tabulate(devs, headers=("order", *(str(t) for t in tols), "isSame?")))



    def analyze(self, it=250, tol=1e-5):
        """run the psi-bar analytical solution on the defined weapon

        key assumptions:
            >burn rate scales linearlly with pressure
            >propellant volume fraction (ψ, psi) is a polynomial
             function of linear burn fraction (Z),
             (this involves rewriting the psi-z polynomial to
             match the starting and end points of the cubic equation,
             which is not exactly rigorous, but is necessary
             for analytically solving the SOE.)
            >l_psi does not vary much for the entire burn duration
             (since with current propellant the increase in free
             chamber volume as a consequence of burning is out-
             weighed by the effect of covolume, this is approxima-
             tely true for moderate chamber loading fractions.)

        when u_1 refers to the linear burn rate model, it is substituted
        for u_0, deviating from the source's use.

        alternate chi, labda for quadratic form equation we use here
        to approximate the more rigorous cubic polynomial form equations.
        values are chosen such that the shot start and fracture points
        match exactly.

        Since this is an approximation, we do not specify an accuracy and
        instead supplies an iteration number. This allows it to be used
        where the calculation finishing and returning is required, such
        as say numerial optimisation routines
        """
        Z_0 = self.Z_0
        psi_s = self.chi * (1 + self.labda + self.mu)

        labda_prime = (psi_s * Z_0 - self.psi_0) / (
            self.psi_0 - Z_0**2 * psi_s
        )
        chi_prime = psi_s / (1 + labda_prime)

        # the equivalent of B used for linear burn rates.
        B_0 = (
            self.S**2
            * self.e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_0**2)
        )

        # labda_prime = self.labda
        # chi_prime = self.chi
        sigma_0 = (1 + 4 * labda_prime / chi_prime * self.psi_0) ** 0.5

        """
        # only applicable in the quadratic form case
        if labda_prime == 0:
            Z_0 = self.psi_0 / chi_prime
        else:
            Z_0 = (sigma_0 - 1) / (2 * labda_prime)
        """
        K_1 = chi_prime * sigma_0
        B_1 = B_0 * self.theta * 0.5 - chi_prime * labda_prime

        v_k = self.S * self.e_1 / (self.u_0 * self.phi * self.m)
        gamma = B_1 * self.psi_0 / K_1**2

        # x = Z - Z_0
        def _v(x):
            return v_k * x

        def _psi(x):
            """
            In effect:
            ψ(x) = ψ_0 + χ * σ_0 * x + λ * χ * x**2
            ψ(x) = ψ_0 + K_1 * x + λ * χ * x**2
            """
            return self.psi_0 + K_1 * x + labda_prime * chi_prime * x**2

        def _l_psi_avg(x_i, x_j):
            psi_i, psi_j = _psi(x_i), _psi(x_j)
            psi_avg = (psi_i + psi_j) / 2
            return self.l_0 * (
                1
                - self.Delta / self.rho_p
                - self.Delta * (self.alpha - 1 / self.rho_p) * psi_avg
            )

        def _Z_x(x):
            # distinct from Z which should be straightforwardly related to x
            beta = B_1 / K_1 * x
            b = (1 + 4 * gamma) ** 0.5
            return (1 - 2 * beta / (b + 1)) ** (0.5 * (b + 1) / b) * (
                1 + 2 * beta / (b - 1)
            ) ** (0.5 * (b - 1) / b)

        B_1_prime = -B_1
        gamma_prime = B_1_prime * self.psi_0 / K_1**2

        def _Z_prime_x(x):
            beta_prime = B_1_prime / K_1 * x
            b_prime = (1 - 4 * gamma_prime) ** 0.5
            return (1 + 2 * beta_prime / (b_prime + 1)) ** (
                0.5 * (b_prime + 1) / b_prime
            ) * (1 - 2 * beta_prime / (b_prime - 1)) ** (
                0.5 * (b_prime - 1) / b_prime
            )

        def _l(x_i, x_j, l_i):
            """
            given a bunch of parameters, calculate l_j
            """
            l_psi_avg = _l_psi_avg(x_i, x_j)
            if B_1 > 0:
                Z_x_i = _Z_x(x_i)
                Z_x_j = _Z_x(x_j)

                return (l_i + l_psi_avg) * (Z_x_j / Z_x_i) ** (
                    -B_0 / B_1
                ) - l_psi_avg

            elif B_1 == 0:
                xi_prime = K_1 / self.psi_0 * x_i
                xj_prime = K_1 / self.psi_0 * x_j

                return (
                    exp(xj_prime - xi_prime) * (1 + xi_prime) / (1 + xj_prime)
                ) ** (B * self.psi_0 / K_1**2) * (l_i + l_psi_avg) - l_psi_avg

            else:
                Z_prime_i = _Z_prime_x(x_i)
                Z_prime_j = _Z_prime_x(x_j)

                return (l_i + l_psi_avg) * (Z_prime_j / Z_prime_i) ** (
                    -B_0 / B_1
                ) - l_psi_avg

        def propagate(x, it):
            # propagate l from x= 0 to x_k (1-Z_0)
            # simulatneousely, also propagate a time.
            if x == 0:
                return 0, 0
            l_i = 0
            t = 0
            step = x / it
            x_i = 0
            for i in range(it):  # 8192
                x_j = x_i + step
                l_j = _l(x_i, x_j, l_i)
                Delta_l = l_j - l_i
                t += 2 * Delta_l / (_v(x_i) + _v(x_j))
                l_i = l_j
                x_i = x_j

            return l_i, t

        def _p(x, l):
            l_psi = _l_psi_avg(x, x)
            return (
                self.f
                * self.omega
                * (self.psi_0 + K_1 * x - B_1 * x**2)
                / (self.S * (l + l_psi))
            )

        delta_1 = 1 / (self.alpha - 1 / self.rho_p)
        x_k = 1 - Z_0  # upper limit of x_m, also known as x_s

        # values reflceting fracture point
        l_k, t_k = propagate(x_k, it=it)
        # l_k = _l(x_k, _l_psi_avg(x_k, 0))

        p_k = _p(x_k, l_k)
        v_k_1 = _v(x_k)

        """
            xi valid range  (0,1)
            psi valid range (0,1)
            Z valid range (0,Z_b)
        """

        psi_k = _psi(x_k)
        xi_k = 1 / self.Z_b

        """
            In the reference, the following terms are chi_s and labda_s.
            Issue is, this definition IS NOT the same as the definition
            self.chi_s and self.labda_s.from the previous chapter.
        """

        chi_k = (psi_k / xi_k - xi_k) / (1 - xi_k)
        labda_k = 1 - 1 / chi_k

        def _psi_fracture(xi):
            """xi as in greek letter ξ
            ψ(ξ) = χ_K * ξ - λ_k * χ_K * ξ^2
            """
            return chi_k * xi * (1 - labda_k * xi)

        I_s = (self.e_1 + self.rho) / (self.u_0)
        xi_0 = Z_0 / self.Z_b

        v_kk = self.S * I_s * (xi_k - xi_0) / (self.m * self.phi)

        def _v_fracture(xi):
            return v_kk * (xi - xi_0) / (xi_k - xi_0)
            # return self.S * I_s / (self.phi * self.m) * (xi - xi_k) + v_kk

        Labda_1 = (
            l_k / self.l_0
        )  # reference is wrong, this is the implied definition

        B_2 = self.S**2 * I_s**2 / (self.f * self.omega * self.phi * self.m)
        B_2_bar = B_2 / (chi_k * labda_k)
        xi_k_bar = labda_k * xi_k

        def _Labda(xi_i, xi_j, Labda_i):  # Labda = l / l_0
            xi_i_bar = labda_k * xi_i
            xi_j_bar = labda_k * xi_j

            Labda_psi_avg = _Labda_psi_avg(xi_i, xi_j)

            r = (
                (1 - (1 + 0.5 * B_2_bar * self.theta) * xi_j_bar)
                / (1 - (1 + 0.5 * B_2_bar * self.theta) * xi_i_bar)
            ) ** (-B_2_bar / (1 + 0.5 * B_2_bar * self.theta))

            Labda_j = r * (Labda_i + Labda_psi_avg) - Labda_psi_avg
            return Labda_j

        def _Labda_psi_avg(xi_i, xi_j):
            psi_i, psi_j = _psi_fracture(xi_i), _psi_fracture(xi_j)
            psi_avg = (psi_i + psi_j) / 2

            return (
                1
                - self.Delta / self.rho_p
                - self.Delta * (self.alpha - 1 / self.rho_p) * psi_avg
            )

        def _p_fracture(xi, Labda):
            """
            Equation is wrong for reference work
            """
            psi = _psi_fracture(xi)
            # v = _v_fracture(xi)
            Labda_psi = _Labda_psi_avg(xi, xi)

            return (
                self.f
                * self.Delta
                * (psi - 0.5 * B_2 * self.theta * (xi - xi_0) ** 2)
                / (Labda + Labda_psi)
            )
            """
            return (
                self.f * self.omega * psi
                - self.theta * self.phi * self.m * v**2 * 0.5
            ) / (self.S * (Labda_psi + Labda) * self.l_0)
            """

        def propagate_fracture(xi, it):
            # propagate Labda in the fracture regime
            Labda_i = Labda_1
            t = t_k
            step = (xi - xi_k) / it
            xi_i = xi_k
            for i in range(it):
                xi_j = xi_i + step
                Labda_j = _Labda(xi_i, xi_j, Labda_i)
                Delta_Labda = Labda_j - Labda_i
                xi_i = xi_j
                t += (
                    2
                    * Delta_Labda
                    * self.l_0
                    / (_v_fracture(xi_i) + _v_fracture(xi_j))
                )

            return Labda_j, t

        # values reflecting end of burn
        Labda_2, t_k_2 = propagate_fracture(1, it=it)
        # Labda_2 = _Labda(1, _Labda_psi_avg(1, xi_k))
        l_k_2 = Labda_2 * self.l_0
        p_k_2 = _p_fracture(1, Labda_2)
        v_k_2 = _v_fracture(1)
        psi_k_2 = _psi_fracture(1)

        """
            Defining equations for post burnout, adiabatic
            expansion phase.
        """
        # l_1 is the the limit for l_psi_avg as psi -> 0
        l_1 = self.l_0 * (1 - self.alpha * self.Delta)

        def _p_adb(l):
            return p_k_2 * ((l_1 + l_k) / (l_1 + l)) ** (self.theta + 1)

        def _v_adb(l):
            return (
                self.v_j
                * (
                    1
                    - ((l_1 + l_k) / (l_1 + l)) ** self.theta
                    * (1 - (v_k_2 / self.v_j) ** 2)
                )
                ** 0.5
            )

        def _t_adb(l, it):
            t = t_k_2
            step = (l - l_k_2) / it
            l_i = l_k_2
            t = t_k_2
            for i in range(it):  # 8192
                l_j = l_i + step
                t += 2 * step / (_v_adb(l_i) + _v_adb(l_j))
                l_i = l_j
            return t

        """
        Solving for muzzle exit condition
        """

        l_e = self.l_g
        Labda_e = l_e / self.l_0

        data = []
        if l_e >= l_k_2:  # shot exit happened after burnout
            v_e = _v_adb(l_e)
            p_e = _p_adb(l_e)
            psi_e = 1
            t_e = _t_adb(l_e, it=it)

        elif l_e >= l_k:  # shot exit happened after fracture
            xi_e = bisect(
                lambda xi: propagate_fracture(xi, it=it)[0] - Labda_e,
                # lambda xi: _Labda(xi, _Labda_psi_avg(xi, xi_k)) - Labda_e,
                xi_k,
                1,
                tol=tol,
            )[0]

            _, t_e = propagate_fracture(xi_e, it=it)
            v_e = _v_fracture(xi_e)
            p_e = _p_fracture(xi_e, Labda_e)
            psi_e = _psi_fracture(xi_e)

        else:  # shot exit happend before fracture
            x_e = bisect(
                lambda x: propagate(x, it=it)[0] - l_e,
                # lambda x: _l(x, _l_psi_avg(x, 0)) - l_e,
                0,
                x_k,
                tol=tol,
            )[0]
            _, t_e = propagate(x_e, it=it)
            v_e = _v(x_e)
            p_e = _p(x_e, l_e)
            psi_e = _psi(x_e)

        def findPeak(x):
            if x < x_k:  # pre burnout
                l, t = propagate(x, it=it)
                p = _p(x, l)
                psi = _psi(x)
                v = _v(x)
            else:
                xi = x / self.Z_b
                Labda, t = propagate_fracture(xi, it=it)
                l = Labda * self.l_0
                p = _p_fracture(xi, Labda)
                psi = _psi_fracture(xi)
                v = _v_fracture(xi)

            return (t, l, psi, v, p)

        x_p_1, x_p_2 = gss(
            lambda x: findPeak(x)[4], 0, self.Z_b - Z_0, tol=tol, findMin=False
        )
        if x_p_2 == self.Z_b - Z_0:  # peak @ burnout
            t_m, l_m, psi_m, v_m, p_m = t_k_2, l_k_2, psi_k_2, v_k_2, p_k_2
        else:
            t_m, l_m, psi_m, v_m, p_m = findPeak((x_p_1 + x_p_2) * 0.5)

        data.append(("SHOT EXIT", t_e, l_e, psi_e, v_e, p_e))
        data.append(("SHOT START", 0.0, 0.0, self.psi_0, 0, self.p_0))
        data.append(("PEAK PRESSURE", t_m, l_m, psi_m, v_m, p_m))
        data.append(("FRACTURE", t_k, l_k, psi_k, v_k_1, p_k))
        data.append(("BURNOUT", t_k_2, l_k_2, psi_k_2, v_k_2, p_k_2))

        data = list(
            (tag, t, l, psi, v, p, self._T(psi, l, p))
            for (tag, t, l, psi, v, p) in data
        )

        data.sort(key=lambda x: x[1].real)

        return data

    def getBP(self, abortLength, abortVel, tol=1e-3):
        """routine, meant exclusively to integrate the ODE to the point
        where burnout has occured, irregardless of barrel length specified.
        Then, the peak-finding routine is ran to find the peak pressure
        point.
        """
        l_g_bar = self.l_g / self.l_0
        t_bar_b = None
        try:
            t_bar_b, l_bar_b, v_bar_b = RKF78(
                self._ode_Z,
                (0, 0, 0),
                self.Z_0,
                self.Z_b,
                tol=tol,
                termAbv=(None, abortLength / self.l_0, abortVel / self.v_j),
            )[0]
            if l_bar_b > abortLength / self.l_0:
                raise AbortedDueToLength(
                    "burnout cannot be found within the abort length specified"
                )
            elif v_bar_b > abortVel / self.v_j:
                raise AbortedDueToVelocity(
                    "burnout cannot be found within the abort velocity specified"
                )

        except ValueError as e:
            raise e

        def f(t_bar):
            Z, l_bar, v_bar = RKF78(
                self._ode_t,
                (self.Z_0, 0, 0),
                0,
                t_bar,
                tol=tol,
            )[0]
            return self._fp_bar(Z, l_bar, v_bar)

        # tolerance is specified a bit differently for gold section search

        t_bar_p_1, t_bar_p_2 = gss(
            f,
            0,
            t_bar_b,
            tol=tol,
            findMin=False,
        )

        t_bar_p = 0.5 * (t_bar_p_1 + t_bar_p_2)

        if t_bar_b is not None and abs(t_bar_p - t_bar_b) / t_bar_p < tol:
            t_bar_p = t_bar_b

        (Z_p, l_bar_p, v_bar_p), (Z_err, l_bar_err, v_bar_err) = RKF78(
            self._ode_t, (self.Z_0, 0, 0), 0, t_bar_p, tol=tol
        )

        p_bar_p = self._fp_bar(Z_p, l_bar_p, v_bar_p)
        p_bar_max = self._fp_bar(
            Z_p + Z_err, l_bar_p - l_bar_err, v_bar_p - v_bar_err
        )
        p_bar_min = self._fp_bar(
            Z_p - Z_err, l_bar_p + l_bar_err, v_bar_p + v_bar_err
        )

        # returns the minimum barrel length required to contain the burnup.
        # the velocity at the burnup point
        # peak pressure
        return (
            l_bar_b * self.l_0,
            v_bar_b * self.v_j,
            p_bar_p * self.f * self.Delta,
            p_bar_max * self.f * self.Delta,
            p_bar_min * self.f * self.Delta,
        )

    def findBL(self, targetVel, abortLength, tol=1e-3):
        """find the necessary barrel length to derive the requisite velocity"""

        try:
            t_bar_e, Z_e, l_bar_e = RKF78(
                self._ode_v,
                (0, self.Z_0, 0),
                0,
                targetVel / self.v_j,
                tol=tol,
                termAbv=(None, None, abortLength / self.l_0),
            )[0]
        except ValueError as e:
            raise e

        return l_bar_e * self.l_0
print("\nanalytical:")
    print(
        tabulate(
            test.analyze(it=100, tol=1e-5),
            headers=("tag", "t", "l", "phi", "v", "p", "T"),
        )
    )





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


    def calculate(self):
        self.intgRecord = []

        compo = self.compositions[self.dropProp.get()]
        # lookup dictionary using the string key to get
        # the requisite object
        geom = self.geometries[self.dropGeom.get()]

        self.tableData = []
        self.errorData = []
        self.intgRecord = []

        if self.prop is None:
            return

        if self.solve_W_Lg.get() == 1:
            try:
                constrained = Constrained(
                    caliber=float(self.calmm.get()) * 1e-3,
                    shotMass=float(self.shtkg.get()),
                    propellant=self.prop,
                    startPressure=float(self.stpMPa.get()) * 1e6,
                    dragCoe=float(self.dgc.get()) * 1e-2,
                    designPressure=float(self.pTgt.get()) * 1e6,
                    designVelocity=float(self.vTgt.get()),
                )

                if self.opt_lf.get() == 0:
                    e_1, l_g = constrained.solve(
                        loadFraction=1e-2 * float(self.ldf.get()),
                        chargeMassRatio=(
                            float(self.chgkg.get()) / float(self.shtkg.get())
                        ),
                        tol=10 ** -(int(self.accExp.get())),
                        minWeb=1e-6 * float(self.minWeb.get()),
                    )

                else:
                    lf, e_1, l_g = constrained.findMinV(
                        chargeMassRatio=(
                            float(self.chgkg.get()) / float(self.shtkg.get())
                        ),
                        tol=10 ** -(int(self.accExp.get())),
                        minWeb=1e-6 * float(self.minWeb.get()),
                    )

                    lfpercent = round(
                        lf * 100, 3 - int(floor(log10(abs(l_g * 100))))
                    )
                    self.ldf.set(lfpercent)

                webmm = round(
                    2000 * e_1, 3 - int(floor(log10(abs(2000 * e_1))))
                )
                lgmm = round(l_g * 1000, 3 - int(floor(log10(abs(l_g * 1000)))))

                # take the 3 most significant digits.
                self.arcmm.set(webmm)
                self.tblmm.set(lgmm)

            except Exception as e:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                self.errorLst.append("Exception in constrained design:")
                if self.DEBUG.get():
                    # self.errorLst.append("".join(traceback.format_exception(e)))
                    self.errorLst.append(
                        "".join(
                            traceback.format_exception(
                                exc_type, exc_value, exc_traceback
                            )
                        )
                    )
                else:
                    self.errorLst.append(str(e))

        try:
            chamberVolume = (
                float(self.chgkg.get())
                / self.prop.rho_p
                / self.prop.maxLF
                / float(self.ldf.get())
                * 100
            )
            self.cv.set(toSI(chamberVolume, useSN=True))
            self.ldp.set(round(self.prop.maxLF * float(self.ldf.get()), 1))

            self.ld.set(
                toSI(
                    self.prop.maxLF
                    * 1e-2
                    * float(self.ldf.get())
                    * self.prop.rho_p,
                    unit=None,  # "g/cm^3",
                    useSN=True,
                )
            )

            self.gun = Gun(
                caliber=float(self.calmm.get()) * 1e-3,
                shotMass=float(self.shtkg.get()),
                propellant=self.prop,
                grainSize=float(self.arcmm.get()) * 1e-3,
                chargeMass=float(self.chgkg.get()),
                chamberVolume=chamberVolume,
                startPressure=float(self.stpMPa.get()) * 1e6,
                lengthGun=float(self.tblmm.get()) * 1e-3,
                chamberExpansion=float(self.clr.get()),
                dragCoe=float(self.dgc.get()) * 1e-2,
            )

            self.lx.set(toSI(float(self.tblmm.get()) / float(self.calmm.get())))
            self.tlx.set(
                toSI(
                    (
                        float(self.tblmm.get())
                        + self.gun.l_0 * 1000 / float(self.clr.get())
                    )
                    / float(self.calmm.get())
                )
            )

            self.va.set(toSI(self.gun.v_j))

        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            self.gun = None
            self.errorLst.append("Exception when defining guns:")
            if self.DEBUG.get():
                # self.errorLst.append("".join(traceback.format_exception(e)))
                self.errorLst.append(
                    "".join(
                        traceback.format_exception(
                            exc_type, exc_value, exc_traceback
                        )
                    )
                )
            else:
                self.errorLst.append(str(e))

        if self.gun is not None:
            try:
                self.tableData, self.errorData = self.gun.integrate(
                    steps=int(self.steps.get()),
                    dom=self.dropOptn.get(),
                    tol=10 ** -(int(self.accExp.get())),
                    record=self.intgRecord,
                )

                i = [i[0] for i in self.tableData].index("SHOT EXIT")
                vg = self.tableData[i][4]
                te, be = self.gun.getEff(vg)
                self.te.set(round(te * 100, 1))
                self.be.set(round(te / self.gun.phi * 100, 1))

                i = [i[0] for i in self.tableData].index("PEAK PRESSURE")
                _, _, lp, _, _, pp, _ = self.tableData[i]
                self.ptm.set(toSI(self.gun.toPt(pp, lp)))
                self.pbm.set(toSI(self.gun.toPb(pp, lp)))

            except Exception as e:
                # for Python 3.9 and below.
                exc_type, exc_value, exc_traceback = sys.exc_info()
                self.gun = None
                self.tableData = []
                self.errorData = []
                self.errorLst.append("Exception while solving gun system:")
                if self.DEBUG.get():
                    self.errorLst.append(
                        "".join(
                            traceback.format_exception(
                                exc_type, exc_value, exc_traceback
                            )
                        )
                    )
                else:
                    self.errorLst.append(str(e))

        self.tv.delete(*self.tv.get_children())
        useSN = (False, False, False, True, False, False, True)
        units = (None, "s", "m", None, "m/s", "Pa", "K")
        tableData = dot_aligned(
            self.tableData,
            units=units,
            useSN=useSN,
        )
        errorData = dot_aligned(self.errorData, units=units, useSN=useSN)
        # negErr, posErr = arrErr(self.errorData, units=units, useSN=useSN)
        i = 0
        for row, erow in zip(tableData, errorData):
            self.tv.insert(
                "", "end", str(i), values=row, tags=(row[0], "monospace")
            )
            self.tv.insert(
                str(i),
                "end",
                str(i + 1),
                values=tuple("±" + e if e != erow[0] else e for e in erow),
                tags="error",
            )

            self.tv.move(str(i + 1), str(i), "end")

            i += 2

        self.updateError()
        self.updateFigPlot()