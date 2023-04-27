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
