"""
Depreceated code that i spent way too much time on and can't be parted
with. 
"""


def analyze(self, tol=1e-5):
    """run the psi-bar analytical solution on the defined weapon

    key assumptions:
        >burn rate scales linearlly with pressure
        >propellant volume fraction (ψ, psi) is a polynomial
         function of linear burn fraction (Z),
         (this involves discarding the 3rd order term entire-
         ly, which is not exactly rigorous, but is necessary
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


    This basically doesn't work since its awkwardly not accurate enough
    to be useful and complicated enough that its too slow for rough
    analysis or optimisation. USELESS.
    """

    def _fpsi(Z):
        return self.chi * Z * (1 + self.labda * Z + self.mu * Z**2)

    Z_0 = bisect(lambda Z: _fpsi(Z) - self.psi_0, 0, 1, tol=tol)[0]

    psi_s = _fpsi(1)
    labda_prime = (psi_s * Z_0 - self.psi_0) / (self.psi_0 - Z_0**2 * psi_s)
    chi_prime = psi_s / (1 + labda_prime)

    # labda_prime = self.labda
    # chi_prime = self.chi

    # the equivalent of B used for linear burn rates.
    B_0 = (
        self.S**2
        * self.e_1**2
        / (self.f * self.phi * self.omega * self.m * self.u_0**2)
    )

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

    def _l(x, l_psi_avg):
        if B_1 > 0:
            return l_psi_avg * (_Z_x(x) ** (-B_0 / B_1) - 1)
        elif B_1 == 0:
            """this case is extremely rare and has not been rigorously
            tested."""
            x_prime = K_1 / self.psi_0 * x
            return (
                exp(B_0 * self.psi_0 / K_1**2 * (x_prime - log(1 + x_prime)))
                - 1
            ) * l_psi_avg
        else:  # this negation was corrected from the original work
            return l_psi_avg * (_Z_prime_x(x) ** (-B_0 / B_1) - 1)

    def propagate(x, tol):
        # propagate l from x= 0 to x_k (1-Z_0) using global adaptive step control
        # simulatneousely, also propagate a time.
        if x == 0:
            return 0, 0
        l = 0
        t = 0
        for i in range(13):  # 8192
            step = x / 2**i
            x_i = 0
            l_i = 0
            t_i = 0
            for j in range(2**i):
                x_j = step * (j + 1)
                l_psi_avg = _l_psi_avg(x_i, x_j)
                Delta_l = _l(x_j, l_psi_avg) - _l(x_i, l_psi_avg)
                l_i += Delta_l
                t_i += 2 * Delta_l / (_v(x_i) + _v(x_j))
                if isinstance(l_i, complex):
                    # complex result! abort now and use a finer step
                    raise ValueError("Analytical model diverged.")
                x_i = x_j
            if abs(l_i - l) > tol * self.l_0 and abs(t_i - t) > tol * (
                self.l_0 / self.v_j
            ):  # change in iteration
                t = t_i
                l = l_i
                continue
            else:
                return l_i, t_i
        raise ValueError(
            "Unable to propagate system in the pre-fracture regime,"
            + " to required accuracy "
            + "in reasonable cycles. If this error remains despite"
            + " reduced accuracy specification, it is likely that "
            + "the analytical model has already diverged."
        )

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
    for i in range(10):  # 8192
        x_m_i = _x_m(p_m_i)

        if x_m_i > x_k:
            x_m_i = x_k
        elif x_m_i < 0:
            x_m_i = 0

        l_m_i, t_m_i = propagate(x_m_i, tol=tol)  # see above
        # l_m_i = _l(x_m_i, _l_psi_avg(x_m_i, 0))
        p_m_j = _p(x_m_i, l_m_i)

        if abs(p_m_i - p_m_j) > self.f * self.Delta * tol:
            p_m_i = p_m_j
        else:
            break
    if i == 10:
        raise ValueError(
            "Unable to iteratively solve for peak pressure"
            + " within reasonable cycles"
        )

    # value reflecting peak pressure point
    p_m = p_m_j
    x_m = x_m_i
    l_m = l_m_i
    t_m = t_m_i
    v_m = _v(x_m)
    psi_m = _psi(x_m)

    # values reflceting fracture point
    l_k, t_k = propagate(x_k, tol=tol / 3)
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

    Labda_1 = (
        l_k / self.l_0
    )  # reference is wrong, this is the implied definition

    B_2 = self.S**2 * I_s**2 / (self.f * self.omega * self.phi * self.m)
    B_2_bar = B_2 / (chi_k * labda_k)
    xi_k_bar = labda_k * xi_k

    def _Labda(xi, Labda_psi_avg):  # Labda = l / l_0
        xi_bar = labda_k * xi

        r = (
            (1 - (1 + 0.5 * B_2_bar * self.theta) * xi_bar)
            / (1 - (1 + 0.5 * B_2_bar * self.theta) * xi_k_bar)
        ) ** (-B_2_bar / (1 + 0.5 * B_2_bar * self.theta))

        Labda = r * (Labda_1 + Labda_psi_avg) - Labda_psi_avg
        return Labda

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
        v = _v_fracture(xi)
        Labda_psi = _Labda_psi_avg(xi, xi)

        return (
            self.f
            * self.Delta
            * (
                psi
                - self.theta
                * self.phi
                * self.m
                * v**2
                / (2 * self.f * self.omega)
            )
            / (Labda + Labda_psi)
        )

    def propagate_fracture(xi, tol):
        # propagate Labda in the fracture regime using global adaptive step control
        Labda = Labda_1
        t = t_k
        for i in range(13):
            step = (xi - xi_k) / 2**i
            xi_i = xi_k
            Labda_i = Labda_1
            t_i = t_k
            for j in range(2**i):
                xi_j = xi_i + step
                Labda_psi_avg = _Labda_psi_avg(xi_i, xi_j)
                Delta_Labda = _Labda(xi_j, Labda_psi_avg) - _Labda(
                    xi_i, Labda_psi_avg
                )
                Labda_i += Delta_Labda
                if isinstance(Labda_i, complex):
                    raise ValueError("Analytical model diverged.")
                xi_i = xi_j
                t_i += (
                    2
                    * Delta_Labda
                    * self.l_0
                    / (_v_fracture(xi_i) + _v_fracture(xi_j))
                )

            if abs(Labda_i - Labda) > tol and abs(t_i - t) > tol * (
                self.l_0 / self.v_j
            ):
                # change in iteration
                Labda = Labda_i
                t = t_i
                continue
            else:
                return Labda_i, t_i
        raise ValueError(
            "Unable to propagate system in the post-fracture regime,"
            + " to required accuracy in reasonable cycles."
            + " If this error remains despite"
            + " reduced accuracy specification, it is likely that "
            + "the analytical model has already diverged."
        )

    # values reflecting end of burn
    Labda_2, t_k_2 = propagate_fracture(1, tol=tol / 3)
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

    def _t_adb(l, tol):
        t = t_k_2
        for i in range(13):  # 8192
            step = (l - l_k_2) / 2**i
            l_i = l_k_2
            t_i = t_k_2
            for j in range(2**i):
                l_j = l_i + step
                t_i += 2 * step / (_v_adb(l_i) + _v_adb(l_j))
                if isinstance(l_i, complex):
                    # complex result! abort now and use a finer step
                    raise ValueError("Analytical model diverged.")
                l_i = l_j
            if abs(t_i - t) > tol * self.l_0 / self.v_j:  # change in iteration
                t = t_i
                continue
            else:
                return t_i
        raise ValueError(
            "Unable to propagate system in the adiabatic expansion regime,"
            + " to required accuracy in reasonable cycles."
            + " If this error remains despite"
            + " reduced accuracy specification, it is likely that "
            + "the analytical model has already diverged."
        )

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
        t_e = _t_adb(l_e, tol=tol / 3)

    elif l_e >= l_k:  # shot exit happened after fracture
        xi_e = bisect(
            lambda xi: propagate_fracture(xi, tol=tol / 3)[0] - Labda_e,
            # lambda xi: _Labda(xi, _Labda_psi_avg(xi, xi_k)) - Labda_e,
            xi_k,
            1,
            tol=tol,
        )[0]

        _, t_e = propagate_fracture(xi_e, tol=tol / 3)
        v_e = _v_fracture(xi_e)
        p_e = _p_fracture(xi_e, Labda_e)
        psi_e = _psi_fracture(xi_e)

    else:  # shot exit happend before fracture
        x_e = bisect(
            lambda x: propagate(x, tol=tol / 3)[0] - l_e,
            # lambda x: _l(x, _l_psi_avg(x, 0)) - l_e,
            0,
            x_k,
            tol=tol,
        )[0]
        _, t_e = propagate(x_e, tol=tol / 3)
        v_e = _v(x_e)
        p_e = _p(x_e, l_e)
        psi_e = _psi(x_e)

    data.append(("SHOT EXIT", t_e, l_e, psi_e, v_e, p_e))
    data.append(("SHOT START", 0.0, 0.0, self.psi_0, 0, self.p_0))
    data.append(("PEAK PRESSURE", t_m, l_m, psi_m, v_m, p_m))
    data.append(("FRACTURE", t_k, l_k, psi_k, v_k_1, p_k))
    data.append(("BURNOUT", t_k_2, l_k_2, psi_k_2, v_k_2, p_k_2))

    data.sort(key=lambda x: x[1])

    return data

    # values


print("\nanalytical")
print(tabulate(test.analyze(1e-3), headers=("tag", "t", "l", "phi", "v", "p")))



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

        print(
            p_1
            / (B * (1 + self.theta))
            * ((1 + 0.5 * self.theta) / (1 + self.theta))
            ** ((2 + self.theta) / self.theta)
        )

        psi_m = 1 / (1 + B * self.theta)
        print(psi_m)
        l_m = _l1_psi(psi_m)
        v_m = _v1_psi(psi_m)
        p_m = _p1_psi(psi_m)

        for i in range(100):
            print((i + 1) / 10, _p1_psi((i + 1) / 10))
            print(_v1_psi((i + 1) / 10), _v1_l(_l1_psi((i + 1) / 10)))
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
