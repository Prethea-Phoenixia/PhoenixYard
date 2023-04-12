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
