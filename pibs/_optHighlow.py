from gun import POINT_PEAK_AVG, POINT_PEAK_SHOT
from highlow import POINT_PEAK_HIGH, POINT_PEAK_BLEED
from math import pi, log
from num import cubic


class constrainedHighlow:
    def __init__(
        self,
        caliber,
        shotMass,
        propellant,
        chargeMass,
        burstPressure,
        startPressure,
        chambrage,
        tol,
        designHighPressure,
        designLowPressure,
        designVelocity,
        dragCoefficient=0,
        minWeb=1e-6,
        maxLength=1e3,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        control=POINT_PEAK_AVG,
        **_
    ):
        if any(
            (
                caliber <= 0,
                shotMass <= 0,
                startPressure <= 0,
                dragCoefficient < 0,
                dragCoefficient >= 1,
                chambrage < 1,
            )
        ):
            raise ValueError("Invalid parameters for constrained design")

        if any(
            (
                designHighPressure <= 0,
                designLowPressure <= designHighPressure,
                designVelocity <= 0,
            )
        ):
            raise ValueError("Invalid design constraint")

        self.propellant = propellant

        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.omega = chargeMass
        self.phi_1 = 1 / (1 - dragCoefficient)  # drag work coefficient

        self.tol = tol

        self.p_0_e = burstPressure
        self.p_0_s = startPressure

        self.p_d_h = designHighPressure
        self.p_d_l = designLowPressure
        self.v_d = designVelocity

        self.chi_k = chambrage

        self.ambientP = ambientP
        self.ambientRho = ambientRho
        self.ambientGamma = ambientGamma
        self.control = control

        self.minWeb = minWeb
        self.maxLength = maxLength

        gamma = self.theta + 1
        self.cfpr = (2 / (gamma + 1)) ** (gamma / self.theta)

    def __getattr__(self, attrName):
        if "propellant" in vars(self) and not (
            attrName.startswith("__") and attrName.endswith("__")
        ):
            try:
                return getattr(self.propellant, attrName)
            except AttributeError:
                raise AttributeError(
                    "%r object has no attribute %r"
                    % (self.__class__.__name__, attrName)
                )
        else:
            raise AttributeError

    def forward(self, lengthGun, grainSize, portArea, chamberVolume, expansionVolume):
        l_g = lengthGun
        e_1 = 0.5 * grainSize
        S_j = portArea
        V_0 = chamberVolume
        V_1 = expansionVolume

        S = self.S
        l_0, l_1 = V_0 / S, V_1 / S

        # loading condition
        m, omega = self.m, self.omega

        Delta = omega / V_0

        # propellant properties
        theta, f, rho_p, alpha = self.theta, self.f, self.rho_p, self.alpha
        maxLF, chi, labda, mu = self.maxLf, self.chi, self.labda, self.mu

        loadFraction = m / (V_0 * rho_p)
        if loadFraction > maxLF:
            raise ValueError("Specified Load Fraction Violates Geometrical Constraint")

        phi_1 = self.phi_1

        p_0_e = self.p_0_e
        p_0_s = self.p_0_s

        # design constraints
        p_d_h, p_d_l, v_d = self.p_d_h, self.p_d_l, self.v_d

        chi_k = self.chi_k
        tol = self.tol

        Labda = l_g / l_1
        cc = 1 - (1 - 1 / chi_k) * log(Labda + 1) / Labda  # chambrage correction factor

        phi = phi_1 + omega / (3 * m) * cc

        v_j = (2 * f * omega / (theta * phi * m)) ** 0.5

        v_bar_d = v_d / v_j

        if self.ambientRho != 0:
            c_a_bar = (self.ambientGamma * self.ambientP / self.ambientRho) ** 0.5 / v_j
            p_a_bar = self.ambientP / (f * Delta)
        else:
            c_a_bar = 0
            p_a_bar = 0

        if v_j < v_d:
            raise ValueError(
                "Propellant load too low to achieve design velocity. "
                + " The 2nd ballistic limit for this loading conditions is"
                + " {:.4g} m/s.".format(v_j)
            )

        psi_0 = (1 / Delta - 1 / rho_p) / (f / p_0_e + alpha - 1 / rho_p)
        Zs = cubic(a=chi * mu, b=chi * labda, c=chi, d=-psi_0)
        # pick a valid solution between 0 and 1
        Zs = tuple(
            Z for Z in Zs if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to "
                + "overcome start pressure, or has burnt to post fracture."
            )
        Z_0 = Zs[0]

        S_j_bar = S_j / S

        gamma = theta + 1
        phi_2 = 0.15

        C_A = (
            phi_2
            * (0.5 * theta * phi * m / omega) ** 0.5
            * gamma**0.5
            * (2 / (gamma + 1)) ** (0.5 * (gamma + 1) / theta)
        )
        C_B = (
            phi_2
            * (0.5 * theta * phi * m / omega) ** 0.5
            * ((2 * gamma) / theta) ** 0.5
        )

        V_bar = V_1 / V_0

        cfpr = self.cfpr
