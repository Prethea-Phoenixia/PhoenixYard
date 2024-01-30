from highlow import Highlow

# from pso import pso
from coordesc import coordesc
from math import pi, log
from num import cubic

# from gun import
from gun import POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_EXIT
from highlow import POINT_PEAK_HIGH, POINT_PEAK_BLEED


class optHighlow:
    def __init__(
        self,
        caliber,
        propellant,
        shotMass,
        burstPressure,
        startPressure,
        dragCoefficient,
        chambrage,
        tol,
        designHighPressure,
        designLowPressure,
        designVelocity,
        minWeb=1e-6,
        maxLength=1e3,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        control=POINT_PEAK_AVG,
        **_,
    ):
        # cache the constants:
        self.S = (0.5 * caliber) ** 2 * pi
        self.propellant = propellant
        self.m = shotMass
        self.p_0_e = burstPressure
        self.p_0_s = startPressure

        self.p_d_h = designHighPressure
        self.p_d_l = designLowPressure
        self.v_d = designVelocity

        self.chi_k = chambrage

        self.phi_1 = 1 / (1 - dragCoefficient)  # drag work coefficient

        self.tol = tol

        self.maxLength = maxLength
        self.minWeb = minWeb

        self.ambientRho = ambientRho
        self.ambientP = ambientP
        self.ambientGamma = ambientGamma
        self.control = control

        theta = self.theta
        gamma = theta + 1
        phi_2 = 0.15  # for small ports 1.5mm to 2mm in size

        self.cfpr = (2 / (gamma + 1)) ** (gamma / theta)
        self.K_0 = gamma**0.5 * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * theta))

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

    # def constrained(
    #     self,
    #     minWeb,
    #     maxWeb,
    #     minLF,
    #     minEV,
    #     maxEV,
    #     minPortRatio,
    #     maxPortRatio,
    #     minLength,
    #     maxLength,
    #     control,  # targeted pressure
    #     designHighPressure,
    #     designLowPressure,
    #     designVelocity,
    # ):
    #     self.control = control

    #     self.designHighPressure = designHighPressure
    #     self.designLowPressure = designLowPressure
    #     self.designVelocity = designVelocity

    #     self.show = False

    #     if self.propellant.maxLF < minLF:
    #         raise ValueError

    #     minCV = self.shotMass / self.propellant.rho_p / self.propellant.maxLF
    #     maxCV = self.shotMass / self.propellant.rho_p / minLF

    #     S = self.caliber**2 * pi * 0.25
    #     minPortA = S * minPortRatio
    #     maxPortA = S * maxPortRatio

    #     constraints = [
    #         (minWeb, maxWeb),
    #         (minCV, maxCV),
    #         (minEV, maxEV),
    #         (minPortA, maxPortA),
    #         # (minLength, maxLength),
    #     ]
    #     self.show = True
    #     solution, dev = coordesc(
    #         self._f, constraints, x_rel_tol=self.tol, y_abs_tol=self.tol
    #     )
    #     self._f(*solution)

    #     for (s_min, s_max), p in zip(constraints, solution):
    #         print(s_min, p, s_max)

    # def _f(self, grainSize, chamberVolume, expansionVolume, portArea):
    #     highlow = Highlow(
    #         caliber=self.caliber,
    #         shotMass=self.shotMass,
    #         propellant=self.propellant,
    #         grainSize=grainSize,
    #         chargeMass=self.chargeMass,
    #         chamberVolume=chamberVolume,
    #         expansionVolume=expansionVolume,
    #         burstPressure=self.burstPressure,
    #         startPressure=self.startPressure,
    #         portArea=portArea,
    #         chambrage=self.chambrage,
    #         lengthGun=10,
    #         dragCoefficient=self.dragCoefficient,
    #     )

    #     data, _, _, _ = highlow.integrate(
    #         step=0,
    #         tol=self.tol,
    #         ambientRho=self.ambientRho,
    #         ambientP=self.ambientP,
    #         ambientGamma=self.ambientGamma,
    #         peaks=[self.control, POINT_PEAK_HIGH],
    #     )

    #     tags = [line[0] for line in data]

    #     i_high = tags.index(POINT_PEAK_HIGH)
    #     i_low = tags.index(self.control)
    #     i_vel = tags.index(POINT_EXIT)

    #     p_high = data[i_high][5]

    #     if self.control == POINT_PEAK_BLEED:
    #         p_low = data[i_low][6]
    #     elif self.control == POINT_PEAK_AVG:
    #         p_low = data[i_low][7]
    #     elif self.control == POINT_PEAK_SHOT:
    #         p_low = data[i_low][8]

    #     v_g = data[i_vel][4]

    #     if self.show:
    #         print(f"p_high={p_high}, p_low={p_low}, v={v_g}")

    #     dev = sum(
    #         v**2
    #         for v in (
    #             (p_high - self.designHighPressure) / self.designHighPressure,
    #             (p_low - self.designLowPressure) / self.designLowPressure,
    #             # (v_g - self.designVelocity) / self.designVelocity,
    #         )
    #     )

    #     return dev

    def solve(
        self,
        loadFraction,
        chargeMassRatio,
        lengthGun=None,
        known_bore=False,
        suppress=False,
        progressQueue=None,
    ):
        if any((chargeMassRatio <= 0, loadFraction <= 0, loadFraction > 1)):
            raise ValueError("Invalid parameters to solve constrained design problem")

        p_max = 1e9  # 1GPa

        m = self.m
        rho_p = self.rho_p
        theta = self.theta
        f = self.f
        chi = self.chi
        chi_k = self.chi_k
        labda = self.labda
        mu = self.mu
        S = self.S
        maxLF = self.maxLF
        phi_1 = self.phi_1
        p_0 = self.p_0
        v_d = self.v_d
        p_d_l = self.p_d_l
        p_d_h = self.p_d_h

        u_1 = self.u_1
        n = self.n
        alpha = self.alpha
        Z_b = self.Z_b

        f_psi_Z = self.f_psi_Z
        f_sigma_Z = self.f_sigma_Z

        p_0_e = self.p_0_e
        p_0_s = self.p_0_s

        cfpr = self.cfpr
        K_0 = self.K_0

        Z_b = self.Z_b
        Z_0 = self.Z_0s

        p_a = self.ambientP
        if self.ambientRho != 0:
            c_a = (self.ambientGamma * self.ambientP / self.ambientRho) ** 0.5
        else:
            c_a = 0
        k_1 = self.ambientGamma

        Sb = S * chi_k

        if loadFraction > maxLF:
            raise ValueError("Specified Load Fraction Violates Geometrical Constraint")

        omega = m * chargeMassRatio
        V_0 = omega / (rho_p * loadFraction)
        Delta = omega / V_0
        l_0 = V_0 / S

        def _f_p_1(self, Z, eta, tau_1, psi=None):
            psi = psi if psi else f_psi_Z(Z)
            V_psi = V_0 - omega / rho_p * (1 - psi) - alpha * omega * (psi - eta)
            return f * omega * tau_1 / V_psi * (psi - eta)

        def func_v(lengthGun):
            def func_p(expansionVolume, portAreaRatio):
                V_1 = expansionVolume

                def _f_p_2(self, l, eta, tau_2):
                    l_star = (V_1 - alpha * omega * eta) / S
                    return f * omega * tau_2 * eta / (S * (l_star + l))

                l_1 = V_1 / S
                S_j_bar = portAreaRatio
                S_j = S * S_j_bar
                l_g = lengthGun
                Labda = l_g / l_1
                cc = 1 - (1 - 1 / chi_k) * log(Labda + 1) / Labda
                phi = phi_1 + omega / (3 * m) * cc

                v_j = (2 * f * omega / (theta * phi * m)) ** 0.5

                psi_0 = (
                    p_0_e
                    * (V_0 - omega / rho_p)
                    / (f * omega - p_0_e * (omega / rho_p - alpha * omega))
                )

                Zs = cubic(chi * mu, chi * labda, chi, -psi_0)
                # pick a valid solution between 0 and 1
                Zs = sorted(
                    Z
                    for Z in Zs
                    if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
                )  # evaluated from left to right, guards against complex >/< float
                if len(Zs) < 1:
                    raise ValueError(
                        "Propellant either could not develop enough pressure to overcome"
                        + " port open pressure, or has burnt to post fracture."
                    )
                Z_0 = Zs[0]

                t_0, eta_0, tau_1_0, tau_2_0 = 0, 0, 1, 1 + theta


if __name__ == "__main__":
    from tabulate import tabulate
    from math import pi
    from prop import GrainComp, Propellant

    compositions = GrainComp.readFile("data/propellants.csv")

    M17 = compositions["M17"]
    M1 = compositions["M1"]
    from prop import SimpleGeometry

    # M17C = Propellant(M17, SimpleGeometry.CYLINDER, None, 2.5)
    # M1C = Propellant(M1, SimpleGeometry.CYLINDER, None, 10)
    # test = optHighlow(
    #     caliber=0.082,
    #     propellant=M1C,
    #     shotMass=5,
    #     chargeMass=0.3,
    #     burstPressure=10e6,
    #     startPressure=30e6,
    #     dragCoefficient=3e-3,
    #     chambrage=1.5,
    #     tol=1e-3,
    # )

    # result = test.constrained(
    #     minWeb=1e-6,
    #     maxWeb=1e-2,
    #     minLF=0.1,
    #     minPortRatio=0.1,
    #     maxPortRatio=1 / 0.15,
    #     minLength=0.01,
    #     maxLength=1,
    #     minEV=0,
    #     maxEV=2e-3,
    #     control=POINT_PEAK_AVG,  # targeted pressure
    #     designHighPressure=80e6,
    #     designLowPressure=40e6,
    #     designVelocity=75,
    # )
