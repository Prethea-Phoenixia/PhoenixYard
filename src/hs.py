"""
Hirschfelder-Sherman calculation
additive approximation for propellant thermochemical properties,
by tallying propellant heat of explosion and calculating the
deviaiton to the reference temperature of 2500K, assuming 
constant heat capacity between 2500K and 3000K, propellant isochoric
adiabatic flame temperature is then calculated. 

Combustion product species are assumed to involve CO, CO2, H2O, H2,
N2, HCl, and linear equations are solved to estimate product gas 
volume (the inverse of molecular weight), and adiabatic ratios.

Results are within few (up to 3) percent of accurate thermalchemical
calculation in the range of 2000K to 4000K. Should not be applied to
cases >4000K where dissociation of of free radicals like -H, -OH and
-Cl would happen. Should not be applied to propellants generating 
substantial amount of consdensed exhaust.

Significant amount of propellant research nowadays are
published with incomplete data to be used in the inetrnal
ballistic formulation that is present in this program. Ideally
accurate thermochemcial programs like BLAKE should be consulted.
In the absence of more sophisticated programs, the Hirschfelder
-Sherman is an approximate, and simple method appropriate for
propellant screening and program guidance purposes.

Bd-MVP copolymer is a 90% butadiene; 10% 2-methyl-5-vinylpyridine copolymer

Metriol trinitrate,MTN is alternatively known as Trimethylolethane trinitrate,
TMETN, two entries are used to represent the two possible naming scheme.

Picrite is an alternative name for Nitroguanidine.

Some entries are removed due to inherent uncertainty of composition
Asphalt                      ,            ,0.2179 ,-2305  ,0.09450,  ,  , ,  ,        
Bd-MVP copolymer             ,            ,0.4132 ,-3183  ,0.11544,  ,  , ,  ,  
Polyester                    ,            ,0.3552 ,-2620  ,0.09123,  ,  , ,  ,  
Polyurethane,PU,0.4073,-3773,0.10796,27,36,2,10,


Other entries are removed due to impossibility to research what material this is.
Petrin                       ,            ,0.3703 ,374    ,0.04109,  ,  , ,  ,  

NC entries are removed, since these can be programmatically generate 
Nitrocellulose 12.20% N      ,NC1220      ,0.3478 ,137.7  ,0.04127,  ,  , ,  ,        
Nitrocellulose 12.60% N      ,NC1260      ,0.3454 ,198.9  ,0.04040,  ,  , ,  ,        
Nitrocellulose 13.15% N      ,NC1315      ,0.3421 ,283.1  ,0.03920,  ,  , ,  , 

These entries are removed due to large mismatch between values computed from 
composition and tabulated ones.
Metriol trinitrate,MTN,0.3052,377,0.04313,5,9,3,9,
Trimethylolethane trinitrate,TMETN,0.3052,377,0.04313,5,9,3,9,


"""
import csv
import difflib
from num import quadratic
from math import exp, log10


class Ingredient:
    """
    Ingredient class

    """

    allDict = {}
    altDict = {}

    def __init__(self, name, alt, Cvi, Ei, invMi, C, H, N, O, Ext):
        self.name = name
        self.alt = alt
        self.Cvi = Cvi
        self.Ei = Ei
        self.invMi = invMi

        self.C = C
        self.H = H
        self.N = N
        self.O = O
        self.Ext = Ext

        A = 12.01 * C + 1.008 * H + 14.008 * N + 16.00 * O + Ext  # g/mol
        self.A = A

    @classmethod
    def readFile(cls, fileName):
        ingredients = []
        with open(fileName, newline="") as csvfile:
            ingrReader = csv.reader(
                csvfile,
                delimiter=",",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )

            skipFirstLine = True

            for ingr in ingrReader:
                if skipFirstLine:
                    skipFirstLine = False
                    continue
                (name, alt, Cvi, Ei, invMi, C, H, N, O, Ext) = ingr

                C, H, N, O, Ext = (
                    float(v) if v != "" else 0 for v in (C, H, N, O, Ext)
                )

                newIngr = cls(
                    name=name,
                    alt=alt,
                    Cvi=float(Cvi),
                    Ei=float(Ei),
                    invMi=float(invMi),
                    C=C,
                    H=H,
                    N=N,
                    O=O,
                    Ext=Ext,
                )

                ingredients.append(newIngr)

        ingrDict = {i.name: i for i in ingredients}
        cls.allDict.update(ingrDict)
        cls.altDict.update({i.alt: i for i in ingredients if i.alt != ""})
        return ingrDict

    @classmethod
    def fromElement(
        cls,
        name,
        C=0,
        H=0,
        N=0,
        O=0,
        HoC=0,
        # Q=None,
        alt="",
        u="kJ/mol",
        keep=True,
    ):
        """
        Given the molecular formula of a chemical, estimate factors necessary for
        use in the Hirschfelder-Sherman calculation, and add the newly created
        ingredient into the class.
        """
        # accurate molecular mass here to account for natural abundance of isotopes
        A = 12.01 * C + 1.008 * H + 14.008 * N + 16.00 * O  # g/mol

        """
        convert element nbr. (mol/mol) into nbr. mol of element per unit mass (mol/g)
        """
        Ci = C / A
        Hi = H / A
        Ni = N / A
        Oi = O / A

        # if HoC is not None:
        if u == "kJ/mol":
            HoC /= 4.184  # to kcal/mol
            HoC /= A  # to kcal/g
            HoC *= 1000  # to cal/g
        elif u == "kJ/kg":
            Q /= 4.184  # to kcal/kg
        elif u == "kcal/mol":
            HoC /= A
            HoC *= 1000
        elif u == "kcal/kg" or u == "cal/g":
            pass  # kcal/kg = cal/g
        else:
            raise ValueError("Unknown unit ", u)

        invM = Ci + 0.5 * Hi + 0.5 * Ni
        # isochoric heat capacity from 2000-3000K
        Cv = 1.620 * Ci + 3.265 * Hi + 3.384 * Ni + 5.193 * Oi
        # HoC: heat of combustion, +: energy is released, -: energy is consumed
        # this is the opposite of the usual convention of enthalpy of combustion.
        E = HoC - 132771 * Ci - 40026 * Hi - 6724 * Ni + 51819 * Oi
        # heat required to bring combustion product to 2500k, +: product hotter than 2.5kK
        # -: product cooler than 2.5kK

        newIngr = cls(name, alt, Cv, E, invM, Ci, Hi, Ni, Oi, 0)
        if keep:
            cls.allDict.update({newIngr.name: newIngr})
            if newIngr.alt != "":
                cls.altDict.update({newIngr.alt: newIngr})

        return newIngr

    @classmethod
    def nitrocellulose(cls, Nfraction, keep=True):
        """
        Initialize nitrocellulose ingredient with arbitrary nitrogen mass fraction
        Ref: table 2.09 from Hunter, 1951, values updated with AMCP 706-175.
            When values conflict, the AMCP version is adopted.

        If the an entry at the same nitration level is found, the keep parameter is
        ignored and the class dictionary is not updated.

        Typically N level is between 11.50% and 13.50%
        """
        f = Nfraction
        x = 162.14 * f / (14.008 - 45 * f)
        C = 6
        H = 10 - x
        O = 5 + 2 * x
        N = x

        name = "Nitrocellulose {:.2%} N".format(f)
        alt = "NC{:d}".format(int(f * 10000))

        Cv = 0.3421 + 0.006 * (13.15 - f * 100)
        E = 283.1 - 153 * (13.15 - f * 100)
        invM = 0.03920 + 0.00218 * (13.15 - f * 100)

        newIngr = cls(
            name=name,
            alt=alt,
            Cvi=Cv,
            Ei=E,
            invMi=invM,
            C=C,
            H=H,
            O=O,
            N=N,
            Ext=0,
        )
        if keep and name not in cls.allDict:
            cls.allDict.update({newIngr.name: newIngr})
            cls.altDict.update({newIngr.alt: newIngr})

        return newIngr

    @classmethod
    def find(cls, name):
        """
        return an Ingredient object that is closest to the supplied common name,
        or alternative (abbreviated) name.
        """
        closeIngrs = difflib.get_close_matches(
            name,
            list(cls.allDict.keys()) + list(cls.altDict.keys()),
            n=3,
            cutoff=0.6,
        )
        n = 0
        thisIngr = None
        ingrs = []

        if len(closeIngrs) > 0:
            print("Found candidate:")
            for iname in closeIngrs:
                if iname in cls.allDict:
                    ingr = cls.allDict[iname]
                else:
                    ingr = cls.altDict[iname]
                if ingr not in ingrs:
                    ingrs.append(ingr)
                if n == 0:
                    thisIngr = ingr
                n += 1

            for ingr in ingrs:
                print("-" + str(ingr))

            print("returning " + str(thisIngr) + "\n")
            return thisIngr
        else:
            print('Unknown ingredient description "{:}"'.format(name) + "\n")
            return None

    @classmethod
    def check(cls):
        print(
            "{:_^30}|{:_^15}|{:_>5}{:_>5}{:_>5}{:_>5}|{:_>10}{:_>10}{:_>10}|{:_>10}{:_>10}|".format(
                "Name",
                "Alt",
                "C",
                "H",
                "N",
                "O",
                "Cvi",
                "Ei",
                "ni",
                "%Cvi",
                "%ni",
            )
        )
        for ingr in cls.allDict.values():
            A = ingr.A

            Ci = ingr.C / A
            Hi = ingr.H / A
            Ni = ingr.N / A
            Oi = ingr.O / A

            invM = Ci + 0.5 * Hi + 0.5 * Ni
            # isochoric heat capacity from 2000-3000K
            Cv = 1.620 * Ci + 3.265 * Hi + 3.384 * Ni + 5.193 * Oi
            # HoC: heat of combustion, +: energy is released, -: energy is consumed
            # this is the opposite of the usual convention of enthalpy of combustion.

            print(
                "{:^30}|{:^15}|{:5.3g}{:5.3g}{:5.3g}{:5.3g}|{:10.4g}{:10.4g}{:10.4g}|{:10.1%}{:10.1%}|".format(
                    ingr.name,
                    ingr.alt,
                    ingr.C,
                    ingr.H,
                    ingr.N,
                    ingr.O,
                    ingr.Cvi,
                    ingr.Ei,
                    ingr.invMi,
                    abs(ingr.Cvi - Cv) / ingr.Cvi,
                    abs(ingr.invMi - invM) / ingr.invMi,
                )
            )

    def prettyPrint(self):
        print(str(self))

        print("Hirschfelder-Sherman factors:")
        print("{:>6}  {:>7}  {:>7}".format("Cvi", "Ei", "1/Mi"))
        print(
            "{:>6.4f}  {:>7.1f}  {:>7.5f}".format(self.Cvi, self.Ei, self.invMi)
        )

        print("Chemical Composition:")
        print(
            "C{:.1f} H{:.1f} N{:.1f} O{:.1f}, A={:.1f} g/mol".format(
                self.C, self.H, self.N, self.O, self.A
            )
        )

        # print("Estimated Covolume: {:} m^3/kg".format(self.b))

    def __str__(self):
        if self.alt != "":
            return "{:} ({:})".format(self.name, self.alt)
        else:
            return self.name


class Mixture:
    def __init__(self, name, compoDict, Delta=0.2):
        self.name = name
        self.Delta = Delta
        """
        Normalize the given composition such that the fractions sums to 1
        """
        total = 0
        for ingr, fraction in compoDict.items():
            total += fraction

        self.compoDict = {
            ingr: fraction / total for ingr, fraction in compoDict.items()
        }

        """
        tally the releavnt factors according to their mass fraction
        """
        Cv = 0  # heat capacity at constant volume, presumably in cal/g K
        E = 0  # heat of formation, also presumably in cal/g
        invM = 0  # specific gas volume, mol/g

        Ci, Hi, Ni, Oi = 0, 0, 0, 0

        for ingr, fraction in self.compoDict.items():
            Cv += fraction * ingr.Cvi
            E += fraction * ingr.Ei
            invM += fraction * ingr.invMi

            Ci += fraction * ingr.C / ingr.A
            Hi += fraction * ingr.H / ingr.A
            Ni += fraction * ingr.N / ingr.A
            Oi += fraction * ingr.O / ingr.A

        M = 1 / invM  # g/mol
        Tv = 2500 + E / Cv
        if Tv > 3000:
            Tv = 3000 + 6046 * (
                -Cv
                - 0.01185
                + ((Cv + 0.01185) ** 2 + 3.308e-4 * (E - 500 * Cv)) ** 0.5
            )  # in Kelvin

        gamma = 1 + 1.987 / (Cv * M)
        f = 8.314 * Tv / M * 1e3  # J/kg
        Cp = Cv * gamma

        self.Tv = Tv
        self.gamma = gamma
        self.f = f

        self.speciesList, self.b, self.p = balance(
            Tv, Ci, Hi, Oi, Ni, V=1 / Delta
        )
        """
        covolume estimate suggested by Cook:
        self.b = 1e-3 * (1.18 + 6.9 * Ci - 11.5 * Oi)
        
        linear regression fit from original data:
        self.b = 1e-3 * (
            1.3372800285521016 + 15.00623286 * Ci - 20.81820076 * Oi
        )
        linear regression fit from original data using all elements:
        self.b = 1e-3 * (
            -15.125770768899
            + 202.80374114 * Ci
            + 27.53266578 * Hi
            + 222.42131159 * Ni
            + 242.02396858 * Oi
        )

        These aren't good enough, esp. not for the modern propellants like JA2,
        the deviaiton is up to 0.2 (0.89 vs 0.98)
        """
        self.C = Ci * 12.01
        self.H = Hi * 1.008
        self.O = Oi * 16.00
        self.N = Ni * 14.008

    def prettyPrint(self):
        print("Mixture: {:}".format(self.name))
        print("Specified Composition:---------------------------")
        for ingr, fraction in self.compoDict.items():
            print("--{:-<30}, {:>6.2%}".format(str(ingr), fraction))

        print("")
        print("Elemental Fractions:-----------------------------")
        print(
            "C {:.2%} H {:.2%} N {:.2%} O {:.2%}".format(
                self.C, self.H, self.N, self.O
            )
        )

        print("")
        print("Hirschfelder-Sherman Estimations:----------------")
        print("Isochoric Adb. Temp: {:>6.1f}K".format(self.Tv))
        print("Adiabatic Index    : {:>6.3f}".format(self.gamma))
        print("Specific Force     : {:>6.4f} MJ/kg".format(self.f / 1e6))
        print("")

        print("Species--%wt.---%mol-----------------------------")
        print(
            *("{:^5}  {:>6.2%}  {:>6.2%}".format(*v) for v in self.speciesList),
            sep="\n"
        )
        print("")
        print("Further thermalchemical Calculations:------------")

        print("Covolume           : {:>6.4g} cc/g".format(self.b))
        print(" @Temperature      : {:>6.1f} K".format(self.Tv))
        print(" @Load Density     : {:>6.3g} g/cc".format(self.Delta))
        print(" @Pressure         : {:>6.4g} MPa".format(self.p))
        print("")


"""
TABLE 2.07 from Hunt
T in Kelvin, -DeltaB in cc/(gm.mol), -DeltaC/2 in (cc/(gm.mol))^2
"""

negDeltaTable = [
    [1600, 34.2, 490],
    [1700, 33.8, 460],
    [1800, 33.5, 435],
    [1900, 33.2, 410],
    [2000, 33.0, 390],
    [2100, 32.8, 370],
    [2200, 32.6, 355],
    [2300, 32.5, 340],
    [2400, 32.4, 325],
    [2500, 32.2, 310],
    [2600, 32.1, 300],
    [2700, 32.0, 290],
    [2800, 31.9, 280],
    [2900, 31.8, 270],
    [3000, 31.7, 260],
    [3100, 31.6, 255],
    [3200, 31.5, 245],
    [3300, 31.4, 235],
    [3400, 31.4, 230],
    [3500, 31.3, 225],
    [3600, 31.2, 215],
    [3700, 31.1, 210],
    [3800, 31.1, 205],
    [3900, 31.0, 200],
    [4000, 30.9, 195],
]
"""
TABLE 2.04 from Hunt
T in Kelvin, K0(T)..... K6(T)
"""
equilibriumKT = [
    [1000, 0.7185, 3e-11, 1e-14, 8e-21, 2e-10, 3e-9, 5e-16],
    [1200, 1.406, 7e-9, 1e-11, 2e-16, 3e-8, 2e-7, 7e-13],
    [1400, 2.212, 3e-7, 1e-9, 2e-13, 1e-6, 5e-6, 1e-10],
    [1600, 3.043, 6e-6, 5e-8, 5e-11, 2e-5, 6e-5, 5e-9],
    [1700, 3.438, 2e-5, 2e-7, 5e-10, 6e-5, 0.0002, 2e-8],
    [1800, 3.832, 6e-5, 8e-7, 3e-9, 0.0001, 0.0004, 1e-7],
    [1900, 4.206, 0.0002, 3e-6, 2e-8, 0.0003, 0.0008, 4e-7],
    [2000, 4.574, 0.0004, 8e-6, 9e-8, 0.0007, 0.0017, 1e-6],
    [2100, 4.909, 0.0009, 2e-5, 4e-7, 0.0015, 0.0033, 3e-6],
    [2200, 5.235, 0.0018, 5e-5, 1e-6, 0.0029, 0.0060, 8e-6],
    [2300, 5.533, 0.0035, 0.0001, 5e-6, 0.0053, 0.0103, 2e-5],
    [2400, 5.821, 0.0063, 0.0002, 1e-5, 0.0091, 0.0168, 4e-5],
    [2500, 6.081, 0.0109, 0.0004, 5e-5, 0.0152, 0.0265, 9e-5],
    [2600, 6.331, 0.0182, 0.0009, 0.0001, 0.0244, 0.0403, 0.0002],
    [2700, 6.552, 0.0292, 0.0016, 0.0002, 0.0378, 0.0596, 0.0004],
    [2800, 6.764, 0.0453, 0.0026, 0.0005, 0.0567, 0.0856, 0.0006],
    [2900, 6.949, 0.0680, 0.0043, 0.0011, 0.0828, 0.1201, 0.0010],
    [3000, 7.127, 0.0995, 0.0070, 0.0023, 0.1179, 0.1648, 0.0017],
    [3100, 7.281, 0.1419, 0.0108, 0.0044, 0.1639, 0.2217, 0.0027],
    [3200, 7.428, 0.1981, 0.0163, 0.0081, 0.2233, 0.2926, 0.0042],
    [3300, 7.552, 0.2709, 0.0240, 0.0144, 0.2989, 0.3799, 0.0063],
    [3400, 7.670, 0.3638, 0.0345, 0.0248, 0.3931, 0.4858, 0.0092],
    [3500, 7.764, 0.4801, 0.0486, 0.0415, 0.5094, 0.6130, 0.0133],
    [3600, 7.854, 0.6236, 0.0671, 0.0674, 0.6503, 0.7631, 0.0188],
    [3700, 7.923, 0.799, 0.0911, 0.1065, 0.820, 0.940, 0.0261],
    [3800, 7.989, 1.009, 0.1216, 0.1642, 1.021, 1.144, 0.0357],
    [3900, 8.037, 1.260, 0.1599, 0.2478, 1.257, 1.380, 0.0480],
    [4000, 8.082, 1.556, 0.2073, 0.3663, 1.531, 1.648, 0.0637],
]
""" TABLE 2.06 from Hunt
            B                    C
T in K, H2, N2/CO, CO2, H2O, H2, N2/CO, CO2, H2O
"""
BCTable = [
    [1600, 16.4, 32.1, 45.7, -4.2, 20, 210, 1385, 220],
    [1700, 16.3, 32.3, 47.3, -2.5, 20, 200, 1305, 210],
    [1800, 16.2, 32.4, 48.7, -1.1, 20, 190, 1235, 195],
    [1900, 16.1, 32.6, 49.9, 0.2, 20, 180, 1170, 185],
    [2000, 16.0, 32.6, 50.9, 1.2, 15, 170, 1110, 175],
    [2100, 15.9, 32.7, 51.8, 2.2, 15, 160, 1055, 170],
    [2200, 15.8, 32.7, 52.6, 3.0, 15, 155, 1010, 160],
    [2300, 15.7, 32.8, 53.2, 3.7, 15, 150, 965, 155],
    [2400, 15.6, 32.8, 53.8, 4.4, 15, 140, 925, 145],
    [2500, 15.6, 32.8, 54.4, 5.0, 15, 135, 885, 140],
    [2600, 15.5, 32.7, 54.8, 5.5, 15, 130, 855, 135],
    [2700, 15.4, 32.7, 55.3, 6.0, 15, 125, 825, 130],
    [2800, 15.3, 32.7, 55.6, 6.4, 10, 120, 795, 125],
    [2900, 15.3, 32.6, 56.0, 6.8, 10, 120, 765, 120],
    [3000, 15.2, 32.6, 56.2, 7.1, 10, 115, 740, 120],
    [3100, 15.1, 32.6, 56.5, 7.5, 10, 110, 720, 115],
    [3200, 15.0, 32.5, 56.7, 7.7, 10, 105, 695, 110],
    [3300, 15.0, 32.4, 56.9, 8.0, 10, 105, 675, 105],
    [3400, 14.9, 32.4, 57.1, 8.3, 10, 100, 650, 105],
    [3500, 14.8, 32.3, 57.3, 8.5, 10, 95, 635, 100],
    [3600, 14.8, 32.3, 57.4, 8.7, 10, 95, 615, 100],
    [3700, 14.7, 32.2, 57.5, 8.9, 10, 90, 600, 95],
    [3800, 14.7, 32.2, 57.6, 9.1, 10, 90, 585, 95],
    [3900, 14.6, 32.1, 57.7, 9.3, 10, 85, 570, 90],
    [4000, 14.5, 32.0, 57.8, 9.4, 10, 85, 555, 90],
]


"""
(X): nbr. molecule per unit mass of propellant
{Y}: nbr. atom per unit mass of propellant

unless otherwise stated, mol/gram is assumed for
propellant work derived form Hunt.

Eq. 2.04
(CO) * (H2O) / [(CO2) * (H2)] = K0

Eq. 2.06
(N2)                    = 0.5 * {N}
(CO) + (CO2)            = {C}
(H2) + (H2O)            = 0.5 * {H}
(CO) + 2(CO2) + (H2O)   = {O}

From Eq.2.06 and Eq.2.04
Major                                   Minor
(N2)    = 0.5  {N}                      | - 0.5 (N) - 0.5 (NO)
(CO)    = {C} - (CO2)                   |
(H2O)   = {O} - {C} - (CO2)             | - (OH) - (NO) - 2 [ (O2) + .5 (O) ]
(H2)    = 0.5  {H} + {C} - {O} + (CO2)  | - 0.5 (H) + 2 (O2) + (o) + (NO) + 0.5 (OH)



n = {C} + 0.5 * {H} + 0.5 * {N}         | + (OH) + (H) + (NO) + (O2) + (O) + (N)
   CO/CO2     H2O/H2        N2

dissociaiton considered ,in descending order of significance:

2 H2O <-> 2 OH + H2             (OH)
H2 <-> 2 H                      (H)
2 H2O + N2 <-> 2 H2 + 2 NO      (NO)
2 H2O <-> 2 H2 + O2             (O2)
O2 <-> 2 O                      (O)
N2 <-> 2 N                      (N)

"""


def balance(T, Ci, Hi, Oi, Ni, V=1 / 0.1, tol=1e-5):  # in kelvin  # mol/g
    """
    Ci, Hi, Oi, Ni are in mol/g.
    Consequently,
    CO2j, OHj, Hj, NOj, O2j, Oj, Nj are also in mol/g

    T: temperature in Kelvin
    V: volume per unit mass of gas, in cc/g
    1/V is also commonly known as Delta, or load density, if the propellant solid
    covolume can be ignored.
    in literature, usually 1/V = 0.2 g/cc = 200 kg/m^3 is adopted for IB purposes.
    However, very often covolume values are given for a specific load condition.

        > a value of 0.089(0.084) g/cc is adopted for propellant referenced from
        ADA043778. When this value is adopted, prediciton is good to within 1% of
        tabulated.

        >It is possible that for modern software, even if the user sets a 1/V of
        0.2 g/cc, the program will apply its own averaging routine resulting in
        different result than is calculated here, given the apparent accuracy of
        all other parameter estimation performed here, and the closeness the com-
        position can be solved to published data.

    Covolume is is suitably invariant over typical internal ballistic tempearture,
    range, however it varies not-insignificantly by load density in test chambers,
    and thus must be chosen with reference to gas density achieved in actual gun
    system.

    n: number of gram-molecules per unit mass of gas, in mol/g
    """
    if T > 4000 or T < 1000:
        raise ValueError("T Not supported")
    for i in range(len(negDeltaTable) - 1):
        Tlow, negDeltaBlow, neghalfDeltaClow = negDeltaTable[i]
        Thigh, negDeltaBhigh, neghalfDeltaChigh = negDeltaTable[i + 1]
        if Tlow <= T <= Thigh:
            k = (T - Tlow) / (Thigh - Tlow)
            negDeltaB = negDeltaBlow * (1 - k) + negDeltaBhigh * k
            neghalfDeltaC = neghalfDeltaClow * (1 - k) + neghalfDeltaChigh * k
            break

    for i in range(len(equilibriumKT) - 1):
        Tlow, *Klow = equilibriumKT[i]
        Thigh, *Khigh = equilibriumKT[i + 1]
        if Tlow <= T <= Thigh:
            k = (1 / T - 1 / Tlow) / (1 / Thigh - 1 / Tlow)
            K = [
                10 ** (log10(Ki) * (1 - k) + log10(Kj) * k)
                for Ki, Kj in zip(Klow, Khigh)
            ]
            break

    for i in range(len(BCTable) - 1):
        Tlow, *BClow = BCTable[i]
        Thigh, *BChigh = BCTable[i + 1]
        if Tlow <= T <= Thigh:
            k = (T - Tlow) / (Thigh - Tlow)
            BC = [vi * (1 - k) + vj * k for vi, vj in zip(BClow, BChigh)]
            break

    BH2, BN2CO, BCO2, BH2O, CH2, CN2CO, CCO2, CH2O = BC

    R = 82.06  # in cc atm /(K mol)
    sqrtVdivRT = (V / (R * T)) ** 0.5

    OHj, Hj, NOj, O2j, Oj, Nj = 0, 0, 0, 0, 0, 0

    G = Ci
    H = Oi - Ci
    I = 0.5 * Hi + Ci - Oi

    oldCO2j = -1
    CO2j = 0.5 * (
        min(G, H) - I
    )  # initiate with a value that would not result in a negative fraction

    while True:
        n = Ci + 0.5 * Hi + 0.5 * Ni + OHj + Hj + NOj + O2j + Oj + Nj

        N2j = 0.5 * Ni - 0.5 * Nj - 0.5 * NOj
        G = Ci
        COj = G - CO2j
        H = Oi - Ci - OHj - NOj - 2 * O2j - Oj
        H2Oj = H - CO2j
        I = 0.5 * Hi + Ci - Oi - 0.5 * Hj + 2 * O2j + Oj + NOj + 0.5 * OHj
        H2j = I + CO2j

        K0 = K[0] * exp(n / V * negDeltaB + (n / V) ** 2 * neghalfDeltaC)
        oldCO2j = CO2j
        CO2j = quadratic((1 - K0), -(G + H + K0 * I), G * H)[1]
        CO2j = max(CO2j, 0)

        if abs(oldCO2j - CO2j) * 12.01 < tol:
            break

        OHj = H2Oj / (H2j) ** 0.5 * sqrtVdivRT * K[1] * exp(-20 * n / V)
        Hj = H2j**0.5 * sqrtVdivRT * K[5]
        NOj = H2Oj * N2j**0.5 / H2j * sqrtVdivRT * K[2] * exp(-20 * n / V)
        O2j = (H2Oj / H2j) ** 2 * sqrtVdivRT**2 * K[3]
        Oj = O2j**0.5 * sqrtVdivRT * K[4]
        Nj = N2j**0.5 * sqrtVdivRT * K[6]

    a = N2j + COj + CO2j + H2j + H2Oj + OHj + Hj + NOj + O2j + Oj + Nj

    speciesList = [
        ("N2", N2j * 28.016, N2j / a),
        ("CO", COj * 28.01, COj / a),
        ("CO2", CO2j * 44.01, CO2j / a),
        ("H2", H2j * 2.016, H2j / a),
        ("H2O", H2Oj * 18.016, H2Oj / a),
        ("OH", OHj * 17.008, OHj / a),
        ("H", Hj * 1.008, Hj / a),
        ("NO", NOj * 30.008, NOj / a),
        ("O2", O2j * 32.00, O2j / a),
        ("O", Oj * 16.00, Oj / a),
        ("N", Nj * 14.008, Nj / a),
    ]

    speciesList.sort(key=lambda x: x[1], reverse=True)

    n0 = Ci + 0.5 * Hi + 0.5 * Ni

    B = BCO2 * CO2j + BN2CO * (N2j + COj) + BH2O * H2Oj + BH2 * H2j
    C = CCO2 * CO2j + CN2CO * (N2j + COj) + CH2O * H2Oj + CH2 * H2j

    # b = V * (1 - n / n0 * (1 - B / V - n0 * C / V**2))

    b = (B * V**2 + n * C * V) / (V**2 + B * V + n * C)
    # rint("Calculated Conditions: {:>6.1f} K ".format(T))

    p = n * R * T / (V - b) / 9.869

    # p = n * R * T / V * (1 + B / V + n * C / V**2)

    return speciesList, b, p


if __name__ == "__main__":
    ingredients = Ingredient.readFile("data/hs.csv")
    """
    for T in (1600, 2000, 2400, 2800, 3100):
        print(balance(T=T, Ci=2232e-5, Hi=3010e-5, Ni=1046e-5, Oi=3469e-5))
    """
    # input()
    NC1260 = Ingredient.nitrocellulose(0.1260)
    RDX = Ingredient.find("RDX")
    NG = Ingredient.find("NG")
    C = Ingredient.find("Graphite")

    EtNENA = Ingredient.fromElement(
        name="2-(ethylnitroamino)ethyl nitrate",
        alt="Ethyl NENA",
        C=4,
        H=9,
        N=3,
        O=5,
        HoC=644.43,
        u="kcal/mol",
    )

    MeNENA = Ingredient.fromElement(
        name="2-(methylnitroamino)ethyl nitrate",
        alt="Methyl NENA",
        C=3,
        H=7,
        N=3,
        O=5,
        HoC=485.46,
        u="kcal/mol",
    )
    UREA = Ingredient.fromElement(
        name="Akardite II",
        alt="Urea",
        C=14,
        H=14,
        N=2,
        O=1,
        HoC=1784.2,
        u="kcal/mol",
    )

    PRD20 = Mixture(
        name="ATK RPD20",
        compoDict={
            NC1260: 0.4190,
            RDX: 0.2571,
            MeNENA: 0.14,
            EtNENA: 0.10,
            NG: 0.0769,
            UREA: 0.0070,
        },
        Delta=0.1,
    )
    PRD20.prettyPrint()

    PRD21 = Mixture(
        name="ATK RPD(S)21",
        compoDict={
            NC1260: 0.3648,
            RDX: 0.3033,
            MeNENA: 0.1344,
            EtNENA: 0.0957,
            NG: 0.0946,
            UREA: 0.0072,
        },
        Delta=0.1,
    )
    PRD21.prettyPrint()

    PRD22 = Mixture(
        name="ATK RPD(S)22",
        compoDict={
            NC1260: 0.3111,
            RDX: 0.3408,
            MeNENA: 0.1257,
            EtNENA: 0.0894,
            NG: 0.1258,
            UREA: 0.0072,
        },
        Delta=0.1,
    )
    PRD22.prettyPrint()

    # Ingredient.check()

    NC1320 = Ingredient.nitrocellulose(0.1320)
    # NC1320.prettyPrint()
    DEGDN = Ingredient.find("DEGDN")

    JA2 = Mixture(
        name="JA2",
        compoDict={
            NC1320: 0.5950,
            DEGDN: 0.2480,
            NG: 0.1490,
            UREA: 0.0070,
            C: 0.0010,
        },
        Delta=0.1,
    )
    JA2.prettyPrint()

    NC1315 = Ingredient.nitrocellulose(0.1315)
    DNT = Ingredient.find("DNT")
    DBP = Ingredient.find("Dibutylphathalate")
    DPA = Ingredient.find("Diphenlyamine")
    M1 = Mixture(
        name="M1",
        compoDict={NC1315: 0.8500, DNT: 0.1000, DBP: 0.0500, DPA: 0.01},
        Delta=0.089,
    )
    M1.prettyPrint()
