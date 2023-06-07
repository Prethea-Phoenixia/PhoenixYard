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
from num import *


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
    def __init__(self, name, compoDict):
        self.name = name
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
        # self.b = 1e-3 * (1.18 + 6.9 * Ci - 11.5 * Oi)
        """
        self.b = 1e-3 * (
            1.3372800285521016 + 15.00623286 * Ci - 20.81820076 * Oi
        )
        self.b = 1e-3 * (
            -15.125770768899
            + 202.80374114 * Ci
            + 27.53266578 * Hi
            + 222.42131159 * Ni
            + 242.02396858 * Oi
        )
        """
        #  A = 12.01 * C + 1.008 * H + 14.008 * N + 16.00 * O  # g/mol
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
        print("Isochoric Adb. Temp: {:>4.1f}K".format(self.Tv))
        print("Adiabatic Index    : {:>4.3f}".format(self.gamma))
        print("Specific Force     : {:>4.4f} MJ/kg".format(self.f / 1e6))
        print("")


# cc/(gm.mol)
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

EquilibriumKT = [[1000, 0.7185, 3e-11, 1e-14, 8e-21, 2e-10, 3e-9, 5e-16]]

if __name__ == "__main__":
    ingredients = Ingredient.readFile("data/hs.csv")

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
    )
    JA2.prettyPrint()

    NC1315 = Ingredient.nitrocellulose(0.1315)
    DNT = Ingredient.find("DNT")
    DBP = Ingredient.find("Dibutylphathalate")
    DPA = Ingredient.find("Diphenlyamine")
    M1 = Mixture(
        name="M1",
        compoDict={NC1315: 0.8500, DNT: 0.1000, DBP: 0.0500, DPA: 0.01},
    )
    M1.prettyPrint()
