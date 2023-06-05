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

    def __init__(self, name, alt, Qi, Cvi, Ei, invMi):
        self.name = name
        self.alt = alt
        self.Qi = Qi
        self.Cvi = Cvi
        self.Ei = Ei
        self.invMi = invMi

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
                (name, alt, Qi, Cvi, Ei, invMi) = ingr

                newIngr = cls(
                    name, alt, float(Qi), float(Cvi), float(Ei), float(invMi)
                )

                ingredients.append(newIngr)

        ingrDict = {i.name: i for i in ingredients}
        cls.allDict.update(ingrDict)
        cls.altDict.update({i.alt: i for i in ingredients if i.alt != ""})
        return ingrDict

    @classmethod
    def fromElement(cls, name, C, H, O, N, HoC, alt="", u="kJ/mol", keep=True):
        """
        Given the molecular formula of a chemical, estimate factors necessary for
        use in the Hirschfelder-Sherman calculation, and add the newly created
        ingredient into the class.
        """

        # accurate molecular mass here to account for natural abundance of isotopes
        A = 12.01 * C + 1.008 * H + 16.00 * O + 14.008 * N  # g/mol

        """
        convert element nbr. (mol/mol) into nbr. mol of element per unit mass (mol/g)
        """
        C /= A
        H /= A
        O /= A
        N /= A

        if u == "kJ/mol":
            HoC /= 4.184  # to kcal/mol
            HoC /= A  # to kcal/g
            HoC *= 1000  # to cal/g
        elif u == "kcal/mol":
            HoC /= A
            HoC *= 1000
        elif u == "kcal/kg" or u == "cal/g":
            pass  # kcal/kg = cal/g
        else:
            raise ValueError("Unknown unit ", u)

        invM = C + 0.5 * H + 0.5 * N
        # isochoric heat capacity from 2000-3000K
        Cv = 1.620 * C + 3.265 * H + 5.193 * O + 3.384 * N
        # HoC: heat of combustion, +: energy is released, -: energy is consumed
        # this is the opposite of the usual convention of enthalpy of combustion.
        Q = HoC - 67421 * (2 * C + 0.5 * H - O)
        # heat of explosion, +: heat is released, -: heat is consumed, unknown unit.
        E = HoC - 132771 * C - 40026 * H + 51819 * O - 6724 * N
        # heat required to bring combustion product to 2500k, +: product hotter than 2.5kK
        # -: product cooler than 2.5kK

        newIngr = cls(name, alt, Q, Cv, E, invM)
        if keep:
            cls.allDict.update({newIngr.name: newIngr})
            if newIngr.alt != "":
                cls.altDict.update({newIngr.alt: newIngr})

        return newIngr

    @classmethod
    def nitrocellulose(cls, Nfraction, keep=True):
        """
        Initialize nitrocellulose ingredient with arbitrary nitrogen mass fraction,
        by calculating the corresponding elemental makeup and heat of combustion.

        If the an entry at the same nitration level is found, the keep parameter is
        ignored and the class dictionary is not updated.
        """

        if Nfraction < 0.1150 or Nfraction > 0.1360:
            raise ValueError(
                "Nitrogen fraction between 11.50% and 13.50% are supported,"
                " data is unavailable outside this range."
            )
        f = Nfraction
        x = 162.14 * f / (14.008 - 45 * f)
        C = 6
        H = 10 - x
        O = 5 + 2 * x
        N = x

        name = "Nitrocellulose {:.2%} N".format(f)
        alt = "NC{:d}".format(int(f * 10000))

        HoC = -4176.70 + 14126 * f  # cal/g @ 30 deg
        newIngr = cls.fromElement(
            name=name, alt=alt, C=C, H=H, O=O, N=N, HoC=-HoC, u="cal/g"
        )
        if keep and name not in cls.allDict:
            cls.allDict.update({newIngr.name: newIngr})
            cls.altDict.update({newIngr.alt: newIngr})

        """
        Custom, bespoke factors, extrapolated from the available NC entries at 12.20%,
        12.60% and 13.15% nitrogen level.
        """
        newIngr.invMi += 0.00002
        newIngr.Ei += (Nfraction - 0.1260) * 1075 + 13
        newIngr.Qi += 13

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
        for ingr in cls.allDict.values():
            ingr.prettyPrint()
            print()

    def prettyPrint(self):
        print(str(self))

        print("Hirschfelder-Sherman factors:--------------------")
        print("{:>5}  {:>6}  {:>7}  {:>7}".format("Qi", "Cvi", "Ei", "1/Mi"))
        print(
            "{:>5.0f}  {:>6.4f}  {:>7.1f}  {:>7.5f}".format(
                self.Qi, self.Cvi, self.Ei, self.invMi
            )
        )

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
        Q = 0  # heat of explosion, presumably in cal/g
        Cv = 0  # heat capacity at constant volume, presumably in cal/g K
        E = 0  # heat of formation, also presumably in cal/g
        invM = 0  # specific gas volume, mol/g

        for ingr, fraction in self.compoDict.items():
            Q += fraction * ingr.Qi
            Cv += fraction * ingr.Cvi
            E += fraction * ingr.Ei
            invM += fraction * ingr.invMi

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

    def prettyPrint(self):
        print("Mixture: {:}".format(self.name))
        print("Specified Composition:---------------------------")
        for ingr, fraction in self.compoDict.items():
            print("--{:-<30}, {:>6.1%}".format(str(ingr), fraction))

        print("")
        print("Calcualted Properties:---------------------------")
        print("Isochoric Adb. Temp: {:>4.1f}K".format(self.Tv))
        print("Adiabatic Index    : {:>4.3f}".format(self.gamma))
        print("Specific Force     : {:>4.3f} MJ/kg".format(self.f / 1e6))


if __name__ == "__main__":
    ingredients = Ingredient.readFile("data/hs.csv")
    DNT = Ingredient.fromElement(
        "2,4-Dinitrotoluene",
        alt="DNT",
        C=7,
        H=6,
        N=2,
        O=4,
        HoC=3560,
        u="kJ/mol",
    )
    """
    NG = Ingredient.find("NG")
    NG.prettyPrint()
    NGTEST = Ingredient.fromElement(
        "NitroglycerinTest",
        alt="NGT",
        C=3,
        H=5,
        N=3,
        O=9,
        HoC=365,  # from NIST
        u="kcal/mol",
        keep=False,
    )
    NGTEST.prettyPrint()



    for f in (0.1220, 0.1260, 0.1315):
        alt = "NC{:d}".format(int(f * 10000))
        NC = Ingredient.find(alt)
        NC.prettyPrint()

        NCTEST = Ingredient.nitrocellulose(f, keep=False)
        NCTEST.prettyPrint()
        print("")
    """
    NC1315 = Ingredient.find("Nitrocellulose 13.15% N")
    DBP = Ingredient.find("Dibutyl Phthalate")
    DPA = Ingredient.find("DPA")

    testMix = Mixture(
        "M1",
        {NC1315: 0.85, DNT: 0.10, DBP: 0.05, DPA: 0.01},
    )

    testMix.prettyPrint()
