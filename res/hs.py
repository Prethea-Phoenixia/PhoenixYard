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

Results are within few percent of accurate thermalchemical calculation
in the range of 2000K to 4000K. Should not be applied to cases >4000K
where dissociation of of free radicals like -H, -OH and -Cl would happen.
Should not be applied to propellants generating substantial amount of 
consdensed exhaust.


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
    def find(cls, name):
        """
        return an Ingredient object that is closest to the supplied common name,
        or alternative (abbreviated) name.
        """
        closeIngrs = difflib.get_close_matches(
            name,
            list(cls.allDict.keys()) + list(cls.altDict.keys()),
            n=5,
            cutoff=0.5,
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
        Q = 0  # heat of explosion
        Cv = 0  # heat capacity at constant volume4
        E = 0  # heat of formation
        invM = 0  # gas volume, mol/g

        for ingr, fraction in self.compoDict.items():
            Q += fraction * ingr.Qi
            Cv += fraction * ingr.Cvi
            E += fraction * ingr.Ei
            invM += fraction * ingr.invMi

        Cp = Cv * gamma
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

        self.Tv = Tv
        self.gamma = gamma
        self.f = f


if __name__ == "__main__":
    ingredients = Ingredient.readFile("hs.csv")
    NC = Ingredient.find("Nitrocellulose 12.6%")
    NG = Ingredient.find("NG")
    EC = Ingredient.find("Centraltie")

    testMix = Mixture("test", {NC: 0.5, NG: 0.49, EC: 0.01})

    print(testMix.Tv)
    print(testMix.gamma)
    print(testMix.f)
