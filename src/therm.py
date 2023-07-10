import difflib

from hunt import balance


# the E stand for electron and D for deuterium
# fmt: off
molarMasses = {
    "H": 1.00794, "HE": 4.002602, "LI": 6.941, "BE": 9.012182,
    "B": 10.811, "C": 12.0107, "N": 14.00674, "O": 15.9994,
    "F": 18.9984032, "NE": 20.11797, "NA": 22.989770, "MG": 24.305,
    "AL": 26.981538, "SI": 28.0855, "P": 30.973761, "S": 32.066,
    "CL": 35.4527, "AR": 39.948, "K": 39.0983, "CA": 40.078,
    "SC": 44.95591, "TI": 47.88, "V": 50.9415, "CR": 51.996,
    "MN": 54.938, "FE": 55.847, "CO": 58.9332, "NI": 58.6934,
    "CU": 63.546, "ZN": 65.39, "GA": 69.723, "GE": 72.61,
    "AS": 74.9216, "SE": 78.96, "BR": 79.904, "KR": 83.80,
    "RB": 85.4678, "SR": 87.62, "Y": 88.9059, "ZR": 91.224,
    "NB": 92.9064, "MO": 95.94, "TC": 98.0, "RU": 101.07,
    "RH": 102.9055, "PD": 106.42, "AG": 107.868, "CD": 112.41,
    "IN": 114.82, "SN": 118.71, "SB": 121.757, "TE": 127.60,
    "I": 126.9045, "XE": 131.29, "CS": 132.9054, "BA": 137.33,
    "LA": 138.9055, "CE": 140.12, "PR": 140.9077, "ND": 144.24,
    "PM": 145.0, "SM": 150.36, "EU": 151.965, "GD": 157.25,
    "TB": 158.9253, "DY": 162.50, "HO": 164.9303, "ER": 167.26,
    "TM": 168.9342, "YB": 173.04, "LU": 174.967, "HF": 178.49,
    "TA": 180.9479, "W": 183.85, "RE": 186.207, "OS": 190.2,
    "IR": 192.22, "PT": 195.08, "AU": 196.9665, "HG": 200.59,
    "TL": 204.383, "PB": 207.2, "BI": 208.9804, "PO": 209.0,
    "AT": 210.0, "RN": 222.0, "FR": 223.0, "RA": 226.0254,
    "AC": 227.0, "TH": 232.0381, "PA": 231.0359, "U": 238.029,
    "NP": 237.0482, "PU": 244.0, "FM": 257.0, "E": 0, "D": 2,
}
# fmt: on
class Ingredient:
    allIngr = {}

    def __init__(self, name, elements, Hf, rho, flag=""):
        self.name = name
        self.flag = flag
        self.elements = elements
        self.Hf = Hf
        self.rho = rho
        A = 0
        for element, num in self.elements.items():
            print(element, num)

            A += molarMasses[element] * num
        self.A = A

    @classmethod
    def readFile(cls, fileName):
        # read data from PEP database
        with open("data/PEPCODED.DAF", "r", encoding="ascii") as file:
            fileIngr = []
            for line in file:
                # print(line, end="")
                flag = line[0:2].strip()
                lineNo = int(line[2:9])

                if flag == "*":  # this line is a comment
                    pass
                elif flag == "+":
                    # multiline entry of ingredient name
                    if len(fileIngr) > 0:
                        newIngr = fileIngr[-1]
                        name = line[9:80].strip()
                        newIngr.name = newIngr.name + name
                else:
                    name = line[9:39].strip()
                    elements = dict()
                    for i in range(6):
                        nbr, element = (
                            int(line[39 + i * 5 : 42 + i * 5].strip()),
                            line[42 + i * 5 : 44 + i * 5].strip(),
                        )
                        if element != "":
                            elements.update({element: nbr})

                    Hf = int(line[69:74])
                    rho = float(line[74:80])

                    newIngr = Ingredient(
                        name=name, elements=elements, Hf=Hf, rho=rho, flag=flag
                    )

                    fileIngr.append(newIngr)

        for ingr in fileIngr:
            cls.allIngr.update({ingr.name: ingr})

    @classmethod
    def find(cls, name):
        closeIngrs = difflib.get_close_matches(
            name,
            list(cls.allIngr.keys()),
            n=3,
            cutoff=0.6,
        )
        n = 0
        thisIngr = None
        ingrs = []

        if len(closeIngrs) > 0:
            print("Found candidate:")
            for iname in closeIngrs:
                if iname in cls.allIngr:
                    ingr = cls.allIngr[iname]
                if ingr not in ingrs:
                    ingrs.append(ingr)
                if n == 0:
                    thisIngr = ingr
                n += 1

            for ingr in ingrs:
                print("-" + ingr.name)

            print("returning " + thisIngr.name + "\n")
            return thisIngr
        else:
            print('Unknown ingredient description "{:}"'.format(name) + "\n")
            return None


class Mixture:
    def __init__(self, name, compoDict, Delta=0.2, tol=1e05):
        self.name = name
        self.Delta = Delta  # load density in g/cc

        # Normalize the given composition such that the fractions sums to 1

        total = 0
        for ingr, fraction in compoDict.items():
            total += fraction

        self.compoDict = {
            ingr: fraction / total for ingr, fraction in compoDict.items()
        }

        """
        tally the releavnt factors according to their mass fraction
        """
        invRho = 0
        Ci, Hi, Ni, Oi = 0, 0, 0, 0

        for ingr, fraction in self.compoDict.items():
            invRho += fraction / ingr.rho

            if "C" in ingr.elements:
                Ci += fraction * ingr.elements["C"]
            if "H" in ingr.elements:
                Hi += fraction * ingr.elements["H"]
            if "O" in ingr.elements:
                Oi += fraction * ingr.elements["O"]
            if "N" in ingr.elements:
                Ni += fraction * ingr.elements["N"]

        self.rho = 1 / invRho * 27680  # convert to kg/m^3

        def f(T):
            DeltaE, self.speciesList, self.b, self.p = balance(
                T, Ci, Hi, Oi, Ni, V=1 / Delta, tol=tol
            )

        self.C = Ci * 12.01
        self.H = Hi * 1.008
        self.O = Oi * 16.0
        self.N = Ni * 14.01

    def prettyPrint(self):
        print("Mixture: {:}".format(self.name))
        print("Specified Composition:---------------------------")
        for ingr, fraction in self.compoDict.items():
            print("--{:-<30}, {:>6.2%}".format(ingr.name, fraction))

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
        print("Heat of Explosion  : {:>6.1f} cal/g".format(self.Q))
        print("")

        print("Species--%wt.-----mol/g-----------------------------")
        print(
            *("{:^5}  {:>6.2%}  {:>8.6f}".format(*v) for v in self.speciesList),
            sep="\n"
        )
        print("")
        print("Further thermalchemical Calculations:------------")

        print("Covolume           : {:>6.4g} cc/g".format(self.b))
        print(" @Temperature      : {:>6.1f} K".format(self.Tv))
        print(" @Load Density     : {:>6.3g} g/cc".format(self.Delta))
        print(" @Pressure         : {:>6.4g} MPa".format(self.p))
        print("")


if __name__ == "__main__":
    Ingredient.readFile("data/PEPCODED.DAF")
    # print(*Ingredient.allIngr.keys(), sep="\n")
