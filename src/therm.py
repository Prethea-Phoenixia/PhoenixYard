import difflib

from corner import balance
from num import secant
from periodic import molarMasses


class Ingredient:
    allIngr = {}
    lastLine = 0
    lineIngr = {}

    def __init__(
        self,
        name,
        elements,
        Hf,
        rho,
        rho_u="lb/cu.in",
        Hf_u="cal/g",
        flag="",
        lineNo=None,
    ):
        if lineNo is not None:
            if lineNo > Ingredient.lastLine:
                Ingredient.lastLine = lineNo
            else:
                if lineNo not in Ingredient.lineIngr:
                    pass
                else:
                    raise ValueError("Line Number Collision")
        else:
            lineNo = Ingredient.lastLine + 1
            Ingredient.lastLine += 1

        self.lineNo = lineNo
        self.name = name
        self.flag = flag
        self.elements = elements

        if rho_u == "lb/cu.in":
            self.rho = rho * 27.680  # to g/cc
        elif rho_u == "g/cc":
            self.rho = rho
        else:
            raise ValueError("Unknown Unit For Density")

        A = 0
        for element, num in self.elements.items():
            A += molarMasses[element] * num

        if "C" in self.elements:
            self.Ci = self.elements["C"] / A
        else:
            self.Ci = 0
        if "H" in self.elements:
            self.Hi = self.elements["H"] / A
        else:
            self.Hi = 0
        if "O" in self.elements:
            self.Oi = self.elements["O"] / A
        else:
            self.Oi = 0
        if "N" in self.elements:
            self.Ni = self.elements["N"] / A
        else:
            self.Ni = 0

        self.A = A  # g/mol

        if Hf_u == "cal/g":
            self.Hf = Hf
        elif Hf_u == "cal/mol":
            self.Hf = Hf / A
        elif Hf_u == "J/g":
            self.Hf = Hf / 4.184
        elif Hf_u == "J/mol":
            self.Hf = Hf / (4.184 * A)
        else:
            raise ValueError("Unknown Enthalpy Unit")

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
                        name=name,
                        elements=elements,
                        Hf=Hf,
                        rho=rho,
                        flag=flag,
                        lineNo=lineNo,
                    )

                    fileIngr.append(newIngr)

        for ingr in fileIngr:
            cls.allIngr.update({ingr.name: ingr})
            cls.lineIngr.update({ingr.lineNo: ingr})

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

    @classmethod
    def getLine(cls, lineNo):
        if lineNo in cls.lineIngr:
            print(
                "Returning line {:} : {:}".format(
                    lineNo, cls.lineIngr[lineNo].name
                )
            )
            return cls.lineIngr[lineNo]
        else:
            print("No such line as {:}\n".format(lineNo))
            return None


class Mixture:
    def __init__(self, name, compoDict, Delta=0.2, tol=1e-5):
        self.name = name
        self.Delta = Delta  # load density in g/cc

        # Normalize the given composition such that the fractions sums to 1

        total = 0
        for ingr, fraction in compoDict.items():
            total += fraction

        self.compoDict = {
            ingr: fraction / total for ingr, fraction in compoDict.items()
        }

        # tally the releavnt factors according to their mass fraction

        invRho = 0
        Ci, Hi, Ni, Oi = 0, 0, 0, 0
        Hf = 0

        for ingr, fraction in self.compoDict.items():
            if ingr.rho == 0:
                raise ValueError(
                    "{:} is not provided with density data".format(ingr.name)
                )
            invRho += fraction / ingr.rho
            Ci += fraction * ingr.Ci  # mol/g
            Hi += fraction * ingr.Hi
            Oi += fraction * ingr.Oi
            Ni += fraction * ingr.Ni
            Hf += fraction * ingr.Hf

        self.rho = 1 / invRho

        def f(T):
            DeltaE, _, _, _, _, _, _ = balance(
                T, Ci, Hi, Oi, Ni, V=1 / Delta, tol=tol
            )

            return DeltaE - Hf

        Tv, _ = secant(f, 2500, 3500, x_min=1600, x_max=4000, tol=tol)

        _, _, n, self.speciesList, self.b, self.p, self.f = balance(
            Tv, Ci, Hi, Oi, Ni, V=1 / Delta, tol=tol
        )
        # see Hunt ss 2.13
        _, E1, _, _, _, _, _ = balance(Tv, Ci, Hi, Oi, Ni, V=1 / 0.2, tol=tol)
        _, E2, _, _, _, _, _ = balance(
            0.6 * Tv, Ci, Hi, Oi, Ni, V=1 / 0.1, tol=tol
        )
        sigma_v = (E1 - E2) / (0.4 * Tv)
        self.gamma = (n * 1.987 / sigma_v) + 1
        self.C = Ci * molarMasses["C"]  # weight fraction
        self.H = Hi * molarMasses["H"]
        self.O = Oi * molarMasses["O"]
        self.N = Ni * molarMasses["N"]

        self.Tv = Tv
        self.Hf = Hf

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
        print("Calculated Properties:---------------------------")
        print("Density            : {:>6.4g} g/cc".format(self.rho))
        print("Heat of Formation  : {:>6.3g} cal/g".format(self.Hf))
        print(
            "Flame Temperature  : {:>6.4g} K (Isochoric Adiabatic)".format(
                self.Tv
            )
        )
        print(" @ Product  %mass  mol/g")
        print(
            *[
                "{:>2} : {:^6} {:<6.1%} {:<6.4f}".format(i, name, mass, num)
                for i, (name, mass, num) in enumerate(self.speciesList)
            ],
            sep="\n"
        )
        print("Impetus / Force    : {:>6.4g} J/g".format(self.f))
        print("Covolume           : {:>6.4g} cc/g".format(self.b))
        # print(" @Temperature      : {:>6.0f} K".format(self.Tv))
        print(" @ Load Density    : {:>6.3g} g/cc".format(self.Delta))
        print(" @ Pressure        : {:>6.4g} MPa".format(self.p))
        print("avg Adb. index     : {:>6.4g}".format(self.gamma))
        print("")


if __name__ == "__main__":
    Ingredient.readFile("data/PEPCODED.DAF")
    NC1260 = Ingredient.getLine(683)
    RDX = Ingredient.getLine(847)
    # EC = Ingredient.getLine(397)

    EC = Ingredient(
        name="Ethyl Centralite",
        elements={"C": 17, "H": 20, "O": 1, "N": 2},
        rho=1.140,
        rho_u="g/cc",
        Hf=-391.5,
        Hf_u="J/g",
    )

    ATEC = Ingredient(
        name="Acetyl triethyl citrate",
        elements={"C": 14, "H": 22, "O": 8},
        rho=1.136,
        rho_u="g/cc",
        # Hf=-1257,
        Hf=-5459.6,
        Hf_u="J/g",
    )

    CAB = Ingredient(
        "Cellulose Acetate Butyrate",
        elements={"C": 15, "H": 22, "O": 8},
        Hf=-4933.76,
        Hf_u="J/g",
        rho=1.220,
        rho_u="g/cc",
    )

    BDNPA = Ingredient.getLine(189)
    BDNPF = Ingredient.getLine(190)
    """
    XM39 = Mixture(
        "XM39",
        compoDict={RDX: 76, CAB: 12, NC1260: 4, ATEC: 7.6, EC: 0.4},
    )

    XM39.prettyPrint()
    """
    M43 = Mixture(
        name="M43",
        compoDict={
            RDX: 76,
            CAB: 12,
            NC1260: 4,
            BDNPA: 7.6,
            # BDNPF: 7.6,
            EC: 0.4,
        },
        Delta=0.2,
    )
    M43.prettyPrint()
