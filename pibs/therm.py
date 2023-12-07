import difflib

from corner import balance
from num import secant
from periodic import molarMasses
from math import exp


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
        elif Hf_u == "kJ/mol":
            self.Hf = Hf * 1000 / (4.184 * A)
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

    @classmethod
    def nitrocellulose(cls, nitration):
        y = nitration * 100
        x = 162.14 * y / (1400.8 - 45 * y)
        elements = {"C": 6, "H": 10 - x, "O": 5 + 2 * x, "N": x}
        # see hunt SS 2.02

        return cls(
            name="Nitrocellulose ({:.2f}% N)".format(y),
            elements=elements,
            rho=0.0560,  # copied directly from database.
            Hf=-1417.029
            + 6318.3 * nitration,  # fit from Tab.3 in R.S.Jessup & E.J.Prosen
        )


class Mixture:
    def __init__(
        self, name, compoDict, Delta=0.2, tol_z=1e-3, tol_b=1e-9, its=100
    ):
        self.name = name
        self.Delta = Delta  # load density in g/cc
        self.tol_z = tol_z  # tolerance for zeta
        self.its = its
        self.tol_b = tol_b

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
        self.Hf = Hf

        def f(T):
            zeta, _, _, _, _, _, _ = balance(
                self.Hf, T, Ci, Hi, Oi, Ni, V=1 / Delta, its=its, tol=tol_b
            )

            return zeta

        # zeta on the order of 0.5 per degree
        Tv = 0.5 * sum(
            secant(f, 2500, 3500, x_min=1600, x_max=4000, x_tol=tol_z)
        )

        _, self.speciesList, n, E, self.b, self.p, self.f = balance(
            self.Hf, Tv, Ci, Hi, Oi, Ni, V=1 / Delta, its=its, tol=tol_b
        )
        # see Corner ss 3.4
        _, _, _, E1, _, _, _ = balance(
            self.Hf, Tv, Ci, Hi, Oi, Ni, V=1 / Delta, its=its, tol=tol_b
        )
        _, _, _, E2, _, _, _ = balance(
            self.Hf, 0.7 * Tv, Ci, Hi, Oi, Ni, V=1 / Delta, its=its, tol=tol_b
        )
        C_v = (E1 - E2) / (0.3 * Tv)
        # gas constant: 1.987 cal/(mol K)
        self.gamma = (n * 1.987 / C_v) + 1
        self.n = n

        self.Ci, self.Hi, self.Oi, self.Ni = Ci, Hi, Oi, Ni

        self.Tv = Tv
        self.Hf = Hf

    def balanceAt(self, T, verbose=True, its=100, tol=1e-9):
        Delta, speciesList, n, E, b, p, f = balance(
            self.Hf,
            T,
            self.Ci,
            self.Hi,
            self.Oi,
            self.Ni,
            V=1 / self.Delta,
            its=its,
            tol=tol,
        )
        print(Delta)
        if verbose:
            print("Mixture: {:} At: {:}K".format(self.name, T))

            print(" @ Product  %mass  mol/g")
            print(
                *[
                    "{:>2} : {:^6} {:<6.1%} {:<6.4f}".format(i, name, mass, num)
                    for i, (name, mass, num) in enumerate(speciesList)
                ],
                sep="\n",
            )
            print(
                "Average.Mol.Weight : {:>6.4g} g/mol Δ={:>6.1%}".format(
                    1 / n, (1 / n - 1 / self.n) / (1 / self.n)
                )
            )
            print(
                "Covolume           : {:>6.4g} cc/g  Δ={:>6.1%}".format(
                    b, (b - self.b) / self.b
                )
            )

            C_v = E / (T - 300)
            gamma = (n * 1.987 / C_v) + 1
            print("Adiabatic index (T): {:>6.4g}".format(gamma))

        return speciesList

    def prettyPrint(self):
        C, H, O, N = (
            self.Ci * molarMasses["C"],
            self.Hi * molarMasses["H"],
            self.Oi * molarMasses["O"],
            self.Ni * molarMasses["N"],
        )
        print("Mixture: {:}".format(self.name))
        print("Specified Composition:---------------------------")
        for ingr, fraction in self.compoDict.items():
            print("--{:-<30}, {:>6.2%}".format(ingr.name, fraction))

        print("")
        print("Elemental Fractions:-----------------------------")
        print("C {:.2%} H {:.2%} N {:.2%} O {:.2%}".format(C, H, N, O))
        print("")
        print("Calculated Properties:---------------------------")
        print("Density            : {:>8.4g} g/cc".format(self.rho))
        print("Heat of Formation  : {:>8.0f} cal/g".format(self.Hf))
        print(
            "Flame Temperature  : {:>8.6g} K (Isochoric Adiabatic)".format(
                self.Tv
            )
        )
        print(" @ Product  %mass  mol/g")
        print(
            *[
                "{:>2} : {:^6} {:<6.1%} {:<6.4f}".format(i, name, mass, num)
                for i, (name, mass, num) in enumerate(self.speciesList)
            ],
            sep="\n",
        )
        print("Impetus / Force    : {:>8.5g} J/g".format(self.f))
        print("Covolume           : {:>8.4g} cc/g".format(self.b))
        # print(" @Temperature      : {:>6.0f} K".format(self.Tv))
        print(" @ Load Density    : {:>8.6g} g/cc".format(self.Delta))
        print(" @ Pressure        : {:>8.1f} MPa".format(self.p))
        print("avg Adb. index     : {:>8.5g}".format(self.gamma))
        print("")
        print("Estimated Burnrates:-----------------------------")
        b = -0.8340 + 8.3956e-4 * self.Tv
        print("after Ludwig Stiefel (ADA086093):")
        print(f"r <mm/s> = {b:.3f} (P <MPa>) ^ 0.8053")
        print(f"r <m/s> = {b * 1e-3 / 1e6**0.8053:.4g} (P <Pa>) ^ 0.8053")
        print("")


def estU1(nitration, ng, cen, dbp, dnt, vas):
    """
    Estimate the burn velocity assuming a linear realtionship
    between pressure and burn rate, i.e. u = u1 * P

    nitration:
        N, nitration % of nitrocellulose
        硝化棉中氮的百分含量
    ng:
        H, percentage of nitroglycerin
        硝化甘油的百分含量

    cent:
        E, percentage of centralite
        中定剂的百分含量
    dbp:
        F, percentage of dibutyl phthalate
        苯二甲酸二丁酯的百分含量
    dnt:
        G, percentage of dinitrotoluene
        二硝基甲苯的百分含量
    vaseline:
        I, 凡士林的百分含量

    The result is converted into m/s/Pa
    """
    N, H, E, F, G, I = nitration, ng, cen, dbp, dnt, vas

    return (
        10
        ** (
            0.7838
            - 9
            + 0.1366 * (N - 11.8)
            + 0.008652 * H
            - 0.02620 * E
            - 0.02235 * F
            - 0.007355 * G
            - 0.03447 * I
        )  # m^3/s-kg
        / 9.8
    )


def estU8053(Tv):
    return (-0.8340 + 8.3956e-4 * Tv) * 1e-3 / 1e6**0.8053


if __name__ == "__main__":
    Ingredient.readFile("data/PEPCODED.DAF")

    NC1260 = Ingredient.getLine(683)
    RDX = Ingredient.getLine(847)

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

    XM39 = Mixture(
        "XM39",
        compoDict={RDX: 76, CAB: 12, NC1260: 4, ATEC: 7.6, EC: 0.4},
    )

    XM39.prettyPrint()

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

    NG = Ingredient.getLine(693)

    MeNENA = Ingredient(
        "Methyl-NENA",
        elements={"C": 3, "H": 7, "N": 3, "O": 5},
        Hf=-106.50,  # Burcat, 2015
        Hf_u="kJ/mol",
        rho=1.53,  # a.la ADA377866
        rho_u="g/cc",
    )

    EtNENA = Ingredient(
        "Ethyl-NENA",
        elements={"C": 4, "H": 9, "N": 3, "O": 5},
        Hf=-133.90,  # Burcat, 2015
        Hf_u="kJ/mol",
        rho=1.32,  # a.la ADA377866
        rho_u="g/cc",
    )

    ATKPRD20 = Mixture(
        name="ATK PRD20",
        compoDict={
            NC1260: 41.90,
            RDX: 25.71,
            MeNENA: 14.00,
            EtNENA: 10.00,
            NG: 7.69,
        },
        Delta=0.2,
    )

    ATKPRD20.prettyPrint()

    ATKPRDS21 = Mixture(
        name="ATK PRD(S)21",
        compoDict={
            NC1260: 36.48,
            RDX: 30.33,
            MeNENA: 13.44,
            EtNENA: 9.57,
            NG: 9.46,
        },
        Delta=0.2,
    )

    ATKPRDS21.prettyPrint()

    ATKPRDS22 = Mixture(
        name="ATK PRD(S)22",
        compoDict={
            NC1260: 31.11,
            RDX: 34.08,
            MeNENA: 12.57,
            EtNENA: 8.94,
            NG: 12.58,
        },
        Delta=0.2,
    )

    ATKPRDS22.prettyPrint()

    """

    import matplotlib.pyplot as plt
    from labellines import labelLines

    fig, ax = plt.subplots(1, 1)

    speciesDict = {}

    print("HERE")

    Ts = list(range(1600, 4000, 100))
    for T in Ts:
        speciesList = ATKPRDS22.balanceAt(T, False)

        for entry in speciesList:
            specie, pop = entry[0], entry[1]
            if specie in speciesDict:
                speciesDict[specie].append(pop)
            else:
                speciesDict.update({specie: [pop]})

    for specie, pops in speciesDict.items():
        ax.plot(Ts, pops, label=specie)

    ax.set_yscale("log")
    ax.set_ylim((1e-6, 1))
    labelLines(ax.get_lines())

    plt.show()

    """

    """
    NG = Ingredient.getLine(693)
    NC1316 = Ingredient.nitrocellulose(0.1316)
    NC1315 = Ingredient.nitrocellulose(0.1315)
    NC1311 = Ingredient.nitrocellulose(0.1311)
    DBP = Ingredient.find("DIBUTYL PHTHALATE")
    DNT = Ingredient.find("DINITRO TOLUENE")
    DPA = Ingredient.find("2 NITRO DIPHENYL AMINE")

    KNO3 = Ingredient.find("POTASSIUM NITRATE")
    K2SO4 = Ingredient.find("POTASSIUM SULFATE")

    ETC = Ingredient(
        name="Ethyl Centralite",
        elements={"C": 17, "H": 20, "O": 1, "N": 2},
        rho=1.140,
        rho_u="g/cc",
        Hf=-391.5,
        Hf_u="J/g",
    )
    MEC = Ingredient(
        name="Methyl Centralite",
        elements={"C": 15, "H": 16, "O": 1, "N": 2},
        rho=1.2,
        rho_u="g/cc",
        Hf=-73.1,
        Hf_u="kJ/mol",
    )

    WC846 = Mixture(
        name="WC846",
        compoDict={NC1316: 81.40, NG: 10.39, DBP: 5.61, DPA: 0.97, DNT: 0.06},
        Delta=0.2,
    )

    WC846.prettyPrint()

    WC870 = Mixture(
        name="WC870",
        compoDict={
            NC1311: 79.70,
            NG: 9.94,
            DBP: 5.68,
            DPA: 0.95,
            KNO3: 0.78,
            DNT: 0.69,
        },
        Delta=0.2,
    )
    WC870.prettyPrint()

    IMR4227 = Mixture(
        name="IMR4227",
        compoDict={NC1315: 89.72, DNT: 6.74, DPA: 0.84, K2SO4: 0.73},
        Delta=0.2,
    )
    IMR4227.prettyPrint()

    IMR4350 = Mixture(
        name="IMR4350",
        compoDict={NC1315: 91.56, DNT: 4.80, DPA: 0.85, K2SO4: 0.73},
        Delta=0.2,
    )
    IMR4350.prettyPrint()

    IMR8138M = Mixture(
        name="IMR8138M",
        compoDict={NC1315: 92.96, ETC: 3.60, DPA: 0.89, K2SO4: 0.69},
        Delta=0.2,
    )
    IMR8138M.prettyPrint()

    IMR4895 = Mixture(
        name="IMR4895",
        compoDict={NC1315: 88.99, DNT: 7.55, DPA: 0.78, K2SO4: 0.93},
        Delta=0.2,
    )
    IMR4895.prettyPrint()

    IMR5010 = Mixture(
        name="IMR5010",
        compoDict={NC1315: 88.27, DNT: 8.76, K2SO4: 0.64, DPA: 0.55},
        Delta=0.2,
    )
    IMR5010.prettyPrint()

    CMR160 = Mixture(
        name="CMR160",
        compoDict={NC1315: 95.77, MEC: 2.43, DPA: 0.73, K2SO4: 0.84},
    )
    CMR160.prettyPrint()

    for T in (2855, 2782, 2878, 2927, 2816, 2874, 2869, 2776):
        print(f"{estU8053(T):.5g}")
"""
