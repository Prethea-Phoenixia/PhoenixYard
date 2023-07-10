import difflib

from hunt import balance


class Ingredient:
    allIngr = {}

    def __init__(self, name, elements, Hf, rho, flag=""):
        self.name = name
        self.flag = flag
        self.elements = elements
        self.Hf = Hf
        self.rho = rho

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
                        element, nbr = (
                            int(line[39 + i * 5 : 42 + i * 5].strip()),
                            line[42 + i * 5 : 44 + i * 5].strip(),
                        )
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
        """
        return an Ingredient object that is closest to the supplied common name,
        or alternative (abbreviated) name.
        """

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

if __name__ == "__main__":
    Ingredient.readFile("data/PEPCODED.DAF")
    # print(*Ingredient.allIngr.keys(), sep="\n")


"""
Major chemical equilibriums include:

0.5 O2 -> O
Kp1 = pO / pO2 ^ (1/2)

0.5 H2 -> H
Kp2 = pH / pH2 ^ (1/2)

0.5 H2 + 0.5 O2 -> OH
Kp3 = pOH / (pO2 ^ 1/2 pH2 ^ 1/2)

H2 + 0.5 O2 -> H2O
Kp4 = pH2O / (pH2 pO2 ^ 1/2)

0.5 N2 -> N
Kp5 = pN / pN2 ^ (1/2)

0.5 N2 + 0.5 O2 -> NO
Kp6 = pNO / (pN2 pO2) ^ 1/2

C -> C
Kp7 = pC(c,diamond) / pC(c,graphite)

C -> C
Kp8 = pC (g) / pC

C + 0.5 O2 -> CO
Kp9 = pCO / (pC pO2 ^ 1/2)

C + O2 -> CO2
Kp10 = pCO2 / (pC pO2)

C + H4 -> CH4
Kp11 = pCH4 / (pC (pH2)^2)

2C + H2 -> C2H2
Kp12 = pC2H2 / ((pC)^2 pH2)

0.5 CL2 -> CL
Kp13 = pCl / pCl2 ^ 1/2

0.5 H2 + 0.5 CL2 -> HCL
Kp14 = pHCl / (pH2 ^ 1/2 pCl2 ^ 1/2)

0.5 N2 + 0.5 H2 -> NH
Kp15 = pNH / (pN2 ^ 1/2 pH2 ^ 1/2)

0.5 N2 + 1.5 H2 -> NH3
Kp16 = pNH3 / (pN2 ^ 1/2 pH2 ^ 3/2)

C + H2 -> CH2
Kp17 = pCH2 / (pC pH2)

NO + 0.5 O2 -> NO2
Kp18 = pNO2 / (pNO pO2 ^ 1/2)

0.5 N2 + 0.5 O2 -> NO
Kp19 = pNO / ( pN2 ^ 1/2 pO2 ^ 1/2)
"""
