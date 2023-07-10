import difflib

# fmt: off
#the E stand for electron and D for deuterium*/
molarMasses = {
  "H": 1.00794,   "HE": 4.002602,"LI": 6.941,    "BE": 9.012182,"B": 10.811,    "C": 12.0107,
  "N": 14.00674,  "O": 15.9994,  "F": 18.9984032,"NE": 20.11797,"NA": 22.989770,"MG": 24.305,
  "AL": 26.981538,"SI": 28.0855, "P": 30.973761, "S": 32.066,   "CL": 35.4527,  "AR": 39.948,
  "K": 39.0983,   "CA": 40.078,  "SC": 44.95591, "TI": 47.88,   "V": 50.9415,   "CR": 51.996,
  "MN": 54.938,   "FE": 55.847,  "CO": 58.9332,  "NI": 58.6934, "CU": 63.546,   "ZN": 65.39,
  "GA": 69.723,   "GE": 72.61,   "AS": 74.9216,  "SE": 78.96,   "BR": 79.904,   "KR": 83.80,
  "RB": 85.4678,  "SR": 87.62,   "Y": 88.9059,   "ZR": 91.224,  "NB": 92.9064,  "MO": 95.94,
  "TC": 98.0,     "RU": 101.07,  "RH": 102.9055, "PD": 106.42,  "AG": 107.868,  "CD": 112.41,
  "IN": 114.82,   "SN": 118.71,  "SB": 121.757,  "TE": 127.60,  "I": 126.9045,  "XE": 131.29,
  "CS": 132.9054, "BA": 137.33,  "LA": 138.9055, "CE": 140.12,  "PR": 140.9077, "ND": 144.24,
  "PM": 145.,     "SM": 150.36,  "EU": 151.965,  "GD": 157.25,  "TB": 158.9253, "DY": 162.50,
  "HO": 164.9303, "ER": 167.26,  "TM": 168.9342, "YB": 173.04,  "LU": 174.967,  "HF": 178.49,
  "TA": 180.9479, "W": 183.85,   "RE": 186.207,  "OS": 190.2,   "IR": 192.22,   "PT": 195.08,
  "AU": 196.9665, "HG": 200.59,  "TL": 204.383,  "PB": 207.2,   "BI": 208.9804, "PO": 209.,
  "AT": 210.,     "RN": 222.,    "FR": 223.,     "RA": 226.0254,"AC": 227.,     "TH": 232.0381,
  "PA": 231.0359, "U": 238.029,  "NP": 237.0482, "PU": 244.,    "FM": 257.0,    "E": 0,
  "D": 2};
# fmt: on


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
