"""
Jinpeng Zhai, 翟锦鹏
2023 07 08
Contact: 914962409@qq.com

Periodic table, coded for use with the PEPCODED.DAF and NASA
thermo databases.
"""

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
    "U1": 4.002602, "U2": 9.012182, "U3": 24.305, "U4": 26.981538,
    "U5": 12.0107
}
# fmt: on
