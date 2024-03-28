from ballistics.therm import Ingredient, Mixture

if __name__ == "__main__":

    """Material Nitrocellulose
    Elements                        C   H   N   O
    ElementCount                    6   7   3   11
    Density_kg__m3                  1490
    EnthalpyOfFormation_kJ__mol     -534.732
    SpecificHeat_J__kg_K            1406.7
    ThermalConductivity_W__m_K      .424
    YoungsModulus_GPa               .029
    ShearModulus_GPa                .010

    """
    NC = Ingredient(
        name="Nitrocellulose",
        elements={"C": 6, "H": 7, "N": 3, "O": 11},
        rho=1.490,
        rho_u="g/cc",
        Hf=-534.732,
        Hf_u="kJ/mol",
    )
    Nitrocellulose = Mixture("Nitrocellulose", compoDict={NC: 1}, Delta=0.2)
    Nitrocellulose.prettyPrint()
    """
    Material Octogen
    Elements                        C   H   N   O
    ElementCount                    4   8   8   8
    Density_kg__m3                  1910
    EnthalpyOfFormation_kJ__mol     104.77
    MeltingPoint_K                  549
    SpecificHeat_J__kg_K            1260
    ThermalConductivity_W__m_K      .379
    YoungsModulus_GPa               15.3
    ShearModulus_GPa                5.795
    RefractiveIndex                 Octogen
    RoughnessCoefficient            .9
    """
    HMX = Ingredient(
        name="HMX",
        elements={"C": 4, "H": 8, "N": 8, "O": 8},
        rho=1.91,
        rho_u="g/cc",
        Hf=104.77,
        Hf_u="kJ/mol",
    )
    Octogen = Mixture("Octogen", compoDict={HMX: 1}, Delta=0.2)
    Octogen.prettyPrint()

    """
    Material Cyclonite
    Elements                        C   H   N   O
    ElementCount                    3   6   6   6
    Density_kg__m3                  1820
    EnthalpyOfFormation_kJ__mol     66.5
    MeltingPoint_K                  478.6
    BoilingPoint_K                  507
    SpecificHeat_J__kg_K            1120.56
    ThermalConductivity_W__m_K      .227
    YoungsModulus_GPa               11.99
    ShearModulus_GPa                9.26
    RefractiveIndex                 Octogen
    RoughnessCoefficient            .9
    """

    RDX = Ingredient(
        name="RDX",
        elements={"C": 3, "H": 6, "N": 6, "O": 6},
        rho=1.82,
        rho_u="g/cc",
        Hf=66.5,
        Hf_u="kJ/mol",
    )
    Cyclonite = Mixture("Cyclonite", compoDict={RDX: 1}, Delta=0.2)
    Cyclonite.prettyPrint()

    """

    Material Nitroglycerin
    Elements                        C   H   N   O
    ElementCount                    3   5   3   9
    Density_kg__m3                  1600
    EnthalpyOfFormation_kJ__mol     -364
    MeltingPoint_K                  287
    BoilingPoint_K                  323
    SpecificHeat_J__kg_K            1778.6
    ThermalConductivity_W__m_K      .21964
    YoungsModulus_GPa               35
    ShearModulus_GPa                16
    Viscosity_Pa_s                  .036
    RefractiveIndex                 Octogen
    RoughnessCoefficient            .9

    """

    # NG = Ingredient(
    #     name="Nitroglycerin",
    #     elements={"C": 3, "H": 5, "N": 3, "O": 9},
    #     rho=1.6,
    #     rho_u="g/cc",
    #     Hf=-364,
    #     Hf_u="kJ/mol",
    # )
    # Nitroglycerin = Mixture("Nitroglycerin", compoDict={NG: 1}, Delta=0.2)
    # Nitroglycerin.prettyPrint()

    """
    Material TNT
    Elements                        C   H   N   O
    ElementCount                    7   5   3   6
    Density_kg__m3                  1654
    EnthalpyOfFormation_kJ__mol     -65.5
    MeltingPoint_K                  353.50
    BoilingPoint_K                  513
    SpecificHeat_J__kg_K            1377
    ThermalConductivity_W__m_K      .219
    YieldStrength_MPa               23.1
    YoungsModulus_GPa               2.551
    ShearModulus_GPa                .8963
    """

    TNT = Ingredient(
        name="Trinitrotoluene",
        elements={"C": 7, "H": 5, "N": 3, "O": 6},
        rho=1.654,
        rho_u="g/cc",
        Hf=-65.5,
        Hf_u="kJ/mol",
    )
    Trinitrotoluene = Mixture("Trinitrotoluene", compoDict={TNT: 1}, Delta=0.2)
    Trinitrotoluene.prettyPrint()

    """
    Material PETN
    Elements                        C   H   N   O
    ElementCount                    5   8   4   12
    Density_kg__m3                  1770
    EnthalpyOfFormation_kJ__mol     -538.5
    SpecificHeat_J__kg_K            1000
    ThermalConductivity_W__m_K      .143
    MeltingPoint_K                  414.4
    YoungsModulus_GPa               12.82
    ShearModulus_GPa                9.37
    RoughnessCoefficient            .9
    """

    # PETN = Ingredient(
    #     "Pentaerythritol tetranitrate",
    #     elements={"C": 5, "H": 8, "N": 4, "O": 12},
    #     rho=1.770,
    #     rho_u="g/cc",
    #     Hf=-538.5,
    #     Hf_u="kJ/mol",
    # )

    # Penthrite = Mixture("Penthrite", compoDict={PETN: 1}, Delta=0.2)
    # Penthrite.prettyPrint()
