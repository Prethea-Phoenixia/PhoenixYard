from enum import Enum


def _30SIMN2MOVA_Y(T):
    return (
        -0.047 * (T / 293) ** 4
        + 0.201 * (T / 293) ** 3
        - 0.15 * (T / 293) ** 2
        - 0.376 * (T / 293)
        + 1.375
    ) * 900e6  # Circumferential Yield Strength in Pa


def _30SIMN2MOVA_E(T):
    return (
        -0.025 * (T / 293) ** 2 + 0.04 * (T / 293) + 0.985
    ) * 213e9  # Youngs Modulus in Pa,


class Material(Enum):
    def __init__(self, desc, rho, Y, E):
        self.desc = desc
        self.rho = rho
        self.Y = Y
        self.E = E

    def __str__(self):
        return self.desc

    _30SIMN2MOVA = (
        "30SiMn2MoVA",
        7801,
        _30SIMN2MOVA_Y,  # Circumferential Yield Strength in Pa
        _30SIMN2MOVA_E,
    )


MATERIALS = {i.desc: i for i in Material}


if __name__ == "__main__":
    print(Material._30SIMN2MOVA)
    print(Material._30SIMN2MOVA.Y(293))
    print(Material._30SIMN2MOVA.E(293))
