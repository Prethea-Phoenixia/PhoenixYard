# low: 200-1000K
# high: 1000K - 6000K
CO2_a_low = [
    0.46365111e01,
    0.27414569e-02,
    -0.99589759e-06,
    0.16038666e-09,
    -0.91619857e-14,
    -0.49024904e05,
    -0.19348955e01,
]

CO2_a_high = [
    0.23568130e01,
    0.89841299e-02,
    -0.71220632e-05,
    0.24573008e-08,
    -0.14288548e-12,
    -0.48371971e05,
    0.99009035e01,
]

CO2_dfH298 = -0.47328105e0


if __name__ == "__main__":
    R = 1.987207  # cal/molK

    def f(T):
        if T < 200:
            raise ValueError("Temp too low (<200K)")
        elif T < 1000:
            a1, a2, a3, a4, a5, a6, a7 = CO2_a_low
        elif T < 6000:
            a1, a2, a3, a4, a5, a6, a7 = CO2_a_high
        else:
            raise ValueError("Temp too high (>6000K)")

        return (R * T) * (
            a1
            + a2 * T / 2
            + a3 * T**2 / 3
            + a4 * T**3 / 4
            + a5 * T**4 / 5
            + a6 / T
        )

    print((f(3000) - f(300)) / (3000 - 300))
