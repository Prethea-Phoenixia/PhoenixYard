from sklearn import linear_model
import numpy as np

acCHNO = [
    [21.43, 31.18, 10.04, 35.65],
    [21.56, 28.00, 10.36, 35.48],
    [21.57, 30.33, 10.12, 35.53],
    [22.32, 30.10, 10.46, 34.69],
    [19.56, 27.21, 10.78, 36.66],
    [24.65, 32.36, 10.08, 33.13],
    [24.65, 32.36, 10.08, 33.13],
    [24.28, 32.99, 9.57, 33.81],
    [16.45, 35.71, 25.93, 25.20],
    [14.74, 33.84, 26.07, 26.47],
    [16.52, 36.10, 25.91, 25.14],
    [22.27, 27.91, 8.97, 36.17],
    [24.66, 28.88, 9.23, 34.09],
    [25.38, 29.93, 9.05, 33.65],
]

acCO = [[v[0], v[-1]] for v in acCHNO]
cook = [[6.9], [-11.5]]

bt = [
    940,
    910,
    934,
    957,
    866,
    1026,
    1018,
    1021,
    1058,
    1001,
    1066,
    893,
    978,
    1005,
]

if __name__ == "__main__":
    X = np.array(acCHNO)
    Y = np.array(bt)
    lr = linear_model.LinearRegression().fit(X, Y)
    print(lr.score(X, Y))
    print(lr.coef_)
    print(lr.intercept_)

    print(lr.predict(X))
    """
    X = np.array(acCO)
    Y = np.array(bt)
    lr = linear_model.LinearRegression().fit(X, Y)
    print(lr.score(X, Y))
    print(lr.coef_)
    print(lr.intercept_)
    """
