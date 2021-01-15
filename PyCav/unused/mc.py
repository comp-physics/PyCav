import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt


class integrate:
    def __init__(self):
        pass


def fun(t, y):
    return np.array([y[1] * t, 2 * t * y[0]])


if __name__ == "__main__":
    t0 = 0.0
    T = 1.0
    y0 = np.array([1.0, 0.0])
    ret = sp.solve_ivp(fun, (t0, T), y0, method="RK45")
    print(ret.t)
    print(ret.y)

    # plt.plot(ret.t,ret.y[0])
    # plt.plot(ret.t,ret.y[1])
    # plt.show()
