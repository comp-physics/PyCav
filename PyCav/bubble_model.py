import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt


class bubble_model:
    def __init__(self, config={}, R0=1.0):

        self.config = config
        self.R0 = R0
        self.parse_config()

        if self.model == "RPE":
            self.num_RV_dim = 2
            self.state = np.array([self.R, self.V])
        else:
            raise NotImplementedError

    def parse_config(self):
        if "model" in self.config:
            self.model = self.config["model"]
        else:
            self.model = "RPE"

        if "R" in self.config:
            self.R = self.config["R"]
        else:
            self.R = self.R0

        if "V" in self.config:
            self.V = self.config["V"]
        else:
            self.V = 0.0

        if "gamma" in self.config:
            self.gamma = self.config["gamma"]
        else:
            self.gamma = 1.4

        if "Ca" in self.config:
            self.Ca = self.config["Ca"]
        else:
            self.Ca = 1.0

        if "Re_inv" in self.config:
            self.Re_inv = self.config["Re_inv"]
            if self.Re_inv <= 0.0:
                raise ValueError(self.Re_inv)
            else:
                self.viscosity = True
        else:
            self.viscosity = False
            self.Re_inv = 0

        if "Web" in self.config:
            self.Web = self.config["Web"]
            if self.Web <= 0.0:
                raise ValueError(self.Web)
            else:
                self.tension = True
        else:
            self.tension = False
            self.Web = 0

    def pbw(self):
        self.cpbw = self.Ca * ((self.R0 / self.R) ** (3.0 * self.gamma)) - self.Ca + 1.0
        if self.tension:
            self.cpbw += (
                2.0 / (self.Web * self.R0) * (self.R0 / self.R) ** (3.0 * self.gamma)
            )

    def rpe(self, p):
        self.pbw()
        rhs = -1.5 * self.V ** 2.0 + (self.cpbw - p) / self.R
        if self.viscosity:
            rhs -= 4.0 * self.Re_inv * self.V / (self.R ** 2.0)
        return [self.V, rhs]

    def rhs(self, p):
        # print(self.state)
        self.update_state()
        if self.model == "RPE":
            rhs = self.rpe(p)
        else:
            raise NotImplementedError
        return rhs

    def update_state(self):
        if self.model == "RPE":
            self.R = self.state[0]
            self.V = self.state[1]
        else:
            raise NotImplementedError

    def wrap(self, t, y):
        self.R = y[0]
        self.V = y[1]
        return np.array(self.rpe(self.p(t)))

    def solve(self, T=0, Ro=1.0, Vo=0.0, p=1.0, ts=None):
        self.p = p
        y0 = np.array([Ro, Vo])
        if ts is None:
            ret = sp.solve_ivp(self.wrap, (0.0, T), y0, method="LSODA", rtol=1e-3)
        else:
            ret = sp.solve_ivp(
                self.wrap, (0.0, T), y0, method="LSODA", rtol=1e-3, t_eval=ts
            )

        return ret


if __name__ == "__main__":

    config = {}
    config["model"] = "RPE"
    # config["R"] = 1.0
    config["V"] = 0.0
    config["gamma"] = 1.4
    config["Ca"] = 1.0
    # config["Re_inv"] = 0
    # config["Web"] = 0

    R0 = 0.3
    mybub = bubble_model(config=config, R0=R0)
    # rhs = mybub.rhs(p=1.1)

    sol = mybub.solve(T=10.0, p=1.1)
    plt.plot(sol.t, sol.y[0])
    plt.show()
