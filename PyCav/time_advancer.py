import bubble_model as bm
import bubble_state as bs
import numpy as np


class time_advancer:
    def __init__(self, config={}):

        if "T" in config:
            self.T = config["T"]
            if self.T <= 0:
                raise ValueError(self.T)
        else:
            raise Exception("No final time T")

        if "dt" in config:
            self.dt = config["dt"]
            if self.dt <= 0.0:
                raise ValueError(self.dt)
        else:
            raise Exception("No time step dt")

        if "method" in config:
            self.method = config["method"]
        else:
            self.method = "Euler"

        if "p" in config:
            self.p = config["p"]
        else:
            raise Exception("No pressure p")

        if self.method == "Euler":
            self.advance = self.euler
            self.n_stages = 1
        else:
            raise NotImplementedError

        return

    def initialize_state(self, pop_config=None, model_config=None):

        self.state = bs.bubble_state(pop_config=pop_config, model_config=model_config)

    def euler(self):

        # self.stage_state[0] = self.state.copy()
        # self.bubble_mgr.compute_rhs(self.stage_state[0], self.stage_k[0])
        # self.state = self.stage_state[0] + self.time_step * self.stage_k[0]

        self.state.get_rhs(self.p)
        self.state.vals += self.dt * self.state.rhs
        # print('vals: ', self.state.vals)
        # print('rhs:  ', self.state.rhs )

        # print(self.state.vals)

        return

    def run(self):
        self.time = 0.0
        i_step = 0
        step = True
        self.save = []
        self.times = []
        # its=10
        while step == True:
            self.times.append(self.time)
            self.save.append(self.state.vals[0].copy())

            # self.save.append(self.state.vals.copy())
            self.advance()
            i_step += 1
            self.time += self.dt

            if self.time >= self.T:
                step = False
            # if i_step % its == 0:

        # self.times.append(self.time)
        # self.save.append(self.state.vals[0].copy())
        # print('save: ', self.save)
        self.save = np.array(self.save, dtype=np.float32)
        self.plot()

    def plot(self):
        import matplotlib.pyplot as plt

        # plt.plot(self.save[:,0],self.save[:,1])
        plt.plot(self.times, self.save[:, 0])
        plt.xlabel("$t$")
        plt.ylabel("$R(t)$")
        plt.show()


if __name__ == "__main__":

    config = {}
    config["method"] = "Euler"
    config["dt"] = 1.0e-5
    config["T"] = 30.0
    config["p"] = 1.0
    myadv = time_advancer(config=config)
