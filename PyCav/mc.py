import bubble_model as bm
import bubble_state as bs
import time_advancer as adv
import numpy as np


class mc:
    def __init__(self, config=None):

        self.adv_config = config["advancer"]
        self.pop_config = config["pop"]
        self.model_config = config["model"]

        self.advancer = adv.time_advancer(self.adv_config)
        self.state = bs.bubble_state(
            pop_config=self.pop_config, model_config=self.model_config
        )

    def moment(self, sample=[], mom=[]):
        N = len(sample)
        if self.state.num_RV_dim == 2:
            return 1.0 / N * np.sum((sample[:, 0] ** mom[0]) * (sample[:, 1] ** mom[1]))
        else:
            raise NotImplementedError

    def get_sample(self):
        if self.state.shape == "lognormal":
            self.sample = np.random.lognormal(
                np.log(self.state.muR0), self.state.sigR0, self.N
            )
        elif self.state.shape == "normal":
            self.sample = np.random.normal(self.state.muR0, self.state.sigR0, self.N)
        else:
            raise NotImplementedError

    def simulate_sample(self, N=0):
        self.N = N
        T = self.advancer.T
        p = self.advancer.p

        self.get_sample()
        for R0 in self.sample:
            bubble = bm.bubble_model(config=self.model_config, R0=R0)
            bubble.solve(T=T, p=p, Ro=R0, Vo=0.0)


# if __name__ == "__main__":
#     monte = mc()

# plt.plot(sol.t,sol.y[0])
# plt.show()
