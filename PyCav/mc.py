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

    def get_sample(self):
        if self.state.shape == "lognormal":
            self.sample = np.random.lognormal(
                np.log(self.state.muR0), self.state.sigR0, self.Nmc
            )
        elif self.state.shape == "normal":
            self.sample = np.random.normal(self.state.muR0, self.state.sigR0, self.Nmc)
        else:
            raise NotImplementedError

    def simulate_sample(self, Nmc=1, Nt=10):
        self.Nmc = Nmc
        self.Nt = Nt
        T = self.advancer.T
        p = self.advancer.p
        ts = np.linspace(0, T, num=Nt)

        self.get_sample()
        # print('sample = ',self.sample)
        sols = []
        for R0 in self.sample:
            bubble = bm.bubble_model(config=self.model_config, R0=R0)
            sol = bubble.solve(T=T, p=p, Ro=R0, Vo=0.0, ts=ts)
            sols.append(sol)

        return sols
