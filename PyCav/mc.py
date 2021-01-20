import bubble_model as bm
import bubble_state as bs
import time_advancer as adv
import waveforms as wf
import numpy as np
import matplotlib.pyplot as plt


class mc:
    def __init__(self, config=None):
        self.adv_config = config["advancer"]
        self.pop_config = config["pop"]
        self.model_config = config["model"]
        self.wave_config = config["wave"]
        self.mc_config = config["mc"]

        if "Ntimes" in self.mc_config:
            self.Nt = self.mc_config["Ntimes"]
        else:
            raise Exception("No number of output times")

        if "Nsamples" in self.mc_config:
            self.Nmc = self.mc_config["Nsamples"]
        else:
            raise Exception("No number of samples")

        self.advancer = adv.time_advancer(self.adv_config)
        self.state = bs.bubble_state(
            pop_config=self.pop_config, model_config=self.model_config
        )
        self.wave = wf.waveforms(config=self.wave_config)

    def get_sample(self):
        if self.state.shape == "lognormal":
            self.sample = np.random.lognormal(
                np.log(self.state.muR0), self.state.sigR0, self.Nmc
            )
        elif self.state.shape == "normal":
            self.sample = np.random.normal(self.state.muR0, self.state.sigR0, self.Nmc)
        else:
            raise NotImplementedError

    def moment(self, sample=[]):
        ret = np.zeros((self.state.Nmom, self.Nt))
        for k, mom in enumerate(self.state.moments):
            if self.state.num_RV_dim == 2:
                for samp in sample:
                    ret[k, :] += samp.y[0] ** mom[0] * samp.y[1] ** mom[1]
            else:
                raise NotImplementedError

        return ret / self.Nmc

    def run(self):
        T = self.advancer.T
        p = self.wave.p
        ts = np.linspace(0, T, num=self.Nt)

        self.get_sample()
        sols = []
        for R0 in self.sample:
            bubble = bm.bubble_model(config=self.model_config, R0=R0)
            sol = bubble.solve(T=T, p=p, Ro=R0, Vo=0.0, ts=ts)
            sols.append(sol)

        moments = self.moment(sols)
        fig, ax = plt.subplots(1, self.state.Nmom)
        for i in range(self.state.Nmom):
            ax[i].plot(sols[i].t, moments[i])
            ax[i].set(xlabel="$t$", ylabel="$M$" + str(self.state.moments[i]))
