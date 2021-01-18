import bubble_model as bm
import numpy as np


class bubble_state:
    def __init__(self, pop_config={}, model_config={}):

        self.model_config = model_config

        if "NR0" in pop_config:
            self.NR0 = pop_config["NR0"]
        else:
            self.NR0 = 1

        if "shape" in pop_config:
            self.shape = pop_config["shape"]
        else:
            self.shape = "lognormal"

        if "binning" in pop_config:
            self.binning = pop_config["binning"]
        else:
            self.binning = "Simpson"

        if "sigR0" in pop_config:
            self.sigR0 = pop_config["sigR0"]
        else:
            self.sigR0 = 0.3

        if "muR0" in pop_config:
            self.muR0 = pop_config["muR0"]
        else:
            self.muR0 = 1.0

        if "moments" in pop_config:
            self.moments = pop_config["moments"]
        else:
            self.moments = [[0,0]]

        # Initiate bubbles, weights, abscissas
        if self.NR0 == 1:
            self.init_mono()
        elif self.NR0 > 1:
            if self.binning == "Simpson":
                self.init_simp()
            else:
                raise NotImplementedError
        else:
            raise ValueError(self.NR0)

        # Assume all bubbles have the same model
        self.num_RV_dim = self.bubble[0].num_RV_dim
        self.vals = np.zeros((self.NR0, self.num_RV_dim))
        for i in range(self.NR0):
            self.vals[i, :] = self.bubble[i].state
            # Create view so that bubble state is
            # auto-updated when changing vals
            self.bubble[i].state = self.vals[i, :].view()

        self.rhs = np.zeros((self.NR0, self.num_RV_dim))

    def init_simp(self):
        self.w = np.ones(1)
        self.R0 = np.ones(1)
        self.bubble = [bm.bubble_model(config=self.model_config, R0=1)]

        raise NotImplementedError

    def init_mono(self):
        self.w = np.ones(1)
        # self.R0 = np.ones(1)
        self.bubble = [bm.bubble_model(config=self.model_config, R0=1)]

    def get_rhs(self, p):
        for i in range(self.NR0):
            self.rhs[i, :] = self.bubble[i].rhs(p)

    def quad(self, mom):
        # mom : the moment we want to compute
        # should be of the same length as the
        # number of variables (e.g. R,V)
        ret = 0.0
        for i in range(self.NR0):
            change = self.w[i]
            for j in range(self.num_RV_dim):
                change *= self.vals[i, j] ** mom[j]
            ret += change
        return ret

    def moment(self, sample=[]):
        Nmc = len(sample)
        Nmom = len(self.moments)
        Nt = len(sample[0].y[0])
        ret = np.zeros((Nmom,Nt))
        for k, mom in enumerate(self.moments):
            if self.num_RV_dim == 2:
                for samp in sample:
                    ret[k,:] += samp.y[0] ** mom[0] * samp.y[1] ** mom[1]
            else:
                raise NotImplementedError

        ret = ret / Nmc
        return ret


if __name__ == "__main__":

    state = bubble_state()
    val = state.vals
    p = 1.1
    dt = 0.1

    print("state = ", val)
    for i in range(5):
        state.get_rhs(p)
        val += dt * state.rhs
        print("state = ", val)
