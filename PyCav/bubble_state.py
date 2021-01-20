import bubble_model as bm
import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt
from sys import exit


class bubble_state:

    def __init__(self, pop_config={}, model_config={}):

        self.model_config = model_config
        self.pop_config = pop_config
        self.parse_config()

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

        for mom in self.moments:
            if len(mom) != self.num_RV_dim:
                raise ValueError(mom,self.num_RV_dim)

        self.vals = np.zeros((self.NR0, self.num_RV_dim))
        for i in range(self.NR0):
            self.vals[i, :] = self.bubble[i].state
            # Create view so that bubble state is
            # auto-updated when changing vals
            self.bubble[i].state = self.vals[i, :].view()

        self.rhs = np.zeros((self.NR0, self.num_RV_dim))

    def parse_config(self):

        if "NR0" in self.pop_config:
            self.NR0 = self.pop_config["NR0"]
        else:
            self.NR0 = 1

        if "shape" in self.pop_config:
            self.shape = self.pop_config["shape"]
        else:
            self.shape = "lognormal"

        if "binning" in self.pop_config:
            self.binning = self.pop_config["binning"]
        else:
            self.binning = "Simpson"

        if "sigR0" in self.pop_config:
            self.sigR0 = self.pop_config["sigR0"]
        else:
            self.sigR0 = 0.3

        if "muR0" in self.pop_config:
            self.muR0 = self.pop_config["muR0"]
        else:
            self.muR0 = 1.0

        if "moments" in self.pop_config:
            self.moments = self.pop_config["moments"]
            self.Nmom = len(self.moments)
        else:
            self.moments = [[0, 0]]
            self.Nmom = 1

    def init_pdf(self):
        if self.shape == "lognormal":
            self.f = sp.lognorm.pdf(self.R0,self.sigR0,loc=np.log(self.muR0))
        elif self.shape == "normal":
            self.f = sp.norm.pdf(self.R0,self.sigR0,loc=self.muR0)
        else:
            raise NotImplementedError

    def init_simp(self):
        a = 0.8*np.exp(-2.8*self.sigR0)
        b = 0.2*np.exp( 9.5*self.sigR0)+1.
        # dR0 = (b-a)/(self.NR0-1.)
        self.R0 = np.logspace(np.log10(a),np.log10(b),num=self.NR0)

        self.dR0 = np.zeros(self.NR0)
        for i in range(self.NR0-1):
            self.dR0[i] = self.R0[i+1] - self.R0[i] 
        self.dR0[self.NR0-1] = self.dR0[self.NR0-2]
        
        self.init_pdf()

        self.w = np.zeros(self.NR0)
        for i in range(self.NR0):
            if i == 0 or i == self.NR0-1: 
                self.w[i] = 1./3.
            elif i % 2 == 0:
                self.w[i] = 2./3.
            else:
                self.w[i] = 4./3.

        self.w *= self.dR0
        self.w *= self.f

        # print(np.sum(self.w))
        # print(np.sum(self.f*self.dR0))
        # plt.plot(self.R0,self.f)
        # plt.xscale("log")
        # plt.show()

        self.bubble = [ bm.bubble_model(config=self.model_config, R0=x) for x in self.R0 ]

    def init_mono(self):
        self.w = np.ones(1)
        self.R0 = np.ones(1)
        self.bubble = [bm.bubble_model(config=self.model_config, R0=1)]

    def get_rhs(self, state, p):
        self.vals[:,:] = state  
        for i in range(self.NR0):
            self.rhs[i, :] = self.bubble[i].rhs(p)
        return self.rhs

    def get_quad(self):
        # should be of the same length as the
        # number of variables (e.g. R,V)
        ret = np.zeros(self.Nmom)
        for k, mom in enumerate(self.moments):
            # for i in range(self.NR0):
            if self.num_RV_dim == 2:
                ret[k] = np.sum(self.w[:] *
                        self.vals[:,0] ** mom[0] *
                        self.vals[:,1] ** mom[1])
        return ret

    def moment(self, sample=[]):
        Nmc = len(sample)
        Nt = len(sample[0].y[0])
        ret = np.zeros((self.Nmom, Nt))
        for k, mom in enumerate(self.moments):
            if self.num_RV_dim == 2:
                for samp in sample:
                    ret[k, :] += samp.y[0] ** mom[0] * samp.y[1] ** mom[1]
            else:
                raise NotImplementedError

        ret = ret / Nmc
        return ret


if __name__ == "__main__":

    pop_config = {}
    pop_config["NR0"] = 101
    pop_config["shape"] = "lognormal"
    pop_config["binning"] = "Simpson"
    pop_config["muR0"] = 1.0
    pop_config["sigR0"] = 0.1
    pop_config["moments"] = [[0, 0], [1, 0], [0, 1]]

    model_config = {}
    model_config["model"] = "RPE"
    model_config["V"] = 0.0
    model_config["gamma"] = 1.4
    model_config["Re_inv"] = 0.001

    state = bubble_state(
            pop_config=pop_config,
            model_config=model_config)

    val = state.vals
    # print("state = ", val)
    plt.plot(state.R0,state.w)
    plt.show()

    # p = 1.1
    # dt = 0.001
    # for i in range(5):
    #     state.get_rhs(p)
    #     val += dt * state.rhs
    #     print("state = ", val)
