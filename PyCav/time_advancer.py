import bubble_model as bm
import bubble_state as bs
import numpy as np

class time_advancer:

    def __init__(self,config=None):

        self.method = config["method"]
        self.dt = config["dt"]
        self.T = config["T"]

        if self.method == "Euler":
            self.advance = self.euler
            self.n_stages = 1
        else:
            raise NotImplementedError

        if self.dt <= 0.:
            raise ValueError(self.dt)

        if self.T <= 0:
            raise ValueError(self.T)

        return

    def initialize_state(self,config=None):

        model = config["model"]
        NR0 = config["NR0"]
        shape = config["shape"]
        binning = config["binning"]
        sigR0 = config["sigR0"]

        self.state = bs.bubble_state(
                model=model,
                NR0=NR0,
                shape=shape,
                binning=binning,
                sigR0=sigR0)

        return

    def euler(self):

        # self.stage_state[0] = self.state.copy()
        # self.bubble_mgr.compute_rhs(self.stage_state[0], self.stage_k[0])
        # self.state = self.stage_state[0] + self.time_step * self.stage_k[0]

        p = 1.1
        self.state.get_rhs(p)
        print('vals: ', self.state.vals)
        # print('rhs: ', self.state.rhs)
        # self.state.update_vals(self.dt*self.state.rhs)
        self.state.vals += self.dt * self.state.rhs

        return

    def run(self):

        self.time = 0.0

        i_step = 0
        step = True
        while step == True:

            self.advance()

            i_step += 1
            self.time += self.dt

            if self.time > self.T:
                step = False

        return

if __name__ == "__main__":

    config = {}
    config["advancer"] = {}
    config["advancer"]["method"] = "Euler"
    config["advancer"]["dt"] = 1.0e-5
    config["advancer"]["T"] = 30.0
    myadv = time_advancer(config["advancer"])

