import bubble_model as bm
import numpy as np

class bubble_state:

    def __init__(self,
        state_config={},
        model_config={}
        ):

        self.model_config = model_config

        if "NR0" in state_config:
            self.NR0 = state_config["NR0"]
        else:
            self.NR0 = 1

        if "shape" in state_config:
            self.shape = state_config["shape"]
        else:
            self.shape = "lognormal"

        if "binning" in state_config:
            self.binning= state_config["binning"]
        else:
            self.binning= "Simpson"

        if "sigR0" in state_config:
            self.sigR0 = state_config["sigR0"]
        else:
            self.sigR0 = 0.3

        # Initiate bubbles, weights, abscissas
        if self.NR0 == 1:
            self.init_mono()
        elif self.NR0 > 1:
            if self.binning == "Simpson":
               self.init_simp() 
            else:
                raise NotImplementedError
        else:
            raise ValueError(NR0)

        # Assume all bubbles have the same model
        self.num_RV_dim = self.bubble[0].num_RV_dim
        self.vals = np.zeros((self.NR0,self.num_RV_dim))
        for i in range(self.NR0):
            self.vals[i,:] = self.bubble[i].state
            # Create view so that bubble state is
            # auto-updated when changing vals
            self.bubble[i].state = self.vals[i,:].view()

        self.rhs = np.zeros((self.NR0,self.num_RV_dim))

    def init_simp(self):
        raise NotImplementedError

    def init_mono(self):
        self.w = np.ones(1)
        self.R0 = np.ones(1)
        self.bubble = [ bm.bubble_model(config=self.model_config,R0=1) ]

    def get_rhs(self,p):
        for i in range(self.NR0):
            self.rhs[i,:] = self.bubble[i].rhs(p)

if __name__ == "__main__":

    state = bubble_state()
    val = state.vals
    p = 1.1
    dt = 0.1

    print('state = ',val)
    for i in range(5):
        state.get_rhs(p)
        val += dt * state.rhs
        print('state = ',val)
        # print('state2 = ', state.bubble[0].state)
        # print('rhs    = ',state.rhs)

