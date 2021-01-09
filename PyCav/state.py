import numpy as np
import bubble_model as bub

class state:

    def __init__(self,
        model="RPE",
        NR0=1,
        shape="lognormal",
        sigR0=0.3,
        binning="simpson"
        ):

        self.model = model
        self.NR0 = NR0
        self.shape = shape
        self.sigR0 = sigR0
        self.binning = binning

        # Initiate bubbles, weights, abscissas
        if self.NR0 == 1:
            self.init_mono()
        elif self.NR0 > 1:
            if self.binning == "simpson":
               self.init_simp() 
            else:
                raise NotImplementedError
        else:
            raise ValueError(NR0)

        # Assume all bubbles have the same model
        self.num_RV_dim = self.bubble[0].num_RV_dim
        self.vals = np.zeros((self.NR0,self.num_RV_dim))
        for i in range(self.NR0):
            self.vals[i,:] = self.bubble[i].state()

    def init_simp(self):
        raise NotImplementedError

    def init_mono(self):
        self.w = np.ones(1)
        self.R0 = np.ones(1)
        self.bubble = [ bub.bubble_model(model=self.model,R0=1) ]


if __name__ == "__main__":

    stat = state()
    val = stat.vals
    print('state = ',val)
