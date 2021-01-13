import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt

class bubble_model:

    def __init__(self,
            config={},
            R0=1.):

        if "model" in config:
            self.model = config["model"]
        else:
            self.model = "RPE"

        self.R0 = R0

        if "R" in config:
            self.R = config["R"]
        else:
            self.R = self.R0

        if "V" in config:
            self.V = config["V"]
        else:
            self.V = 0.

        if "gamma" in config:
            self.gamma = config["gamma"]
        else:
            self.gamma = 1.4

        if "Ca" in config:
            self.Ca = config["Ca"]
        else:
            self.Ca = 1.

        if "Re_inv" in config:
            self.Re_inv = config["Re_inv"]
            if self.Re_inv <= 0.:
                raise ValueError(self.Re_inv)
            else:
                self.viscosity = True
        else:
            self.viscosity = False
            self.Re_inv = 0

        if "Web" in config:
            self.Web = config["Web"]
            if self.Web <= 0.:
                raise ValueError(self.Web)
            else:
                self.tension = True
        else:
            self.tension = False
            self.Web = 0

        if self.model == "RPE":
            self.num_RV_dim = 2
            self.state = np.array([self.R,self.V])
        else:
            raise NotImplementedError

    
    def pbw(self):
        self.cpbw = self.Ca*((self.R0/self.R)**(3.0*self.gamma)) - self.Ca + 1.
        if self.tension:
            self.cpbw += 2./(self.Web*self.R0)*(self.R0/self.R)**(3.*self.gamma)

    def rpe(self,p):
        self.pbw()
        rhs = -1.5*self.V**2.0 + (self.cpbw - p)/self.R
        if self.viscosity:
            rhs -= 4.0*self.Re_inv*self.V/(self.R**2.0)
        return [self.V, rhs]


    def rhs(self,p):
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

    def wrap(self,t,y):
        self.R = y[0]
        self.V = y[1]
        return np.array(self.rpe(1.1))
        
    def mc_solve(self):
        t0 = 0.
        T = 10.
        y0 = np.array([1.,0.])
        ret = sp.solve_ivp(self.wrap,(t0,T),y0,method='RK45')
        # print(ret.t)
        # print(ret.y)
        plt.plot(ret.t,ret.y[0])
        plt.show()

if __name__ == "__main__":

    config = {}
    config["model"] = "RPE"
    config["R"] = 1.
    config["V"] = 0.
    config["gamma"] = 1.4
    config["Ca"] = 1.
    # config["Re_inv"] = 0
    # config["Web"] = 0

    mybub = bubble_model(config=config,R0=1.)
    # rhs = mybub.rhs(p=1.1)
    # print('rhs = ',rhs)

    mybub.mc_solve()

