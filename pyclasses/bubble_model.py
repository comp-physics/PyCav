class bubble_model:

    def __init__(self,
        model="RPE",
        R=1.,
        V=0.,
        R0=1.,
        gamma=1.4,
        Ca=1.,
        Re_inv=0,
        Web=0
        ):

        self.model = model
        self.Re_inv = Re_inv
        self.Web = Web
        self.gamma = gamma
        self.Ca = Ca
        self.R0 = R0
        self.R = R
        self.V = V

        self.viscosity = False
        self.tension = False
        if self.Re_inv > 0.:
            self.viscosity == True
        if self.Web > 0.:
            self.tension == True

        if self.model == "RPE":
            self.num_RV_dim = 2
        else:
            raise NotImplementedError

    
    def pbw(self):
        self.cpbw = self.Ca*((self.R0/self.R)**(3.0*self.gamma)) - self.Ca + 1.
        if self.tension:
            self.cpbw += 2./(self.Web*self.R0)*(self.R0/self.R)**(3.*gamma)

    def rpe(self,p):
        rhs = -1.5*self.V**2.0 + (self.cpbw - p)/self.R
        if self.viscosity:
            rhs -= 4.0*self.Re_inv*self.V/(self.R**2.0)
        return rhs

    def rhs(self,p):
        self.pbw()
        if self.model == "RPE":
            rhs = self.rpe(p)
        else:
            raise NotImplementedError
        return rhs

if __name__ == "__main__":

    mybub = bubble_model()
    rhs = mybub.rhs(p=1.0)
    print('rhs = ',rhs)
