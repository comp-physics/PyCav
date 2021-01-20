import bubble_state as bs
import waveforms as wf
import numpy as np
from sys import exit
import matplotlib.pyplot as plt


class time_advancer:
    def __init__(self, config={}):
        if "T" in config:
            self.T = config["T"]
            if self.T <= 0:
                raise ValueError(self.T)
        else:
            raise Exception("No final time T")

        if "dt" in config:
            self.dt = config["dt"]
            if self.dt <= 0.0:
                raise ValueError(self.dt)
        else:
            raise Exception("No time step dt")

        if "method" in config:
            self.method = config["method"]
        else:
            self.method = "Euler"

        if self.method == "Euler":
            self.advance = self.euler
            self.n_stages = 0
        elif self.method == "RK2":
            self.advance = self.rk2
            self.n_stages = 1
        elif self.method == "RK3":
            self.advance = self.rk3
            self.n_stages = 1
        elif self.method == "RK23":
            self.advance = self.rk23
            self.n_stages = 1
        else:
            raise NotImplementedError

        return

    def initialize_state(self, pop_config=None, model_config=None):
        self.state = bs.bubble_state(pop_config=pop_config, model_config=model_config)

    def initialize_wave(self, wave_config=None):
        self.wave = wf.waveforms(config=wave_config)

    def euler(self):
        # print(self.state.vals[0,:])
        # vals = self.state.vals.copy()
        f0 = self.state.vals.copy()

        p = self.wave.p(self.time)
        l1 = self.state.get_rhs(f0, p)

        self.state.vals[:, :] = f0 + self.dt * l1

        # update = self.state.vals.copy() + self.dt*l1.copy()
        # self.state.vals = update

        # self.state.vals += self.dt * l1
        # print(self.state.vals[0,:])

    def rk2(self):
        f0 = self.state.vals.copy()

        p = self.wave.p(self.time)
        l1 = self.state.get_rhs(f0, p)
        f1 = f0 + self.dt * l1

        pdt = self.wave.p(self.time + self.dt)
        L = self.state.get_rhs(f1, pdt)

        self.state.vals[:, :] = 0.5 * f0 + 0.5 * (f1 + self.dt * L)

    def rk3(self):
        f0 = self.state.vals.copy()

        p = self.wave.p(self.time)
        l1 = self.state.get_rhs(f0, p)
        f1 = f0 + self.dt * l1

        pdt = self.wave.p(self.time + self.dt)
        L = self.state.get_rhs(f1, pdt)
        f2 = 0.75 * f0 + 0.25 * (f1 + self.dt * L)

        pdt2 = self.wave.p(self.time + self.dt / 2.0)
        L2 = self.state.get_rhs(f2, pdt2)

        self.state.vals[:, :] = 1.0 / 3.0 * f0 + 2.0 / 3.0 * (f2 + self.dt * L2)

    def rk23(self):

# RK23[mom_,myrhs_,t_,dt_]:=Module[{moms=mom,mome,momstemp1,momstemp2,rhs},
# 	(* SSP-RK2 *)
# 	{moms,rhs}=myrhs[moms,t];
# 	momstemp1=moms+dt rhs;
# 	{momstemp1,rhs}=myrhs[momstemp1,t+dt];
# 	mome=(1/2) moms+(1/2)(momstemp1+dt rhs);

# 	(* SSP-RK3 *)
# 	momstemp2=(3/4)moms+(1/4)(momstemp1+dt rhs);
# 	{momstemp2,rhs}=myrhs[momstemp2,t+dt/2];
# 	moms=(1/3)moms+(2/3)(momstemp2+dt rhs);
# 	Return[{moms,err[moms,mome]},Module];
# ];

        # SSP-RK2
        f0 = self.state.vals.copy()

        p = self.wave.p(self.time)
        l1 = self.state.get_rhs(f0, p)
        f1 = f0 + self.dt * l1

        pdt = self.wave.p(self.time + self.dt)
        L = self.state.get_rhs(f1, pdt)

        mome = 0.5 * f0 + 0.5 * (f1 + self.dt * L)

        # SSP-RK3
        f2 = 0.75 * f0 + 0.25 * (f1 + self.dt * L)
        pdt2 = self.wave.p(self.time + self.dt / 2.0)
        L2 = self.state.get_rhs(f2, pdt2)
        mom = 1.0 / 3.0 * f0 + 2.0 / 3.0 * (f2 + self.dt * L2)
        err = np.norm(mom-mome)

        self.state.vals[:, :] = mom

    def run(self):
        self.time = 0.0
        i_step = 0
        step = True
        self.save = []
        self.times = []
        self.moms = []
        np.set_printoptions(precision=24)

        while step:
            print("step = ", i_step)
            self.times.append(self.time)
            self.save.append(self.state.vals.copy())
            self.moms.append(self.state.get_quad())
            self.advance()
            i_step += 1
            self.time += self.dt

            if self.time >= self.T:
                step = False

        self.save = np.array(self.save, dtype=np.float32)
        self.moms = np.array(self.moms, dtype=np.float32)
        self.plot()

    def plot(self):
        # plot R evolution for all quad points
        # for i in range(self.state.NR0):
        #     plt.plot(self.times, self.save[:,i,0])
        # plt.plot(self.times, self.save[:,0,0])
        # plt.xlabel("$t$")
        # plt.ylabel("$R(t)$")
        # plt.show()

        fig, ax = plt.subplots(1, self.state.Nmom)
        for i in range(self.state.Nmom):
            ax[i].plot(self.times, self.moms[:, i])
            ax[i].set(xlabel="$t$", ylabel="$M$" + str(self.state.moments[i]))
        plt.tight_layout()


if __name__ == "__main__":
    config = {}
    config["method"] = "Euler"
    config["dt"] = 1.0e-5
    config["T"] = 30.0
    config["p"] = 1.0
    myadv = time_advancer(config=config)
