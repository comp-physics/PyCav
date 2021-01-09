import stats_util as stats
import bubble_model as bub
import numpy as np

class time_advancer:

    def __init__(self, 
        method="Euler",
        time_step=0.,
        final_time=0.,
        NR0=1,
        bub_config=None
        ):

        self.method = method
        self.time_step = time_step
        self.final_time = final_time
        self.num_R0_nodes = NR0

        self.bubble = bub.bubble_model()

        self.num_RV_dim = self.bubble.num_RV_dim
        self.num_dim = self.num_R0_nodes + self.num_RV_dim

        self.state = np.zeros(self.num_dim)
        self.rhs = np.zeros(self.num_dim)

        if self.method == "Euler":
            self.advance = self.euler
            self.n_stages = 1
        else:
            raise NotImplementedError

        return

    def initialize_state_gaussian_univar(self, mu, sigma):

        self.R0state = stats.raw_gaussian_moments_univar(self.num_dim, mu, sigma)
        print(self.R0state)

        return

    def euler(self):

        self.stage_state[0] = self.state.copy()
        self.bubble_mgr.compute_rhs(self.stage_state[0], self.stage_k[0])
        self.state = self.stage_state[0] + self.time_step * self.stage_k[0]

        return

    def run(self):

        self.time = 0.0

        i_step = 0
        step = True
        while step == True:

            self.advance()

            i_step += 1
            self.time += self.time_step

            if self.time > self.final_time or i_step >= self.num_steps:
                step = False

        return

if __name__ == "__main__":

    adv = time_advancer()

    rhs = adv.bubble.rhs(p=2.0)
    print('rhs = ',rhs)
