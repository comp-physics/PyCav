import matplotlib.pyplot as plt
import time_advancer as adv
from mc import mc


def inputs():
    config = {}
    config["advancer"] = {}
    config["wave"] = {}
    config["pop"] = {}
    config["model"] = {}
    config["mc"] = {}

    # Advancer parameters
    # config["advancer"]["method"] = "Euler"
    config["advancer"]["method"] = "RK12"
    config["advancer"]["dt"] = 1.0e-5
    config["advancer"]["T"] = 10
    config["advancer"]["error_tol"] = 1.0e-3

    # Acoustic
    config["wave"]["amplitude"] = 3
    config["wave"]["form"] = "sine"
    # config["wave"]["form"] = "constant"
    config["wave"]["period"] = 2.0
    config["wave"]["cycles"] = 2.0

    # Monte Carlo
    config["mc"]["Nsamples"] = 1000
    config["mc"]["Ntimes"] = 100

    # Population properties
    config["pop"]["NR0"] = 201
    config["pop"]["shape"] = "lognormal"
    config["pop"]["binning"] = "Simpson"
    config["pop"]["muR0"] = 1.0
    config["pop"]["sigR0"] = 0.1
    config["pop"]["moments"] = [[0, 0], [1, 0], [0, 1]]

    # Bubble properties
    config["model"]["model"] = "RPE"
    # config["model"]["R"] = 1.0
    config["model"]["V"] = 0.0
    config["model"]["gamma"] = 1.4
    # config["model"]["Ca"] = 0.5
    # config["model"]["Re_inv"] = 0.001
    # config["model"]["Web"] = 10.

    return config


def advance_classes(config):

    myadv = adv.time_advancer(config=config["advancer"])
    myadv.initialize_state(pop_config=config["pop"], model_config=config["model"])
    myadv.initialize_wave(wave_config=config["wave"])
    myadv.run()


def advance_mc(config):

    mymc = mc(config)
    sols = mymc.simulate_sample()
    moments = mymc.state.moment(sols)

    fig, ax = plt.subplots(1, mymc.state.Nmom)
    for i in range(mymc.state.Nmom):
        ax[i].plot(sols[i].t, moments[i])
        ax[i].set(xlabel="$t$", ylabel="$M$" + str(mymc.state.moments[i]))


if __name__ == "__main__":

    config = inputs()

    # Classes
    advance_classes(config)

    # MC
    advance_mc(config)

    plt.show()

    # plt.savefig('out.pdf')
