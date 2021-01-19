import matplotlib.pyplot as plt
import time_advancer as adv
from mc import mc

def inputs():
    config = {}
    config["advancer"] = {}
    config["pop"] = {}
    config["model"] = {}

    # Advancer parameters
    config["advancer"]["method"] = "Euler"
    config["advancer"]["dt"] = 0.001
    config["advancer"]["T"] = 1
    config["advancer"]["p"] = 3

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
    myadv.run()


def advance_mc(config):

    mymc = mc(config)
    sols = mymc.simulate_sample(Nmc=200, Nt=100)
    moments = mymc.state.moment(sols)

    fig, ax = plt.subplots(1,mymc.state.Nmom)
    for i in range(mymc.state.Nmom):
        ax[i].plot(sols[i].t, moments[i])
        ax[i].set(
                xlabel="$t$",
                ylabel="$M$" + str(mymc.state.moments[i])
                )


if __name__ == "__main__":

    config = inputs()

    # Classes
    advance_classes(config)

    # MC
    advance_mc(config)

    plt.show()

    # plt.savefig('out.pdf')
