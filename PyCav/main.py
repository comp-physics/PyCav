import time_advancer as adv
from mc import mc
import matplotlib.pyplot as plt

def inputs():
    config = {}
    config["advancer"] = {}
    config["pop"] = {}
    config["model"] = {}

    # Advancer parameters
    config["advancer"]["method"] = "Euler"
    config["advancer"]["dt"] = 0.001
    config["advancer"]["T"] = 10
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
    sols = mymc.simulate_sample(Nmc=10, Nt=100)
    # for i in range(10):
    #     print(sols[i].y[0])

    moments = mymc.state.moment(sols)
    for i in range(mymc.state.Nmom):
        plt.subplot(1, mymc.state.Nmom, i + 1)
        plt.plot(sols[i].t, moments[i])
        plt.xlabel("$t$")
        plt.ylabel("$M$" + str(mymc.state.moments[i]))

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":

    config = inputs()

    # Classes
    advance_classes(config)

    # MC
    # advance_mc(config)
