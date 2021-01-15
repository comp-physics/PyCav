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
    config["advancer"]["T"] = 10
    config["advancer"]["p"] = 3

    # Population properties
    config["pop"]["NR0"] = 1
    config["pop"]["shape"] = "lognormal"
    config["pop"]["binning"] = "Simpson"
    config["pop"]["muR0"] = 1.0
    config["pop"]["sigR0"] = 0.3

    # Bubble properties
    config["model"]["model"] = "RPE"
    config["model"]["R"] = 1.0
    config["model"]["V"] = 0.0
    config["model"]["gamma"] = 1.4
    # config["model"]["Ca"] = 0.5
    config["model"]["Re_inv"] = 0.001
    # config["model"]["Web"] = 10.
    return config


def advance_classes(config):

    myadv = adv.time_advancer(config=config["advancer"])
    myadv.initialize_state(pop_config=config["pop"], model_config=config["model"])
    myadv.run()

    return


def advance_mc(config):

    mymc = mc(config)
    mymc.simulate_sample(N=100)

    # samp = mymc.random_sample(N=100)

    # state = bubble_state()
    # val = state.vals
    # p = 1.1
    # dt = 0.1

    # print(samp)

    return


if __name__ == "__main__":

    config = inputs()

    # Classes
    # advance_classes(config)

    # MC
    advance_mc(config)
