import matplotlib.pyplot as plt
from time_advancer import time_advancer
from mc import mc
import warnings
warnings.simplefilter('error')

def inputs():
    config = {}
    config["advancer"] = {}
    config["wave"] = {}
    config["pop"] = {}
    config["model"] = {}
    config["mc"] = {}

    # Advancer parameters
    # config["advancer"]["method"] = "Euler"
    config["advancer"]["method"] = "RK23"
    config["advancer"]["dt"] = 1.0e-5
    config["advancer"]["T"] = 40
    config["advancer"]["error_tol"] = 1.0e-3
    config["advancer"]["Nfilt"] = 1

    # Acoustic
    config["wave"]["amplitude"] = 3
    config["wave"]["form"] = "sine"
    # config["wave"]["form"] = "constant"
    # config["wave"]["period"] = 4.0
    # config["wave"]["cycles"] = 2.0

    # Monte Carlo
    config["mc"]["Nsamples"] = 1000
    config["mc"]["Ntimes"] = 100

    # Population properties
    config["pop"]["NR0"] = 501
    config["pop"]["shape"] = "lognormal"
    config["pop"]["binning"] = "Simpson"
    # config["pop"]["binning"] = "GH"
    config["pop"]["muR0"] = 1.0
    config["pop"]["sigR0"] = 0.3
    config["pop"]["moments"] = [ [1, 0], [0, 1], [1, 1] ]
    # config["pop"]["moments"] = [ [3, 2], [2, 1], [3, 0], [ 3*(1-1.4), 0, 3*1.4 ] ]


    # Bubble properties
    config["model"]["model"] = "RPE"
    # config["model"]["model"] = "KM"
    # config["model"]["model"] = "Linear"
    # config["model"]["R"] = 1.0
    config["model"]["V"] = 0.0
    config["model"]["gamma"] = 1.4
    # config["model"]["c"] = 1000.
    # config["model"]["Ca"] = 0.5
    # config["model"]["Re_inv"] = 1/100.
    # config["model"]["Web"] = 13.9

    return config


def advance_classes(config):

    myadv = time_advancer(config=config["advancer"])
    myadv.initialize_state(pop_config=config["pop"], model_config=config["model"])
    myadv.initialize_wave(wave_config=config["wave"])
    myadv.run()


def advance_mc(config):

    mymc = mc(config)
    mymc.run()


if __name__ == "__main__":

    config = inputs()

    config["pop"]["NR0"] = 501
    config["advancer"]["Nfilt"] = 0
    advance_classes(config)
    config["pop"]["NR0"] = 101
    config["advancer"]["Nfilt"] = 0
    advance_classes(config)
    config["pop"]["NR0"] = 101
    config["advancer"]["Nfilt"] = 40
    advance_classes(config)
    # advance_mc(config)

    plt.show()
    # plt.savefig('out.pdf')
