import matplotlib.pyplot as plt
from time_advancer import time_advancer
from mc import mc
from moments import get_moments, write_moments
import warnings
from sys import exit
from plotting import plot_moments, plot_integrands, plot_integrands_sub
import numpy as np

warnings.simplefilter("error")


def inputs():

    # Compute natural frequency of the bubbles
    Ro = 1
    omega = np.sqrt(3*1.4/(Ro**2.))
    T = 2.*np.pi/omega

    config = {}
    config["advancer"] = {}
    config["wave"] = {}
    config["pop"] = {}
    config["model"] = {}
    config["mc"] = {}

    # Advancer parameters
    # config["advancer"]["method"] = "RK3"
    config["advancer"]["method"] = "RK23"
    config["advancer"]["dt"] = 1.0e-4
    config["advancer"]["T"] = 10*T
    config["advancer"]["error_tol"] = 1.0e-4
    config["advancer"]["io"] = 1

    # Acoustic
    # config["wave"]["amplitude"] = 3
    # config["wave"]["amplitude"] = 1.00001
    config["wave"]["form"] = "sine"
    # config["wave"]["form"] = "constant"
    config["wave"]["period"] = 1.2*T
    config["wave"]["cycles"] = 1000.

    # Monte Carlo
    config["mc"]["Nsamples"] = 1000
    config["mc"]["Ntimes"] = 100

    # Population properties
    config["pop"]["NR0"] = 301
    config["pop"]["shape"] = "lognormal"
    config["pop"]["binning"] = "Simpson"
    # config["pop"]["binning"] = "GH"
    config["pop"]["muR0"] = 1.
    config["pop"]["sigR0"] = 0.3
    config["pop"]["moments"] = [[1, 0], [0, 1], [1, 1]]
    # config["pop"]["moments"] = [ [3, 2], [2, 1], [3, 0], [ 3*(1-1.4), 0, 3*1.4 ] ]
    # config["pop"]["Nfilt"] = 100
    # config["pop"]["Tfilt"] = 1


    # Bubble properties
    config["model"]["model"] = "RPE"
    # config["model"]["model"] = "KM"
    # config["model"]["model"] = "Linear"
    # config["model"]["R"] = 1.0
    config["model"]["V"] = 0.0
    config["model"]["gamma"] = 1.4
    # config["model"]["c"] = 1000.
    # config["model"]["Ca"] = 0.97
    # config["model"]["Re_inv"] = 1/1000.
    # config["model"]["Web"] = 13.9

    return config


def advance_classes(config, sols):

    myadv = time_advancer(config=config["advancer"])
    myadv.initialize_state(pop_config=config["pop"], model_config=config["model"])
    myadv.initialize_wave(wave_config=config["wave"])
    myadv.run()
    sols.append(myadv)

    return sols


def advance_mc(config):

    mymc = mc(config)
    mymc.run()


if __name__ == "__main__":

    config = inputs()
    sols = []

    pratios = [ 0.1 ]
    for i, p in enumerate(pratios):
        config["wave"]["amplitude"] = p
        # Name of HDF5 file
        config["advancer"]["file_name"] = "out/" + config["wave"]["form"] + "-p" + str(p)
        advance_classes(config, sols)
        sols = get_moments(sols)
        write_moments(sols[i])
        # plot_moments(sols)

    # plot_integrands(sols)
    # plot_integrands_sub(sols)

    # advance_mc(config)

    # plt.savefig('out.pdf')


    # config["pop"]["NR0"] = 150
    # config["pop"]["Nfilt"] = 0
    # advance_classes(config, sols)

    # print(sols[0].state.R0)

    # config["pop"]["NR0"] = 51
    # config["pop"]["Nfilt"] = 300
    # advance_classes(config, sols)


    # config["pop"]["NR0"] = 51
    # config["pop"]["Nfilt"] = 600
    # advance_classes(config, sols)


    # config["pop"]["NR0"] = 501
    # config["pop"]["Nfilt"] = 0
    # advance_classes(config, sols)
