import scipy.stats as stats
import scipy.special as sc
import math
import numpy as np

def raw_gaussian_moments_univar(num_moments, mu, sigma):
    """
    This function returns raw 1D-Gaussian moments as a function of
    mean (mu) and standard deviation (sigma)
    """
    moments = np.zeros(num_moments)
    for i_moment in range(num_moments):
        moments[i_moment] = stats.norm.moment(i_moment, mu, sigma)
    return moments

