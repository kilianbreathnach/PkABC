import numpy as np 
from scipy.stats import uniform
from scipy.stats import norm

from pmc_abc import PmcAbc
from test_sim import simz

if __name__=="__main__":
    # fake data
    data_x = uniform( -1.0, 2.0).rvs(size=1000)
    data_y = norm(0.0, 1.0).pdf(data_x)
    data = {'input': data_x, 'output': data_y}

    pmcabc_test = PmcAbc(data, simz, 
            prior_dict = {
                'mu': {'shape': 'uniform', 'min': -0.5, 'max': 0.5}, 
                'sigma': {'shape': 'uniform', 'min': 0.5, 'max': 1.5}
                }, 
            N=100, eps0 = 0.5, T = 10, Nthreads=3)

    pmcabc_test.pmc_abc()
