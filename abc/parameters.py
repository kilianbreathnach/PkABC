'''

Parameters

'''
from scipy.stats import uniform
from scipy.stats import norm 
from scipy.stats import multivariate_normal
import numpy as np

class Params(object): 

    def __init__(self, 
            prior_dict
            ): 
        """ Class that describes parameters
        """
        
        self.prior_dict = prior_dict.copy() # dictionary that specifies prior info

    def prior(self): 
        """ Prior from which parameters are drawn
        """

        priorz = []  

        for key in self.prior_dict.keys(): 

            prior_key = self.prior_dict[key]

            if prior_key['shape'] == 'uniform': 
                
                loc = prior_key['min']

                scale = prior_key['max'] - prior_key['min']

                priorz.append( uniform( loc , scale ) ) 

            elif prior_key['shape'] == 'gauss': 

                loc = prior_key['mean']

                scale = prior_key['stddev']

                priorz.append( norm( loc , scale ) )

        return priorz

    def simulator(self, theta): 
        """ Simulator 
        """
       
        loc = theta[0]
        scale = theta[1]

        simz = norm(loc, scale)

        return simz.pdf


if __name__=="__main__": 
    parz = Params({
                'mu': {'shape': 'uniform', 'min': -1.0, 'max': 1.0}, 
                'sigma': { 'shape': 'gauss', 'mean': 0.0, 'stddev': 1.0}
                })
    prz_mu, prz_sigma = parz.prior()
    print prz_mu.rvs(size=100)
    print prz_mu.pdf(0.0)
    print prz_sigma.rvs(size=100)
    print prz_sigma.pdf(0.0)
    w = parz.simulator([0, 1])
    print w(0)
    #theta_t = np.zeros((2, N))


