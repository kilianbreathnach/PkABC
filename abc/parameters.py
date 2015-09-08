'''

Parameters

'''
from scipy.stats import uniform
from scipy.stats import norm 

class Params(object): 

    def __init__(self, 
            prior_dict
            ): 
        """ Class that describes parameters
        """
        
        self.prior_dict = prior_dict.copy() # dictionary that specifies prior info
        self.prior_shape = self.prior_dict['shape']


    def prior(self): 
        """ Prior from which parameters are drawn
        """

        if self.prior_shape == 'uniform': 

            loc = self.prior_dict['min']
            scale = self.prior_dict['max'] - self.prior_dict['min']

            priorz = uniform( loc , scale ) 

        elif self.prior_shape == 'gauss': 

            loc = self.prior_dict['mean']
            scale = self.prior_dict['stddev']

            priorz = norm( loc , scale ) 

        return priorz

    def simulator(self, theta): 
        """ Simulator 
        """
       
        loc = theta

        simz = norm(loc, 1.0)

        return simz.pdf


if __name__=="__main__": 
    parz = Params(
            {'shape': 'gauss', 'mean': 0.0, 'stddev': 1.0}
            )
    prz = parz.prior()
    print prz.rvs(size=100)
    print prz.pdf(0.0)
