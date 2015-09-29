"""

Simulator for PMC_ABC 


"""

from halotools.empirical_models import Zheng07

class Simul(object): 

    def __init__(self): 

        self.model = Zheng07()
    

    def sum_stat(self, theta_star): 
        """ 
        theta_star = [alpha, logM0, logM1, logMmin, sigma_logM]

        Zheng07 parameters = {
        'alpha': 1.06,
        'logM0': 11.38,
        'logM1': 13.31,
        'logMmin': 12.02,
        'sigma_logM': 0.26
        }
        nbar(z) = 0.0047808
        xi[1] = 6.76599366e+02 

        """
        #self.model.param_dict['alpha'] = theta_star[1]
        #self.model.param_dict['logM0'] = theta_star[0]
        #self.model.param_dict['logM1'] = theta_star[3]
        self.model.param_dict['logMmin'] = theta_star[0]
        self.model.param_dict['sigma_logM'] = theta_star[1] 

        self.model.populate_mock()

        r, xi = self.model.compute_galaxy_clustering()
    
        # number density of mock catalog
        nz = self.model.mock.number_density
        # magnitude of the correlation function 
        corr_mag = xi[1]

        return [nz, corr_mag]
