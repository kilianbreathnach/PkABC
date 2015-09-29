"""

Simulator for PMC_ABC 


"""

from halotools.empirical_models import Zheng07


class Simul(object): 

    def __init__(self): 

        self.model = Zheng07()
    

    def nz(self, theta_star): 

        self.model.param_dict['logMmin'] = theta_star[1]
        self.model.param_dict['sigma_logM'] = theta_star[0]
        self.model.param_dict['alpha'] = theta_star[2]

        self.model.param_dict['logM0'] = theta_star[3]
        self.model.param_dict['logM1'] = theta_star[4]
	#'alpha': 1.06,
 	#'logM0': 11.38,
 	#'logM1': 13.31,
 	#'logMmin': 12.02,
 	#'sigma_logM': 0.26} 
        self.model.populate_mock()
	
        self.nz = self.model.mock.number_density
        
        return self.nz
