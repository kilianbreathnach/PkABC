"""

Simulator for PMC_ABC 


"""

from halotools.empirical_models import Zheng07


class Simul(object): 

    def __init__(self): 

        self.model = Zheng07()
    

    def nz(self, theta_star): 

        self.model.param_dict['logMmin'] = theta_star[0]
        self.model.param_dict['sigma_logM'] = theta_star[1] 
        self.model.populate_mock()

        self.nz = self.model.mock.number_density

        return self.nz
