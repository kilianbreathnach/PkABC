"""

Simulator for PMC_ABC 


"""

from halotools.empirical_models import Zheng07

class Simul(object): 

    def __init__(self): 

        self.model = Zheng07()
    

    def nz(self, sigma, m_min): 

        self.model.param_dict['logMmin'] = m_min
        self.model.param_dict['sigma_logM'] = sigma
        self.model.populate_mock()

        self.nz = self.model.mock.number_density

        return self.nz
