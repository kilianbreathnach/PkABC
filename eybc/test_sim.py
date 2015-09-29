'''
    
Test simulator for ABC

'''
from scipy.stats import norm 
from simulator import Simul

def simz(theta): 
    """ Given parameter vector theta, return parameterized method 
    """
       
    loc = theta[0]
    scale = theta[1]

    simz = norm(loc, scale)

    return simz.pdf
    
def halotools_nz_sim(sigma, m_min, model): 
    ''' Number density simulator from halotools

    Parameters
    ----------
    param : dictionary of parameters 
    model : halotools Zheng07 model object
    '''
    
    model.param_dict['logMmin'] = m_min
    model.param_dict['sigma_logM'] = sigma
    model.populate_mock()

    return model.mock.number_density
