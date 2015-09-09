'''
    
Test simulator for ABC

'''
from scipy.stats import norm 

def simz(theta): 
    """ Given parameter vector theta, return parameterized method 
    """
       
    loc = theta[0]
    scale = theta[1]

    simz = norm(loc, 10.0**scale)

    return simz.pdf
