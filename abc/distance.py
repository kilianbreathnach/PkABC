'''

Distance computation for ABC

'''

import numpy as np

def test_dist(
        data, 
        model
        ): 
    ''' simple distance
    '''

    rho_mean = np.abs( (np.mean(data) - np.mean(model))/np.mean(data) ) 
    rho_sigma =  np.abs( (np.std(data) - np.std(model))/np.std(data) ) 
    
    dist = rho_mean + rho_sigma

    return dist
