'''
Distance computation for ABC
'''

import numpy as np
import warnings

def test_dist( 
        data, 
        model, 
        disttype = 'halotool_nz'
        ): 
    ''' Calculate distance between data set D and D_s from model 
    Parameters:
    -----------
    data : dictionary with keys 'input' and 'output' 
    model : method parameterized with theta 
    '''

    if disttype == 'dumb':
        data_out = data['output']
        model_out = model(data['input'])
        rho_mean = np.abs( 
                (np.mean(data_out) - np.mean(model_out))/np.mean(data_out) 
                ) 
        rho_sigma =  np.abs( 
                (np.std(data_out) - np.std(model_out))/np.std(data_out) 
                )

        dist = rho_mean + rho_sigma

    elif disttype == 'halotool_nz': 
        data_out = data['output']
        dist = np.abs(data_out[0] - model[0])

    elif disttype == 'halotool_sum_stat': 
        data_out = data['output']
        dist_nz = np.abs(data_out[0] - model[0])/model[0]
        dist_corr = np.abs(data_out[1] - model[1])/model[1] 
        dist = dist_nz + dist_corr

    return dist
