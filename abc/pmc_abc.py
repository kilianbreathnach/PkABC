'''

PMC-ABC

'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import uniform
from scipy.stats import norm
import time 


from distance import test_dist
from parameters import Params

def pmc_abc(data, 
        prior_shape = 'uniform', 
        N = 1000,
        eps0 = 0.01, 
        T = 20
        ): 

    start_time = time.time()


    toolz = Params({
                'shape': prior_shape, 
                'min': -1.0, 
                'max': 1.0
                })

    t = 0 
    theta_t, w_t = [] , [] 

    for i in range(N): 
        
        theta = toolz.prior().rvs(size=1)[0]
        model = toolz.simulator( theta )
        rho = test_dist(data[1], model(data[0]))

        while rho > eps0: 

            theta = toolz.prior().rvs(size=1)[0]
            rho = test_dist(data[1], model(data[0]))

        theta_t.append(theta)
        w_t.append(1.0/np.float(N))

    theta_t = np.array(theta_t)
    w_t = np.array(w_t)
    
    sig_t = 2.0 * np.var(theta_t)
    
    # write out 
    np.savetxt(
            ''.join(['theta_w_t', str(t), '.dat']), 
            np.c_[ theta_t, w_t ], 
            fmt=['%10.5f', '%10.5f']
            )

    print 'Initial Pool ', time.time() - start_time

    fig = plt.figure(1)
    sub = fig.add_subplot(111)

    theta_dist, theta_binedge = np.histogram(theta_t, bins=20, weights=w_t)
    theta_low = theta_binedge[:-1]
    theta_high = theta_binedge[1:]
    theta_mid = [ 0.5 * (theta_low[it] + theta_high[it]) for it in range(len(theta_low))]

    sub.plot(theta_mid, theta_dist)
    sub.set_xlim([-2.0, 2.0])

    plt.show()
    plt.close()


    start_time = time.time()
    while t < T: 

        eps_t = np.percentile(rho, 75)
        print 'Epsilon t', eps_t

        theta_t_1 = theta_t.copy()
        w_t_1 = w_t.copy()
        sig_t_1 = sig_t.copy()

        for i in range(N): 
            start_time = time.time()

            theta_star = weighted_sampling( theta_t_1, w_t_1 ) 
            theta_starstar = norm(theta_star, np.sqrt(sig_t_1)).rvs(size=1)[0]
            model_starstar = toolz.simulator( theta_starstar )
            rho = test_dist(data[1], model_starstar(data[0])) 
        
            while rho > eps_t: 

                theta_star = weighted_sampling( theta_t_1, w_t_1 ) 
                theta_starstar = norm(theta_star, np.sqrt(sig_t_1)).rvs(size=1)[0]
                model_starstar = toolz.simulator( theta_starstar )
                rho = test_dist(data[1], model_starstar(data[0])) 

            #print theta_star, theta_starstar
            theta_t[i] = theta_starstar
            w_t[i] = toolz.prior().pdf(theta_t[i]) / np.sum(w_t_1 * norm(theta_t[i], np.sqrt(sig_t_1)).pdf(theta_t_1))
            #print test_dist(data[1], model_starstar(data[0])), w_t[i]
            
            #print 'For loop ', time.time() - start_time
        
        fig = plt.figure(1)
        sub = fig.add_subplot(111)
    
        theta_dist, theta_binedge = np.histogram(theta_t, bins=20, weights=w_t)
        theta_low = theta_binedge[:-1]
        theta_high = theta_binedge[1:]
        theta_mid = [ 0.5 * (theta_low[it] + theta_high[it]) for it in range(len(theta_low))]

        sub.plot(theta_mid, theta_dist)
        sub.set_xlim([-2.0, 2.0])

        plt.show()
        plt.close()

        sig_t = 2.0 * np.var(theta_t)
        np.savetxt(
                ''.join(['theta_w_t', str(t), '.dat']), 
                np.c_[ theta_t, w_t ], 
                fmt=['%10.5f', '%10.5f']
                )

        t += 1 
        print t, ' ;D'

def weighted_sampling(theta, w): 
    """ Weighted sampling
    """
    
    w_cdf = w.cumsum()/w.sum()

    rand1 = w.sum() * np.random.random(1)
    #print w.sum(), rand1

    cdf_closest_index = min(range(len(w_cdf)), key = lambda i: abs(w_cdf[i]-rand1[0]))
    closest_theta = theta[cdf_closest_index]

    return closest_theta 


if __name__=='__main__': 
    # fake data
    data_x = uniform( -1.0, 2.0).rvs(size=1000)
    data_y = norm(0.0, 1.0).pdf(data_x)
    data = [data_x, data_y]

    #fig = plt.figure(1)
    #sub = fig.add_subplot(111)

    #sub.scatter(data_x, data_y)
    
    pmc_abc(data, 
        prior_shape = 'uniform', 
        N = 1000,
        eps0 = 10.0, 
        T=20
        )
    
