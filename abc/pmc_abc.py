'''

PMC-ABC

'''
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import uniform
from scipy.stats import norm
from scipy.stats import multivariate_normal
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
                'mu': {'shape': 'gauss', 'mean': .0, 'stddev': 1.0}, 
                'sigma': { 'shape': 'uniform', 'min': 0.0, 'max': 2.0}
                })

    t = 0 
    theta_t = np.zeros((2, N))
    w_t = np.zeros((N))

    for i in xrange(N): 
        
        theta_star = np.zeros(2)
        theta_star[0] = toolz.prior()[0].rvs(size=1)[0]
        theta_star[1] = toolz.prior()[1].rvs(size=1)[0]
	#print theta_star
        model = toolz.simulator( theta_star )
        

        rho = test_dist(data[1], model(data[0]))

        while rho > eps0: 

            theta_star[0] = toolz.prior()[0].rvs(size=1)[0]
            theta_star[1] = toolz.prior()[1].rvs(size=1)[0]

            model = toolz.simulator( theta_star )

            rho = test_dist(data[1], model(data[0]))
        
        theta_t[:,i] = theta_star

        w_t[i] = 1.0/np.float(N)

    sig_t = 2.0 * np.cov( theta_t )    # covariance matrix
    print sig_t
    
    # write out 
    np.savetxt(
            ''.join(['theta_w_t', str(t), '.dat']), 
            np.c_[theta_t[0,:], theta_t[1,:], w_t ]
            )

    print 'Initial Pool ', time.time() - start_time

    fig = plt.figure(1)
    sub = fig.add_subplot(111)
    sub.hist2d(theta_t[0,:], theta_t[1,:], weights=w_t, bins=100, range=[[-2.0, 2.0],[0.0, 2.0]])
    plt.savefig("theta0.png")
    plt.close()

    start_time = time.time()

    while t < T: 

        eps_t = np.percentile(rho, 75)
        print 'Epsilon t', eps_t

        theta_t_1 = theta_t.copy()
        w_t_1 = w_t.copy()
        sig_t_1 = sig_t.copy()

        for i in xrange(N): 
            start_time = time.time()

            theta_star = weighted_sampling( theta_t_1, w_t_1 ) 
            theta_starstar = multivariate_normal( theta_star, sig_t_1 ).rvs(size=1)
	    
	    while theta_starstar[1] < 0 :
		    theta_star = weighted_sampling( theta_t_1, w_t_1 )
		    theta_starstar = multivariate_normal(theta_star, sig_t_1).rvs(size=1)
	    #print theta_starstar
            model_starstar = toolz.simulator( theta_starstar )
            rho = test_dist(data[1], model_starstar(data[0])) 
        
            while rho > eps_t: 

                theta_star = weighted_sampling( theta_t_1, w_t_1 )
		 
                theta_starstar = multivariate_normal(theta_star, sig_t_1).rvs(size=1)
		while theta_starstar[1] < 0 :
			theta_star = weighted_sampling( theta_t_1, w_t_1 )
			theta_starstar = multivariate_normal(theta_star, sig_t_1).rvs(size=1)
		#print theta_starstar
                model_starstar = toolz.simulator( theta_starstar )
                rho = test_dist(data[1], model_starstar(data[0])) 

            #print theta_star, theta_starstar
            theta_t[:,i] = theta_starstar
	    #print sig_t_1
	    p_theta = toolz.prior()[0].pdf(theta_t[0,i]) * toolz.prior()[1].pdf(theta_t[1,i])
	    #print p_theta
	    pos_t = np.dstack((theta_t_1[0,:],theta_t_1[1,:]))
	    #print multivariate_normal(theta_t[:,i], sig_t_1).pdf(pos_t).shape , w_t_1.shape
            w_t[i] = p_theta / np.sum(w_t_1 * multivariate_normal(theta_t[:,i], sig_t_1).pdf(pos_t))
            #print test_dist(data[1], model_starstar(data[0])), w_t[i]
            
            #print 'For loop ', time.time() - start_time
        
        fig = plt.figure(1)
        sub = fig.add_subplot(111)
    	sub.hist2d(theta_t[0,:], theta_t[1,:], weights=w_t, bins=100, range=[[-2.0, 2.0],[0.0, 2.0]])
       
        plt.savefig("theta"+str(t)+".png")

        sig_t = 2.0 * np.cov(theta_t)
        np.savetxt(
		    ''.join(['theta_w_t', str(t), '.dat']), 
		    np.c_[theta_t[0,:], theta_t[1,:], w_t ]
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
    closest_theta = theta[:,cdf_closest_index]

    return closest_theta 


if __name__=='__main__': 
    # fake data
    data_x = uniform( -1.0, 2.0).rvs(size=1000)
    data_y = norm(0.0, 1.0).pdf(data_x)
    data = [data_x, data_y]

    fig = plt.figure(1)
    sub = fig.add_subplot(111)

    sub.scatter(data_x, data_y)
    fig.savefig('data.png')
    
    pmc_abc(data, 
        prior_shape = ['uniform', 'uniform'], 
        N = 1000,
        eps0 = 10.0, 
        T=20
        )
    
