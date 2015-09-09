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

import triangle
from distance import test_dist
from parameters import Params
from test_sim import simz

class PmcAbc(object): 
    
    def __init__(self, data, N = 1000, eps0 = 0.01, T = 20): 
        """ Class taht describes PMC-ABC 
        """
        self.data = data
        self.N = N 
        self.eps0 = eps0
        self.T = T
    
    def prior_param(self, 
            param_dict= {
                    'mu': {'shape': 'uniform', 'min': -1, 'max': 1}, 
                    'sigma': { 'shape': 'uniform', 'min': -2.0, 'max': 2.0}
                    }): 
        """ Pass priors of parameters in theta
        """

        self.param_obj = Params(param_dict)     # parameter object

        self.param_names = param_dict.keys()

        self.n_params = len(param_dict.keys())  # number of parameters in theta

    def priors_sample(self): 
        """ Sample from priors derived from parameter object
        """
        
        theta_star = np.zeros(self.n_params)

        for i in xrange(self.n_params): 
            theta_star[i] = self.param_obj.prior()[i].rvs(size=1)[0]

        return theta_star
                
    
    def prior_of_priors(self, tt): 
        """ Multiply priors of multile dimensions
        p(theta) = p(theta_0) * p(theta_1) * ... * p(theta_n_params)
        """
        for i in xrange(self.n_params): 
            try: 
                p_theta *= self.param_obj.prior()[i].pdf(tt[i])        
            except UnboundLocalError: 
                p_theta = self.param_obj.prior()[i].pdf(tt[i])        
            
        return p_theta
    
    def initial_pool(self): 
        """ Initial pool of pmc_abc
        """

        self.prior_param()  # first run prior parameters

        self.t = 0 
        self.theta_t = np.zeros((self.n_params, self.N))
        self.w_t = np.zeros((self.N))
        rhos = [] 

        for i in xrange(self.N): 
            
            theta_star = self.priors_sample()
    
            model = simz( theta_star )
            rho = test_dist(self.data, model)

            while rho > self.eps0: 

                theta_star = self.priors_sample()
                
                model = simz( theta_star )

                rho = test_dist(self.data, model)
             
            self.theta_t[:,i] = theta_star

            self.w_t[i] = 1.0/np.float(self.N)

            rhos.append(rho)

        self.sig_t = 2.0 * np.cov( self.theta_t )    # covariance matrix

        self.writeout()
        self.plotout()

        return np.array(rhos)
    
    def pmc_abc(self): 
        """
        """

        rhos = self.initial_pool()

        while self.t < self.T: 

            eps_t = np.percentile(rhos, 75)

            print 'Epsilon t', eps_t

            theta_t_1 = self.theta_t.copy()
            w_t_1 = self.w_t.copy()
            sig_t_1 = self.sig_t.copy()
            rhos = [] 

            for i in xrange(self.N): 

                theta_star = weighted_sampling( theta_t_1, w_t_1 ) 
                theta_starstar = multivariate_normal( theta_star, sig_t_1 ).rvs(size=1)
               
                model_starstar = simz( theta_starstar )

                rho = test_dist(data, model_starstar) 
            
                while rho > eps_t: 

                    theta_star = weighted_sampling( theta_t_1, w_t_1 )
                    theta_starstar = multivariate_normal(theta_star, sig_t_1).rvs(size=1)

                    model_starstar = simz( theta_starstar )

                    rho = test_dist(data, model_starstar) 
                
                rhos.append(rho)

                self.theta_t[:,i] = theta_starstar

                p_theta = self.prior_of_priors(theta_starstar)
                
                pos_t = np.dstack(theta_t_1)
                self.w_t[i] = p_theta / np.sum(w_t_1 * multivariate_normal(self.theta_t[:,i], sig_t_1).pdf(pos_t))
        
            self.sig_t = 2.0 * np.cov(self.theta_t)
            self.t += 1 
            
            self.writeout()

            self.plotout()

        return None
    
    def writeout(self): 
        """ Write out theta_t and w_t
        """

        out_file = ''.join(['theta_w_t', str(self.t), '.dat'])

        data_list = [] 
        for i in xrange(self.n_params): 
            data_list.append( self.theta_t[i,:] ) 
        data_list.append(self.w_t)

        np.savetxt(
                out_file, 
                (np.vstack(np.array(data_list))).T, 
                delimiter='\t'
                )

        return None 

    def plotout(self): 
        """ Triangle plot the things 
        """
        figure = triangle.corner(
               (self.theta_t).T, 
               labels = self.param_names, 
               weights = self.w_t, 
               show_titles=True, 
               range = [[-2,2],[0,1]], 
               title_args={"fontsize": 12}
               ) 
    
        figure.gca().annotate(
                str(self.t), 
                xy=(0.5, 1.0), 
                xycoords="figure fraction",
                xytext=(0, -5), 
                textcoords="offset points",
                ha="center", 
                va="top"
                ) 
        figure.savefig("triangle_theta_t"+str(self.t)+".png")
        plt.close()

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
    data = {'input': data_x, 'output': data_y}

    fig = plt.figure(1)
    sub = fig.add_subplot(111)
    sub.scatter(data_x, data_y)
    fig.savefig('data.png')
    plt.close()
    
    pmcabc_test = PmcAbc(data, N=5000, eps0 = 2.0, T = 10)
    pmcabc_test.pmc_abc()
