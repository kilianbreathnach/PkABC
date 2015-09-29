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
import warnings
import seaborn as sns

import corner 
from distance import test_dist
from parameters import Params
from simulator import Simul

from interruptible_pool import InterruptiblePool

def unwrap_self_importance_sampling(arg, **kwarg):

    return PmcAbc.importance_sampling(*arg, **kwarg)

def unwrap_self_initial_sampling(arg, **kwarg):
    
    return PmcAbc.initial_sampling(*arg, **kwarg)

class PmcAbc(object): 
    
    def __init__(self, data, simulator, 
            prior_dict = {}, 
            N = 1000, 
            eps0 = 0.01, 
            T = 20, 
            Nthreads = 10): 
        """ Class taht describes PMC-ABC 
        """
        self.data = data
        self.N = N 
        self.eps0 = eps0
        self.T = T
        self.Nthreads = Nthreads

<<<<<<< HEAD:abc/pmc_abc.py
    def prior_param(self, 
            param_dict= {
                    'sigma': {'shape': 'uniform', 'min': 0.1, 'max': 0.4}, 
                    'm_min': { 'shape': 'uniform', 'min': 11.0, 'max': 13.0}
		               
                   , 'alpha': { 'shape': 'uniform', 'min': 1.059, 'max': 1.061}

                   , 'm0': { 'shape': 'uniform', 'min': 11.37, 'max': 11.39}
                   , 'm1': { 'shape': 'uniform', 'min': 13.30, 'max': 13.32}
         }): 
=======
        # simulator function has to be a function of theta_star
        self.simz = simulator
        
        self.prior_param(param_dict = prior_dict)  # first run prior parameters


    def prior_param(self, param_dict={}): 
>>>>>>> 20e172ad26c43b1397206341941735dd6a9129eb:eybc/pmc_abc.py
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
            #print 1 
            np.random.seed()
            theta_star[i] = self.param_obj.prior()[i].rvs(size=1)[0]

        return theta_star
                
    
    def prior_of_priors(self, tt): 
        """ Multiply priors of multile dimensions
        p(theta) = p(theta_0) * p(theta_1) * ... * p(theta_n_params
        """
        for i in xrange(self.n_params): 
            try: 
                p_theta *= self.param_obj.prior()[i].pdf(tt[i])        
            except UnboundLocalError: 
                p_theta = self.param_obj.prior()[i].pdf(tt[i])        
            
        return p_theta
    
    def initial_sampling(self, params):
        """Wrapper for parallelized initial pool sampling
        """
        i = params
        
        theta_star = self.priors_sample()
        model = self.simz( theta_star )
        rho = test_dist(self.data, model)
        print rho  
        while rho > self.eps0: 

            theta_star = self.priors_sample()
                
            model = self.simz( theta_star )

            rho = test_dist(self.data, model)
            #print 1 
        data_list = [np.int(i)]
        print "done"
        for i_param in xrange(self.n_params): 
            data_list.append(theta_star[i_param])

        data_list.append(1./np.float(self.N))
        data_list.append(rho)

	return np.array(data_list)   

    def initial_pool(self):
        """
        Creating the initial pool
        """
        self.t = 0 
        self.theta_t = np.zeros((self.n_params, self.N))
        self.w_t = np.zeros((self.N))
        self.rhos = np.zeros((self.N)) 

    
        pool = InterruptiblePool(self.Nthreads)    
        mapfn = pool.map
        args_list = [(i) for i in xrange(self.N)]
        unwrap_self_initial_sampling(zip([self]*len(args_list), args_list)[0])
        results = mapfn(unwrap_self_initial_sampling, zip([self]*len(args_list), args_list))
        pool.close()
        pool.terminate()
        pool.join()
        
 	pars = np.array(results).T
        self.theta_t = pars[1:self.n_params+1,:]
        self.w_t     = pars[self.n_params+1,:]
        self.rhos    = pars[self.n_params+2,:]
        
        self.sig_t = 2.0 * np.cov( self.theta_t )    # covariance matrix

        self.writeout()
        self.plotout()

        return np.array(self.rhos)   
        
    def importance_sampling(self, params): 
        """ Wrapper for parallelized importance sampling
        """

        i_part = params

        theta_star = weighted_sampling( self.theta_t_1, self.w_t_1 ) 
        np.random.seed()
        theta_starstar = multivariate_normal( theta_star, self.sig_t_1 ).rvs(size=1)
       
        model_starstar = self.simz( theta_starstar )

        rho = test_dist(self.data, model_starstar) 
    
        while rho > self.eps_t: 

            theta_star = weighted_sampling( self.theta_t_1, self.w_t_1 )
            theta_starstar = multivariate_normal(theta_star, self.sig_t_1).rvs(size=1)

            model_starstar = self.simz( theta_starstar )

            rho = test_dist(self.data, model_starstar) 

        #print theta_star, theta_starstar

        p_theta = self.prior_of_priors(theta_starstar)
        
        pos_t = np.dstack(self.theta_t_1)
        tmp_w_t = p_theta / np.sum(self.w_t_1 * multivariate_normal(self.theta_t[:,i_part], self.sig_t_1).pdf(pos_t))
        
        data_list = [np.int(i_part)]
        for i_param in xrange(self.n_params): 
            data_list.append(theta_starstar[i_param])
        data_list.append(tmp_w_t)
        data_list.append(rho)
        
        return  np.array(data_list) 
    
    def pmc_abc(self): 
        """
        """

        self.rhos = self.initial_pool()

        while self.t < self.T: 

            self.eps_t = np.percentile(self.rhos, 75)

            print 'Epsilon t', self.eps_t

            self.theta_t_1 = self.theta_t.copy()
            self.w_t_1 = self.w_t.copy()
            self.sig_t_1 = self.sig_t.copy()
            
            pool = InterruptiblePool(self.Nthreads)
            mapfn = pool.map
            args_list = [ i for i in xrange(self.N) ] 
            results = mapfn(unwrap_self_importance_sampling, zip([self]*len(args_list), args_list))

            pool.close()
            pool.terminate()
            pool.join()
	    
            pars = np.array(results).T
            self.theta_t = pars[1:self.n_params+1,:].copy()
            self.w_t     = pars[self.n_params+1,:].copy()
            self.rhos    = pars[self.n_params+2,:].copy()

            self.sig_t = 2.0 * np.cov(self.theta_t)
            self.t += 1 
            
            self.writeout()

            self.plotout()

        return None
    
    def writeout(self): 
        """ Write out theta_t and w_t
        """

        out_file = ''.join(['teta_w_t', str(self.t), '.dat'])

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

<<<<<<< HEAD:abc/pmc_abc.py
    def plotout(self, plot_type = 'triangle'): 
        """ Triangle plot the things 
        """
        print self.theta_t.T.shape
        if plot_type == 'triangle': 
=======
    def plotout(self, plot_type = 'seabreeze'): 
        """ Triangle plot the things 
        """
        if plot_type == 'seabreeze':

            figure = sns.jointplot(x = self.theta_t[0,:], y = self.theta_t[1,:], kind = 'kde', 
                    style = 'white', weights = self.w_t, 
                    xlim = [-1.0, 1], 
                    ylim = [0.0, 2.0]
                    )
            
            plt.savefig("seabreeze_theta_t"+str(self.t)+".png")
            plt.close()

        elif plot_type == 'triangle': 
            # Clunky based on which version of corner.py you have
            # Clunky based on which version of corner.py you have
            # Clunky based on which version of corner.py you have
            # Clunky based on which version of corner.py you have
>>>>>>> 20e172ad26c43b1397206341941735dd6a9129eb:eybc/pmc_abc.py
            figure = triangle.corner(
                   (self.theta_t).T, 
                   labels = self.param_names, 
                   weights = self.w_t, 
                   show_titles=False,
                   range = [[11., 13.] , [.1 , .4]], 
                   title_args={"fontsize": 12},
                   smooth=False, 
                   color = "k",
                   scale_hist=False
                   )
            figure.savefig("triangle_theta_t"+str(self.t)+".png")
            plt.close()

        elif plot_type == 'scatter': 
            
            if len(self.theta_t[:,0]) != 2: 
                warnings.warn("Can only plot two axes on scatter plot. No plot generated")
                return 

            figure = plt.figure(1)
            sub = figure.add_subplot(111)
            sub.scatter(self.theta_t[0,:], self.theta_t[1,:]) 
            sub.set_xlim([-1.0, 1.0])
            sub.set_ylim([0.8, 1.5])
            
            figure.savefig("scatter_theta_t"+str(self.t)+".png")
            plt.close()

def weighted_sampling(theta, w): 
    """ Weighted sampling
    """

    w_cdf = w.cumsum()/w.sum()
    
    np.random.seed()
    rand1 = np.random.random(1)
    cdf_closest_index = np.argmin(np.abs(w_cdf - rand1))
    #min(xrange(len(w_cdf)), key = lambda iii: abs(w_cdf[iii]-rand1[0]))
    closest_theta = theta[:,cdf_closest_index]

    return closest_theta 

if __name__=='__main__': 
    pass
    # fake data
    #data_x = uniform( -1.0, 2.0).rvs(size=1000)
    #data_y = norm(0.0, 1.0).pdf(data_x)
    #data = {'input': data_x, 'output': data_y}

    ##fig = plt.figure(1)
    ##sub = fig.add_subplot(111)
    ##sub.scatter(data_x, data_y)
    ##fig.savefig('data.png')
    plt.close()
    data = {'output': [0.0047808, 6.76599366e+02]}

<<<<<<< HEAD:abc/pmc_abc.py
    pmcabc_test = PmcAbc(data, N=10, eps0 = 0.1, T = 20, Nthreads=20)
=======
    modeel = Simul()
    simz = modeel.sum_stat

    pmcabc_test = PmcAbc(data, simz, 
            prior_dict = {
                'logMmin': {'shape': 'uniform', 'min': 11.5, 'max': 12.5}, 
                'sigma_logM': {'shape': 'uniform', 'min': 0.2, 'max': 0.3}
                }, 
            N=100, eps0 = 0.5, T = 10, Nthreads=3)
>>>>>>> 20e172ad26c43b1397206341941735dd6a9129eb:eybc/pmc_abc.py
    pmcabc_test.pmc_abc()
