import numpy as np
import abcpmc
import matplotlib.pyplot as plt
from interruptible_pool import InterruptiblePool
import time
plt.switch_backend("Agg")


from halotools.empirical_models import Zheng07


model = Zheng07(threshold = -20.5)
print 'Data HOD Parameters ', model.param_dict

n_avg = 5
avg_nz, avg_corr = 0., 0.

for i in xrange(n_avg): 
    
    model.populate_mock()
    
    # number density
    avg_nz += model.mock.number_density
    
    # 6th element of xi(r) array
    r, xi_r = model.mock.compute_galaxy_clustering(N_threads=1)

    try:
    	avg_xi += xi_r
    except NameError:
	avg_xi = xi_r

    avg_corr += xi_r[6]

avg_nz /= np.float(n_avg)
avg_corr /= np.float(n_avg)
avg_xi /=np.float(n_avg)

data = [avg_nz, avg_xi]
data_hod = np.array([1.12 , 12.3 , 0.21 , 11.84 , 13.58])

class HODsim(object): 
    
    def __init__(self): 
        self.model = Zheng07(threshold = -20.5)
    
    def sum_stat(self, theta_star):
        
        self.model.param_dict['alpha'] = theta_star[0]
        self.model.param_dict['logMmin'] = theta_star[1]     
        self.model.param_dict['sigma_logM'] = theta_star[2]     
        self.model.param_dict['logM0'] = theta_star[3]
        self.model.param_dict['logM1'] = theta_star[4]     

	#print self.model.param_dict
        
        self.model.populate_mock()
        nz = self.model.mock.number_density
        
        r, xi_r = model.mock.compute_galaxy_clustering(N_threads = 1)
        
        #xi = xi_r[6]
        #xi = xi_r[7]
        
        return [nz, xi_r]

ourmodel = HODsim()
simz = ourmodel.sum_stat


def distance(data, model, type = 'sum_stat'): 
    
    if type == 'sum_stat': 
        dist_nz = np.abs(data[0] - model[0])/data[0]
        dist_xi = np.sum(np.abs(data[1] - model[1])/data[1])
        
        dist = dist_nz + dist_xi 
        
    return dist 


"""Prior"""

from scipy.stats import uniform
from scipy.stats import norm

class Prior(object): 
    def __init__(self, prior_dict): 
        self.prior_dict = prior_dict.copy()
    
    def prior(self): 
        priorz = [] 
        for key in self.prior_dict.keys(): 
            
            prior_key = self.prior_dict[key]
            
            if prior_key['shape'] == 'uniform': 
                
                loc = prior_key['min']
                scale = prior_key['max'] - prior_key['min']
                
                priorz.append( uniform(loc, scale))
            
            elif prior_key['shape'] == 'gauss':
                
                loc = prior_key['mean']
                scale = prior_key['stddev']
                
                priorz.append( norm(loc, scale) )
                
        return priorz


#prior_dict = { 
#    'alpha': {'shape': 'uniform', 'min': 0.8, 'max': 1.2},
#    'm_1'  : {'shape': 'uniform', 'min': 13., 'max': 14.},   
#    'm_min': {'shape': 'uniform', 'min': 11.5, 'max': 12.5}
#}

prior_dict = { 
    'alpha': {'shape': 'uniform', 'min': 1.02,  'max': 1.22},
    'm_min'  : {'shape': 'uniform', 'min': 12.1,  'max': 12.5},   
    'sigma': {'shape': 'uniform', 'min': .19,  'max': .23},
    'm0': {'shape': 'uniform', 'min': 11.64,  'max': 12.04},
    'm1': {'shape': 'uniform', 'min': 13.28,  'max': 13.88}
}
n_params = len(prior_dict.keys())
prior_obj = Prior(prior_dict) 



def prior_sampler(): 
    """ Sample prior distribution and return theta_star 
    """
    theta_star = np.zeros(n_params)
    
    for i in xrange(n_params): 
        np.random.seed()
        theta_star[i] = prior_obj.prior()[i].rvs(size=1)[0]
        
    return theta_star

def pi_priors(tmp_theta): 
    for i in xrange(n_params): 
        try:
            p_theta *= prior_obj.prior()[i].pdf(tmp_theta[i])
        except UnboundLocalError: 
            p_theta = prior_obj.prior()[i].pdf(tmp_theta[i])
            
    return p_theta 



import corner 
import seaborn as seabreeze

def plot_thetas(theta , w , t): 
    fig = corner.corner(
        theta.T, weights = w.flatten() , truths= data_hod,
        truth_color="red", plot_datapoints=True, fill_contours=True, levels=[0.68, 0.95], 
                color='b', bins=80, smooth=1.0, 
        range=[(0.85, 1.35), (12.0 , 12.5), (0.15, 0.25) , (11.5 , 12.05) , (13. , 14.)], 
        labels=[r"$\alpha$", r"$\log M_{min}$" , r"$\sigma$" , r"$\log M_{0}$" , r"$\log M_{1}$"]
        )
    
    plt.savefig("/home/mj/public_html/scatter_hod_flat_t"+str(t)+".png")
    plt.close()
    np.savetxt("/home/mj/public_html/theta_hod_flat_t"+str(t)+".dat" , theta.T)
    
    np.savetxt("/home/mj/public_html/w_hod_flat_t"+str(t)+".dat" , w.T)

#plot_datapoints=True, fill_contours=True, levels=[0.68, 0.95], 
#                color='b', bins=80, smooth=1.0


N_threads = 20
N_particles = 1000
N_iter = 20
eps0 = 2.0

def initial_pool_sampling(i_particle): 
    """ Sample theta_star from prior distribution for the initial pool
    """
    theta_star = prior_sampler()
    model_theta = simz(theta_star)
    
    rho = distance(data, model_theta)
    
    while rho > eps0: 
        
        theta_star = prior_sampler()
        model_theta = simz(theta_star)
        
        rho = distance(data, model_theta)
        
    pool_list = [np.int(i_particle)]
    for i_param in xrange(n_params): 
        pool_list.append(theta_star[i_param])
        
    pool_list.append(1./np.float(N_particles))
    pool_list.append(rho)
    
    return np.array(pool_list)



def initial_pool():

    pool = InterruptiblePool(processes = N_threads)
    mapfn = pool.map
    args_list = [i for i in xrange(N_particles)]
    results = mapfn(initial_pool_sampling, args_list)
    
    pool.close()
    pool.terminate()
    pool.join()
    
    results = np.array(results).T
    theta_t = results[1:n_params+1,:]
    w_t = results[n_params+1,:]
    rhos = results[n_params+2,:]
    sig_t = np.cov(theta_t)
    
    return theta_t, w_t, rhos, sig_t


def weighted_sampling(theta, w): 
    """ Given array of thetas and their corresponding weights, sample
    """
    w_cdf = w.cumsum()/w.sum() # normalized CDF
    
    np.random.seed()
    rand1 = np.random.random(1)
    cdf_closest_index = np.argmin( np.abs(w_cdf - rand1) )
    closest_theta = theta[:, cdf_closest_index]
    
    return closest_theta

def better_multinorm(theta_stst, theta_before, cov): 
    n_par, n_part = theta_before.shape
    
    sig_inv = np.linalg.inv(cov)
    x_mu = theta_before.T - theta_stst

    nrmliz = 1.0 / np.sqrt( (2.0*np.pi)**n_par * np.linalg.det(cov))

    multinorm = nrmliz * np.exp(-0.5 * np.sum( (x_mu.dot(sig_inv[None,:])[:,0,:]) * x_mu, axis=1 ) )

    return multinorm


from scipy.stats import multivariate_normal 

def importance_pool_sampling(args): 
    # args = [i_particle, theta_t_1, w_t_1, sig_t_1, eps_t]
    i_particle = args[0]
    theta_t_1 = args[1]
    w_t_1 = args[2]
    sig_t_1 = args[3]
    eps_t = args[4]
    
    theta_star = weighted_sampling(theta_t_1, w_t_1)
    
    np.random.seed()
    # perturbed theta (Double check)    
    theta_starstar = multivariate_normal( theta_star, sig_t_1 ).rvs(size=1)
    model_starstar = simz(theta_starstar)
    
    rho = distance(data, model_starstar)
    
    while rho > eps_t:
        theta_star = weighted_sampling(theta_t_1, w_t_1)
        theta_starstar = multivariate_normal( theta_star, sig_t_1 ).rvs(size=1)
        model_starstar = simz(theta_starstar)
        
        rho = distance(data, model_starstar)
    
    p_theta = pi_priors(theta_starstar)

    w_starstar = p_theta/np.sum( w_t_1 * better_multinorm(theta_starstar, theta_t_1, sig_t_1) )    
    
    pool_list = [np.int(i_particle)]
    for i_p in xrange(n_params): 
        pool_list.append(theta_starstar[i_p])
    pool_list.append(w_starstar)
    pool_list.append(rho)
    
    return pool_list 
    
def pmc_abc(N_threads = N_threads): 
    
    # initial pool
    theta_t, w_t, rhos, sig_t = initial_pool()
    t = 0 # iternation number
    
    plot_thetas(theta_t , w_t, t)
    
    plt.savefig("/home/mj/public_html/scatter_hod_gaussian_t"+str(t)+".png")
    plt.close()
    
    while t < N_iter: 
        
        eps_t = np.percentile(rhos, 75)
        print 'New Distance Threshold Eps_t = ', eps_t
        
        theta_t_1 = theta_t.copy()
        w_t_1 = w_t.copy()
        sig_t_1 = sig_t.copy()
    
        """these lines are borrowed from initial sampling to double-check multiprocessing"""
        #pool = InterruptiblePool(processes = N_threads)
    	#mapfn = pool.map
    	#args_list = [i for i in xrange(N_particles)]
    	#results = mapfn(initial_pool_sampling, args_list)
    
    	#pool.close()
    	#pool.terminate()
    	#pool.join()

        pool = InterruptiblePool(processes = N_threads)
        mapfn = pool.map
        args_list = [[i, theta_t_1, w_t_1, sig_t_1, eps_t] for i in xrange(N_particles)]
        #results = [] 
        #for args in args_list: 
        #    pool_sample = importance_pool_sampling(args)
        #    results.append( pool_sample )
        results = mapfn(importance_pool_sampling, args_list)
        pool.close()
        pool.terminate()
        pool.join()
        
        sig_t = np.cov(theta_t)
                 
        results = np.array(results).T
        theta_t = results[1:n_params+1,:]
        w_t = results[n_params+1,:]
        rhos = results[n_params+2,:]
        sig_t = np.cov(theta_t)
        
        t += 1
        
        plot_thetas(theta_t, w_t , t)
        plt.savefig("/home/mj/public_html/scatter_hod_gaussian_t"+str(t)+".png")
        plt.close()
pmc_abc()

