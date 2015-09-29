import time
import warnings
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import uniform
from scipy.stats import norm
from scipy.stats import multivariate_normal
from interruptible_pool import InterruptiblePool
import corner
from distance import test_dist
from parameters import Params
from simulator import Simul


def unwrap_self_importance_sampling(arg, **kwarg):

    return PmcAbc.importance_sampling(*arg, **kwarg)


def unwrap_self_initial_sampling(arg, **kwarg):

    return PmcAbc.initial_sampling(*arg, **kwarg)


class PmcAbc(object):
        """ Class taht describes PMC-ABC
        """

    def __init__(self, data, simulator, prior, N = 1000, eps0 = 0.01, T = 20, Nthreads = 10):
        self.data = data
        self.simulator = simulator
        self.prior = prior
        self.n_params = prior.n_params
        self.N = N
        self.eps0 = eps0
        self.T = T
        self.Nthreads = Nthreads


    def priors_sample(self):
        """ Sample from priors derived from parameter object
        """

        theta_star = np.zeros(self.n_params)

        for i in xrange(self.n_params):
            np.random.seed()
            theta_star[i] = self.prior.prior()[i].rvs(size=1)[0]

        return theta_star


    def prior_of_priors(self, tt):
        """ Multiply priors of multile dimensions
        p(theta) = p(theta_0) * p(theta_1) * ... * p(theta_n_params)
        """
        for i in xrange(self.prior.n_params):
            try:
                p_theta *= self.prior.prior()[i].pdf(tt[i])
            except UnboundLocalError:
                p_theta = self.prior.prior()[i].pdf(tt[i])

        return p_theta


    def initial_sampling(self, params):
        """Wrapper for parallelized initial pool sampling
        """
        i = params

        theta_star = self.priors_sample()
        model = simulator.predict( theta_star )
        rho = test_dist(self.data, model)

        while rho > self.eps0:

            theta_star = self.priors_sample()

            model = simulator.predict( theta_star )

            rho = test_dist(self.data, model)

        data_list = [np.int(i)]

        for i_param in xrange(self.n_params):
            data_list.append(theta_star[i_param])

        data_list.append(1./np.float(self.N))
        data_list.append(rho)

	return np.array(data_list)


    def initial_pool(self):
        """
        Creating the initial pool
        """
#         self.prior_param()  # first run prior parameters

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

        model_starstar = simulator.predict( theta_starstar )

        rho = test_dist(data, model_starstar)

        while rho > self.eps_t:

            theta_star = weighted_sampling( self.theta_t_1, self.w_t_1 )
            theta_starstar = multivariate_normal(theta_star, self.sig_t_1).rvs(size=1)

            model_starstar = simulator.predict( theta_starstar )

            rho = test_dist(data, model_starstar)

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


    def plotout(self, plot_type = 'scatter'):
        """ Triangle plot the things
        """
        if plot_type == 'triangle':
            # Clunky based on which version of corner.py you have
            # Clunky based on which version of corner.py you have
            # Clunky based on which version of corner.py you have
            # Clunky based on which version of corner.py you have
            figure = triangle.corner(
                   (self.theta_t).T,
                   labels = self.param_names,
                   weights = self.w_t,
                   show_titles=True,
                   title_args={"fontsize": 12},
                   smooth=False
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

        elif plot_type == 'scatter':

            if len(self.theta_i[:,0]) != 2:
                warnings.warn("Can only plot two axes on scatter plot. No plot generated")
                return

            figure = plt.figure(1)
            sub = figure.add_subplot(111)
            sub.scatter(self.theta_t[0,:], self.theta_t[1,:])

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
