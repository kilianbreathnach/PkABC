import numpy as np

def lnprob(x , params):
    mu , lnsigma = params
    # Trivial improper prior: uniform in the log.
    if np.any((-1 > lnsigma) + (lnsigma > 1) + (-2 > mu) + (mu > 2)):
        return -np.inf
    lnprior = 0.0

    # Update the kernel and compute the lnlikelihood.
    
    return lnprior + lnlike(params , t , 1.)

def model(params, t):
    m , lsig = params
    sig2 = np.exp(lnsig) ** 2.
    amp = (np.sqrt(2. * np.pi * sig2)) ** -1.
    
    return amp * np.exp(-0.5 * (t - m) ** 2 / sig2)

def lnlike(p, t, yerr):

    y = (2. * np.pi) ** -.5 * np.exp (-1. * t ** 2. / 2)

    return -0.5 * np.sum(((y - model1(p, t))/yerr) ** 2)



lnsigma = np.linspace(-3 , 3 , 1000)
mu = np.linspace(-3 , 3 , 1000)
t = np.linspace(-1 , 1 , 100)


lnposterior = lnprob(t , (mu , lnsigma))

print lnposterior
