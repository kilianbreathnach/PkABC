import numpy as np
from hod import N_cen, N_sat
from scipy.integrate import simps as inty


def nbar_g(lnm, hmf, HODparams):

    integrand = np.exp(lnm) * hmf * \
                 (N_cen(np.exp(lnm), HODparams[:2]) + \
                  N_sat(np.exp(lnm), HODparams))

    return inty(integrand, lnm)






