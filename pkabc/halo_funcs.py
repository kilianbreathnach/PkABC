import numpy as np
from scipy.special import erf



def N_cen(m, logMmin, sig_logm):
    """
    Calculates the average number of central galaxies in halos of mass m.
    """
    return 0.5 * (1 + erf((np.log10(m) - logMmin) / sig_logm))


def N_sat(m, logMmin, sig_logm, m_1, alph_sat, m_cut):
    """
    Calculates the average number of satelite galaxies in halos of mass m.
    """
    return N_cen(m, logMmin, sig_logm) * (m / m_1) ** alph_sat * np.exp(- m_cut / m)


def u_g(k, z, m):
    """
    Calculates the Fourier Transform of the NFW profile
    """
    pass

