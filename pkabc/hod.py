import numpy as np


def N_cen(m, logMmin, sig_logm):

    return 0.5 * (1 + erf((np.log10(m) - logMmin) / sig_logm))

def N_sat(m, logMmin, sig_logm, m_1, alph_sat, m_cut):

    return N_cen(m, logMmin, sig_logm) * (m / m_1) ** alph_sat * np.exp(- m_cut / m)
