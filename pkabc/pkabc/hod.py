import numpy as np
from scipy.special import erf


def N_cen(m, hod):
    """
    Calculates the average number of central galaxies in halos of mass m.
    """
    return 0.5 * (1 + erf((np.log10(m) - hod[0]) / hod[1]))


def N_sat(m, hod):
    """
    Calculates the average number of satelite galaxies in halos of mass m.
    """
    return N_cen(m, hod[:2]) * (m / hod[3]) ** hod[2] * np.exp(- hod[4] / m)
