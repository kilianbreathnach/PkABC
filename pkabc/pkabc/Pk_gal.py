import numpy as np
from scipy.integrate import simps as inty
from halofuncs import N_cen, N_sat
from matter import Matter
from nbar import nbar_g

class Pwspec:
    """
    Object class for computing the analytic galaxy power spectrum as in the
    appendix to Schneider et al, 2006.
    """
    def __init__(self, k, zvec, Mmin, Mmax):

        """
        - k is a 1-d numpy array of wavenumbers. (unit = Mpc^-1)
        - zvec is 1-d array of central values of the redshift bins
        - Mmin (Mmax) is the minimum (maximum) halo mass. (unit = Msun/h)
        """
        self.k = k
        self.zvec = zvec
        self.lnm = np.logspace(np.log(Mmin), np.log(Mmax), base=np.exp(1))


    def set_universe(self, Om=0.3, OL=0.7, ns=0.96,
                     sig_8=0.82, h=0.673, T_cmb=2.725
                     k_min = 1.e-3, k_max = 2.e3, dk = 0.05,
                     lnM_min=np.log(1e11), lnM_max=np.log(1e15),
                     dlnM=np.log(5e9),
                     transfer_fit="EH",
                     hmf_fit="Tinker",
                     bias_fit="Tinker"):
        """
        Set the cosmological parameters and fitting functions used for the
        power spectrum calculation. This attribute is an instance of the
        Matter class defined in matter.py, which has functions for computing
        the required dark matter statistics.
        """
        self.universe = Matter(Om=0.3, OL=0.7, ns=0.96,
                               sig_8=0.82, h=0.673, T_cmb=2.725
                               k_min = 1.e-3, k_max = 2.e3, dk = 0.05,
                               lnM_min=np.log(1e11), lnM_max=np.log(1e15),
                               dlnM=np.log(5e9),
                               transfer_fit="EH",
                               hmf_fit="Tinker",
                               bias_fit="Tinker")


    def precompute(self):
        """
        Using the k-vector, redshift bins and mass range given to the Pwspec
        object, this function precomputes a list for the redshift values of
        halo mass functions, NFW fourier transforms and halo biases for use in
        computing the galaxy power spectra for various HOD models. These are
        computed using functions from the Matter class which the universe
        attribute is an instance of.
        """
        if not self.universe:
            print "You must first set the universe parameters"

        self.hmflist = []
        self.uglist = []
        self.biaslist = []

        for z in self.zvec:
            self.hmflist.append(self.universe.hmf(self.lnm, z))
            self.uglist.append(self.universe.u_g(self.k, self.lnm, z))
            self.biaslist.append(self.universe.bias(self.lnm, z))


    def set_HOD(self, params):
        """
        expects iterable of parameter values in the form
        [Mmin, siglogm, alpha, M_1, Mcut]
        where masses are in log10
        """
        self.HOD = params
        self.compute_nbarlist()


    def compute_nbarlist(self):
        """
        Use the halo mass function and HOD parameters to compute the galaxy
        number density at a given redshift.
        """
        self.nbarlist = []

        for z in self.zvec:
            self.nbarlist.append(nbar_g(self.lnm, self.get_hmf(z), self.HOD))


    def zind(self, z):

        return np.where(self.zvec == z)


    def get_nbar(self, z):

        return self.nbarlist[self.zind(z)]


    def get_hmf(self, z):

        return self.hmflist[self.zind(z)]


    def get_ug(self, z):

        return self.uglist[self.zind(z)]


    def get_Plin(self, z1, z2):
        """
        Compute variance per logarithmic interval in k as per Schneider et al
        """
        g1 = self.universe.gf(z1)
        g2 = self.universe.gf(z2)

        return dHsq * g1 * g2 * Pspec_norm


    def get_bias(self, z):

        return self.biaslist[self.zind(z)]


    def _1h_int(self, z):
        """
        The integrand for the 1-halo term of the galaxy power spectrum.
        """
        return np.exp(self.lnm) * self.get_hmf(z) * \
                (N_sat(np.exp(self.lnm), self.HOD) ** 2 * \
                 self.get_ug(z) ** 2 + \
                 2 * N_cen(np.exp(self.lnm), self.HOD) * \
                     N_sat(np.exp(self.lnm), self.HOD) * \
                     self.get_ug(z))


    def P_1h(self, z1, z2):
        """
        1-halo term in the galaxy power spectrum
        """
        if z1 != z2:
            return 0.
        else:
            return (1. /  self.get_nbar(z1) ** 2) * \
                    inty(self._1h_int(z1), self.lnm)


    def _I2inty(self, z):
        """
        The integrand for the 2-halo term of the galaxy power spectrum.
        """
        return np.exp(self.lnm) * self.get_hmf(z) * self.get_bias(z)\
                (N_cen(np.exp(self.lnm), self.HOD) + \
                 N_sat(np.exp(self.lnm), self.HOD) * self.get_ug(z))


    def I_2(self, z):
        """
        Integral 2-halo term in the galaxy power spectrum
        """
        return (1. / self.get_nbar(z)) * \
                inty(self._I2inty, self.lnm)


    def P_2h(self, z1, z2):
        """
        2-halo term in the galaxy power spectrum
        """
        # TODO : after finishing matter.py, modify this accordingly

        return self.get_Plin(z1, z2) * I_2(z1) * I_2(z2)


    def P_g(self, z1, z2):
        """
        Galaxy power spectrum
        Auto (Cross) power spectrum if z1=z2(z1!=z2)
        """
        if not (z1 in self.zvec) * (z2 in self.zvec):
            print "These redshifts are not in the precomputed redshift values"
            break

        return self.P_1h(z1, z2) + self.P_2h(z1, z2)
