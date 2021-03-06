import numpy as np
from scipy.integrate import simps as inty
from hod import N_cen, N_sat
from matter import Matter
from nbar import nbar_g

class Pwspec:
    """
    Object class for computing the analytic galaxy power spectrum as in the
    appendix to Schneider et al, 2006 and Van der Bosch et al, 2013.
    """
    def __init__(self, zvec, k=None, Mmin=1e11, Mmax=1e15):

        """
        - k is a 1-d numpy array of wavenumbers. (unit = Mpc^-1)
        - zvec is 1-d array of central values of the redshift bins
        - Mmin (Mmax) is the minimum (maximum) halo mass. (unit = Msun/h)
        """
        self.k = k
        self.zvec = zvec
#        self.lnm = np.logspace(np.log(Mmin), np.log(Mmax), base=np.exp(1))
        self.set_universe()
	self.precompute()

    def set_universe(self, Om=0.3, OL=0.7, ns=0.96,
                     sig_8=0.81, h=1.0, T_cmb=2.725, h_transf=0.7,
                     k=None, k_min=1.e-3, k_max=2.e3, kbins=1000,
                     lnM_min=np.log(1e11), lnM_max=np.log(1e15),
                     transfer_fit="EH",
                     hmf_fit="Tinker",
                     bias_fit="Tinker"):
        """
        Set the cosmological parameters and fitting functions used for the
        power spectrum calculation. This attribute is an instance of the
        Matter class defined in matter.py, which has functions for computing
        the required dark matter statistics.
        """
        self.universe = Matter(Om=Om, OL=OL, ns=ns,
                               sig_8=sig_8, h=h, T_cmb=T_cmb,
                               h_transf=h_transf,
                               k=k,
                               k_min=k_min, k_max=k_max, kbins=kbins,
                               lnM_min=lnM_min, lnM_max=lnM_max,
                               lnMbins=1000,
                               transfer_fit=transfer_fit,
                               hmf_fit=hmf_fit,
                               bias_fit=bias_fit)

        self.lnm = self.universe.lnM


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
            self.hmflist.append(self.universe.hmf(z))
            self.uglist.append(self.universe.u_g(z))
            self.biaslist.append(self.universe.bias(z))


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
            self.nbarlist.append(nbar_g(self.universe.lnM, self.get_hmf(z), self.HOD))


    def zind(self, z):

        return np.where(self.zvec == z)[0]


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
        psp = self.universe.normal_p0_lin()

        return g1 * g2 * psp


    def get_bias(self, z):

        return self.biaslist[self.zind(z)]


    def _1h_int(self, z, k_i):
        """
        The integrand for the 1-halo term of the galaxy power spectrum.
        """
        return self.get_hmf(z) * \
                (N_sat(np.exp(self.lnm), self.HOD) ** 2 * \
                self.get_ug(z)[k_i, :] ** 2 + \
                 2 * N_cen(np.exp(self.lnm), self.HOD) * \
                     N_sat(np.exp(self.lnm), self.HOD) * \
                     self.get_ug(z)[k_i, :])


    def P_1h(self, z1, z2):
        """
        1-halo term in the galaxy power spectrum
        """
        if z1 != z2:
            return 0.
        else:
            plist = []
            for i in range(len(self.universe.k)):
                plist.append((1. /  self.get_nbar(z1) ** 2) * \
                              inty(self._1h_int(z1, i), self.lnm))

            return np.array(plist)


    def _I2inty(self, z, k_i):
        """
        The integrand for the 2-halo term of the galaxy power spectrum.
        """
        return self.get_hmf(z) * self.get_bias(z) * \
                (N_cen(np.exp(self.lnm), self.HOD) + \
                 N_sat(np.exp(self.lnm), self.HOD) * self.get_ug(z)[k_i, :])


    def I_2(self, z):
        """
        Integral 2-halo term in the galaxy power spectrum
        """
        ilist = []
        for i in range(len(self.universe.k)):
            ilist.append((1. / self.get_nbar(z)) * \
                          inty(self._I2inty(z, i), self.lnm))

        return np.array(ilist)


    def P_2h(self, z1, z2):
        """
        2-halo term in the galaxy power spectrum
        """
        # TODO : after finishing matter.py, modify this accordingly

        return self.get_Plin(z1, z2) * self.I_2(z1) * self.I_2(z2)


    def P_g(self, z1, z2):
        """
        Galaxy power spectrum
        Auto (Cross) power spectrum if z1=z2(z1!=z2)
        """
        if not (z1 in self.zvec) * (z2 in self.zvec):

            raise ValueError("These redshifts are not in the precomputed redshift values")

        return self.P_1h(z1, z2) + self.P_2h(z1, z2)
