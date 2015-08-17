import numpy as np
from scipy.integrate import quad as inty
from scipy.special import erf
from scipy.fftpack import


class Universe:


    def __init__(self, Om, OL, ns, sig_8,
                 hmf_mod="Tinker", transf_mod="Eisenstein"):

        self.Om = Om
        self.OL = OL
        self.ns = ns
        self.sig_8 = sig_8
        self.hmf_mod = hmf_mod
        self.transf_mod = transf_mod

        self.Ok = 1. - Om - OL


    def E(self, z):
        """
        This is a handy Hubble parameter to keep around
        """
        return np.sqrt(self.OL + self.Ok * (1 + z) ** 2 + self.Om * (1 + 3) ** 3)


    def Transf(self, k):
        """
        transfer function to add the effects of gravitational instabilities
        and/or baryonic interactions to the fluctuations generated
        during inflation
        """

        if self.transf_mod == "Eisenstein":
            pass


    def window(kr):

        return 3 * ((np.sin(kr) = kr * np.cos(kr)) / kr ** 3)


    def powsp(self, k):
        """
        your basic power spectrum
        """

        return k ** self.ns * Transf(k) ** 2


    def Plin(self, k):
        """
        function for the linear matter power spectrum
        """

        if not self.pnorm:
            self.get_pnorm()

        return 2 * (np.pi ** 2) * (self.delta_H ** 2) * \
                   (k ** ns) * ((1 / self.H_0) ** (ns + 3)) * (self.Transf(k) ** 2)


    def get_pnorm(self):
        pass


    def _sigint(k, R):

        kr = k * R

        return self.Plin(k) * self.window(kr) * k ** 2


    def sig(self, R):

        return np.sqrt(inty(self._sigint, self.klo, self.khi, args=(R,))[0])

    def rho_bar(z):
        pass

    def hmf(self, m, z):
        """
        Compute the halo mass function for a variety of models
        given the cosmological parameters
        """

        R = 0.6203504908994001 * (m / self.rho_bar(z)) ** (1. / 3)

        if self.hmf_mod == "Tinker":

            f = self.A_t * ((self.sig(R) / self.b_t) ** (- self.a_t) + 1) \
                         * np.exp(- self.c_t / self.sig(R) ** 2)

            n = f * (self.rho_bar(z) / m) * self.dlnsiginvdM

        return n


    def N_cen(self, m):

        return 0.5 * (1 + erf((np.log10(m) - self.logMmin) / self.sig_logm))

    def N_sat(self, m):

        return N_cen * (m / self.m_1) ** self.alph_sat * np.exp(- self.m_cut / m)


    def _nbar_func(self, m, z):

        return self.hmf(m, z) * (self.N_cen(m) + self.N_sat(m))


    def nbar_g(self, z):

        return inty(self._nbar_func, 0, np.inf, args=(z,))[0]





class Pwspec(Universe):


    def u_g(self, k, z, m):




    def P_1h(self, k, z1, z2):

        if z1 != z2:
            return 0.
        else:
            return (1. / self.nbar_g(z1) ** 2) * \
                    inty(self.hmf, 0, np.inf, args=(z,))[0] * \
                    (self.N_sat(m) ** 2 * self.u_g(k, z1, m) ** 2 + \
                    2 * self.N_cen(m) * self.N_sat(m) * self.u_g(k, z1, m))


    def _I2inty(self, m, k, z):


    def I_2(k, z):

        return (1. / self.nbar_g(z1)) * \
                inty(self._I2inty, )


    def P_2h(self, k, z1, z2):


    def P_g(self, k, z1, z2):

        return P_1h(k, z1, z2) + P_2h(k, z1, z2)
