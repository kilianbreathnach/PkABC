import numpy as np
from universe import Universe
from scipy.integrate import quad as inty
#from scipy.fftpack import
#from eisenstein_hu import transfnc_eh
"""sanity check
xniverse = Universe(Om=.3, OL=.7, ns=1., sig_8=.8 ,h=0.673,\
                    T_cmb=2.725 , hmf_mod="Tinker", transf_mod="Eisenstein")
print xniverse.Om
"""
class Matter(Universe):

    

    def Transf(self, k, transf_mod="Eisenstein"):

        """
        transfer function to add the effects of gravitational instabilities
        and/or baryonic interactions to the fluctuations generated
        during inflation
        """

        if self.transf_mod == "Eisenstein":
            transfnc_eh(k, )


    def tophat_kspace(self, kr):
        """
        tophat window function in k space
        note: there could be numerical issues near kr=0,
        this should probably be truncated around kr~10^-6ish
        """

        return 3 * ((np.sin(kr) - kr * np.cos(kr)) / kr ** 3)


    def vec_tophat_kspace(self, kr):

        """
        vectorized form of top-hat window function in k-space,
        note the truncation of the analytic formula at small KR
        """
        W = np.ones(len(kr))
        KR = kr[kr>1.4e-6]
        W[kr>1.4e-6] = (3./KR**3.)*(np.sin(KR) - KR*np.cos(KR))

        return W 

    def unnormal_powsp(self, k):
        """
        unnormalized power spectrum P0(k)=k^nsT(k)
        """
        return k ** self.ns * Transf(k) ** 2



    def _D1int(self, a, Om, OL):

        return (a * self.E(a, var='a')) ** (-3)


    def D1(self, z):

        a = 1. / (1 + z)

        return self.E(z) * inty(self._D1int, 0.5, a, args=(self.Om, self.OL))[0]


    def growthfac(self, z):

        return (self.D1() / self.D1()) ** 2


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


    def _sigint(self, k, R):

        kr = k * R

        return self.Plin(k) * self.window(kr) * k ** 2


    def sig(self, R):

        return np.sqrt(inty(self._sigint, self.klo, self.khi, args=(R,))[0])


    def rho_crit(z):
        """
        Compute critical density of the universe at a given redshift:
        p_c = 3 * H^2(z) / (8 * pi * G)
        and convert to units of solar mass per cubic megaparsec
        """
        Hsq = (100 * self.h * self.E(z)) ** 2
        rhoc = 27746.677665951076 * Hsq  # using G = 4.302e-6 Mpc Msun^-1 (km/s)^2
        return


    def rho_bar(z):
        pass
