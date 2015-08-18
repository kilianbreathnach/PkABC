from numpy import np
from universe import Universe
from scipy.integrate import quad as inty
from scipy.fftpack import
from eisenstein_hu import transfnc_eh


class Matter(Universe):

    def Transf(self, k, transf_mod=self.transf_mod):
        """
        transfer function to add the effects of gravitational instabilities
        and/or baryonic interactions to the fluctuations generated
        during inflation
        """

        if self.transf_mod == "Eisenstein":
            transfnc_eh(k, )


    def window(kr):

        return 3 * ((np.sin(kr) = kr * np.cos(kr)) / kr ** 3)


    def powsp(self, k):
        """
        your basic power spectrum
        """

        return k ** self.ns * Transf(k) ** 2


    def _D1int(self, a, Om, OL):




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
