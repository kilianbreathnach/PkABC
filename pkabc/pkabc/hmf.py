import numpy as np
from Dmstats import Matter
from scipy.integrate import quad as inty

class HMF(Matter):


    def hmf(self, m, z):
        """
        Compute the halo mass function for a variety of models
        given the cosmological parameters
        """

        R = 0.6203504908994001 * (m / self.rho_bar(z)) ** (1. / 3)

        sig = self.sig(R)

        if self.hmf_mod == "Tinker":

            f = 0.186 * ((self.sig(R) / self.b_t) ** (- self.a_t) + 1) \
                         * np.exp(- self.c_t / self.sig(R) ** 2)

            n = f * (self.rho_bar(z) / m) * self.dlnsiginvdM

        return n

