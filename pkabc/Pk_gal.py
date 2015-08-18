import numpy as np
from scipy.integrate import quad as inty
from hod import N_cen, N_sat
from DMstats import Matter


class Pwspec(Matter):


    def _nbar_func(self, m, z):

        return self.hmf(m, z) * \
                (N_cen(m, self.logMmin, self.sig_logm) + \
                 N_sat(m, self.logMmin, self.sig_logm,
                       self.m_1, self.alph_sat, self.m_cut))


    def nbar_g(self, z):
        #TODO: figure out matter limits on integral here
        return inty(self._nbar_func, 0, np.inf, args=(z,))[0]


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
        pass


    def I_2(k, z):

        return (1. / self.nbar_g(z1)) * \
                inty(self._I2inty, )


    def P_2h(self, k, z1, z2):
        pass


    def P_g(self, k, z1, z2):

        return P_1h(k, z1, z2) + P_2h(k, z1, z2)
