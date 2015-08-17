from numpy import np
from universe import Universe


class Pwspec(Universe):


    def _nbar_func(self, m, z):

        return self.hmf(m, z) * (self.N_cen(m) + self.N_sat(m))


    def nbar_g(self, z):

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
