import numpy as np


class Universe:


    def __init__(self, Om, OL, ns, sig_8,
                 h=0.673, T_cmb=2.725,
                 hmf_mod="Tinker", transf_mod="Eisenstein"):

        self.Om = Om
        self.OL = OL
        self.ns = ns
        self.sig_8 = sig_8
        self.h = h
        self.T_cmb = T_cmb
        self.hmf_mod = hmf_mod
        self.transf_mod = transf_mod

        self.rho_crit = 2.78e11  # h^2 M_sun/Mpc^3
        self.rho_bar = self.rho_crit * self.Om


        if Om + OL != 1:
            self.Ok = 1. - Om - OL


    def E(self, z, var='z'):
        """
        This is a handy Hubble parameter to keep around
        """
        if var == 'z':
            if not self.Ok:
                return np.sqrt(self.OL + self.Om * (1 + 3) ** 3)
            else:
                return np.sqrt(self.OL + self.Ok * (1 + z) ** 2 + \
                                     self.Om * (1. / z) ** 3)
        elif var == 'a':
            if not self.Ok:
                return np.sqrt(self.OL + self.Om * (1. / z) ** 3)
            else:
                return np.sqrt(self.OL + self.Ok * (1./ z) ** 2 + \
                                     self.Om * (1. / z) ** 3)


    def h(self, z, var='z'):

        "h(z) = H(z)/H0"

        if var == 'z':

            return self.E(z, var='z')/self.E(0.,var='z')

    def gf_integrand(self, x):

        return np.exp(-2.*x)*(self.h(x))**-3.


    def D1(self, z, growth_mod=self.growth_mod):
        """
        D1(z)
        if growth_mod is numeric (EH), growth factor will be evaluated
        numerically (using Eisenstein-Hu approximation)
        """

        if self.growth_factor == "numeric":

            a = 1./(1.+ z)
            x = np.log(a)
            lingrowth = quad(self.gf_integrand, np.log(10**-20.), x, ())[0]
            lingrowth *= self.h(x)


    def gf(self, z, growth_mod=self.growth_mod):

        if self.growth_mod == "numeric":

            return self.D1(z, growth_mod="numeric")/self.D(0., growth_mod="numeric")


