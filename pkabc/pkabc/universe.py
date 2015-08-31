import numpy as np
from scipy.integrate import quad


class Universe(object):
    """
    Base class for Matter class. An instance of this class will hold the
    cosmological parameters for a LambdaCDM universe and functions to compute
    the Hubble parameter and growth of structure for such a version of the
    universe.
    """

    def __init__(self, Om=0.3, OL=0.7, ns=0.96,
                 sig_8=0.82, h=0.673, T_cmb=2.725):

        self.Om = Om
        self.OL = OL
        self.ns = ns
        self.sig_8 = sig_8
        self.h = h
        self.T_cmb = T_cmb

        if Om + OL != 1:
            self.Ok = 1. - Om - OL


    def Esq(self, z, var='z'):
        """
        The square of the dimensionless Hubble is often used
        """
        if var == 'z':
            if not self.Ok:
                return self.OL + self.Om * (1 + 3) ** 3
            else:
                return self.OL + self.Ok * (1 + z) ** 2 + \
                                 self.Om * (1. + z) ** 3
        elif var == 'a':
            if not self.Ok:
                return self.OL + self.Om * (1. / z) ** 3
            else:
                return self.OL + self.Ok * (1./ z) ** 2 + \
                                 self.Om * (1. / z) ** 3


    def E(self, z, var='z'):
        """
        Dimensionless Hubble parameter.
        """
        return np.sqrt(self.Esq(z, var=var))


    def gf_integrand(self, x):

        return np.exp(-2. * x)*(self.E(x, var='z'))**-3.


    def D1(self, z):
        """
        D1(z)
        if growth_mod is numeric (EH), growth factor will be evaluated
        numerically (using Eisenstein-Hu approximation)
        """
        if not isinstance(z, collections.Iterable):
            z = [z]
        d1 = np.zeros_like(z)

        if self.growth_factor == "numeric":

            a = 1./(1.+ z)
            x = np.log(a)
            for i, xx in enumerate(x):

              d1[i] = quad(self.gf_integrand, np.log(10**-20.), xx, ())[0]
              d1[i] *= self.h(x)

            return d1


    def gf(self, z):

        if self.growth_mod == "numeric":

            return self.D1(z, growth_mod="numeric")/self.D1(0., growth_mod="numeric")
