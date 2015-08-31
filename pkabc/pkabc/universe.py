import numpy as np
from scipy.integrate import quad
import collections

class Universe(object):
    """
    Base class for Matter class. An instance of this class will hold the
    cosmological parameters for a LambdaCDM universe and functions to compute
    the Hubble parameter and growth of structure for such a version of the
    universe.
    """


    def __init__(self, Om=0.3, OL=0.7, Ok= None, ns=0.96,
                 sig_8=0.82, h=0.673, T_cmb=2.725, growth_mod = "numeric"):

        self.Om = Om
        self.OL = OL
        self.ns = ns
        self.sig_8 = sig_8
        self.h = h
        self.T_cmb = T_cmb
        self.Ok = Ok
        if Om + OL != 1:
            self.Ok = 1. - Om - OL
        self.growth_mod = growth_mod


    def Esq(self, z, var='z'):
        """
        The square of the dimensionless Hubble is often used
        """
        if var == 'z':
            if not self.Ok:
                return self.OL + self.Om * (1. + z) ** 3
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


    def omegamz(self, z):
        """
        dimensionless matter density at redshift z
        """
        return (self.Om * (1 + z) ** 3) / self.Esq(z)


    def rho_crit(self, z):

        return 2.77e11 * self.h * self.Esq(z)


    def rho_mean(self, z):

        return self.omegamz(z) * self.rho_crit(z)

    
    def ha(self, a):
  
        "H(a)/H0"

        return (self.Om * (1. / a) ** 3. + self.OL) ** -0.5


    def gf_integrand(self, a):

        return (a * self.ha(a)) ** -3.


    def D1(self, z):
        """
        D1(z)
        if growth_mod is numeric (EH), growth factor will be evaluated
        numerically (using Eisenstein-Hu approximation)
        """

        d1 = np.zeros_like(z)

        if self.growth_mod == "numeric":

            a = 1./(1.+ z)
            for i, xx in enumerate(a):

              d1[i] = quad(self.gf_integrand, 10.**-4., xx, ())[0]
              d1[i] *= self.ha(a)

            return d1


    def gf(self, z):

        if self.growth_mod == "numeric":

            return self.D1(z)/self.D1(0.)
