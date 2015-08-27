import numpy as np
from scipy.optimize import minimize
from universe import Universe
from sigma import Sigma
import functions
import halo_functions
from eisenstein_hu import transferfnc_eh

class Matter(Universe):

   """
   The purpose of this class is to deliver all the
   objects necessary to compute the galaxy-galaxy
   or galaxy-matter power spectra.
   These are:
   ------------------------------------------------
   Inputs :

   Cosmology = list of cosmological parameters
               inherited from the Universe
   k = array of wave numbers (unit = Mpc^-1)
   z = redshift
   M = list of halo masses (unit Msun/h)
   -------------------------------------------------

   Returns:

   Plin(k, z) = 1d Array: linear matter power spectrum
                evaluated at array of k values
                at redshift z,

   dndlnm(z)  = 1d Array: \frac{d\n}{d\ln(m)}
                Halo mass function evaluated at array of
                halo masses (M) and redshift z

   bh(m, z)   = 1d Array: halo bias evaluated at array of halo masses
                (M) and redshift z

   ug(k, m, z)= 2d Array: Fourier transform of dark matter
                halo density profile (here we assume NFW)

   ------------------------------------------------
   """
   def __init__(self,
                 k_min = 1.e-3, k_max = 2.e3, dk = 0.05,
                 lnM_min=np.log(1.e4), lnM_max=np.log(5.e4), dlnM = np.log(1.e1)
                 transfer_fit = "EH",
                 hmf_fit = "Tinker",
                 bias_fit = "Tinker",
                 **kwargs):

           #initializing cosmology inherited from the Universe class
           super(Matter, self).__init__(**kwargs)

           #initializing other stuff
           self.z  = z
           self.k_min = k_min
           self.k_max = k_max
       	   self.dk = dk
           self.lnM_min = lnM_min
           self.lnM_max = lnM_max
       	   self.dlnM = dlnM
   	   self.transfer_fit = transfer_fit
           self.hmf_fit = hmf_fit
           self.bias_fit = bias_fit

           self.trackz = -1


   def omegamz(self, z):
       """
       dimensionless matter density at redshift z
       """
           return self.omegam/(self.omegam + self.omegal*(1.+ z)**-3.)

    def rho_crit(self, z):

        return 2.77e11 * self.h * self.Esq(z)


    def rho_mean(self, z):

        return self.omegamz(z) * self.rho_crit(z)


   def delta_c(self, z):
       """
       critical matter over density for spherical
       collapse at redshift z
       """

	   return 0.15 * ((2.*np.pi)**(2./3))* (self.omegamz(z))**0.0055

   def k(self):
       """
       array of wave numbers
       """

           return np.arange(self.k_min, self.k_max, self.dk)

   def lnM(self):
       """
       array of log halo masses
       """

           return np.arange(self.lnM_min, self.lnM_max, self.lnM)

   def T(self):
       """
       transfer function
       """

           return transfnc_eh(self.k)

   def unnormal_p0_lin(self):
       """
       basic power spectrum
       """

       	   return self.k**self.n * self.T

   def normal_p0_lin(self):
       """
       normalized linear power spectrum at redshift zero
       """
           # TODO : make sure rho_mean is delivered by the Universe

           sig = Sigma(self.rho_mean, self.delta_c, self.k, self.unnormal_p0_lin)
	   norm = sig.normalize_power(self.sigma_8)

           return  norm * self.unnormal_p0_lin

   def normal_pz_lin(self, z):
       """
       linear power spectrum at redshift z
       """
           return self.normal_p0_lin * self.gf(z)


    def reset_mastersig(self, z):
        self.master_sig = Sigma(self.rho_mean(z), self.delta_c(z), self.k,
                     self.normal_pz_lin(z))


    def sigma(self, R):
        return mp.sqrt(self.master_sig.sigma_squared_r(R))


    def check_z(self, z):

        if z != self.trackz:
            reset_mastersig(z)
            self.trackz = z


    def hmf(self, lnm, z):
        """
        Compute the halo mass function for a variety of models
        given the cosmological parameters
        """

        A_0 = 0.186
        a_0 = 1.47
        b_0 = 2.57
        c_0 = 1.19

        alph = np.exp( - (0.75 / np.log(200. / 75)) ** 1.2)

        A_z = A_0 * (1 + z) ** (-0.14)
        a_z = a_0 * (1 + z) ** (-0.06)
        b_z = b_0 * (1 + z) ** ( - alph)

        R = 0.6203504908994001 * (np.exp(lnm) / self.rho_mean(z)) ** (1. / 3)

        self.check_z(z)

        sig = self.sigma(R)

        if self.hmf_mod == "Tinker":

            f = A_z * ((sig / b_z) ** (- a_z) + 1) * np.exp(- c_0 / sig ** 2)

            n = - f * (self.rho_mean(z) / np.exp(lnm) ** 2) * \
                      self.master_sig.dlnsigma_dlnm(np.exp(lnm))

        return n


    def c_minfunc(self, z, Rstar):

        self.check_z(z)

        s = self.sigma(Rstar)

        return np.abs(s - delta_c(z))


    def c200(self, m, z):
        """
        returns halo concentration
        """

        # From v. d. Bosch, Miyatake...
        F = 0.001
        K = 2.7

        Mstar = F * m
        Rstar = 0.6203504908994001 * (Mstar / self.rho_mean(z)) ** (1. / 3)

        zc = minimize(self.c_minfunc, 0.5, args=(Rstar,))['x']

        return K * ((1 + zc) / (1 + z))


    def u_g(self, k, lnm, z):
        """
        returns a matrix of u_g(k|m, z) in k and m for each redshift
        """

        umat = np.zeros((k.shape[0], lnm.shape[0]))


        umat =





