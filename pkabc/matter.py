import numpy as np
from universe import Universe
import DMstats
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
   def __init__(self, z=0.0, 
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

   def growth(self):

           return Universe.gf(self.z)

   def omegamz(self):

           return self.omegam/(self.omegam + self.omegal*(1.+self.z)**-3.)

   def delta_c(self):

	   return 0.15 * ((2.*np.pi)**(2./3))* (self.omegamz)**0.0055

   def k(self):
          
           return np.arange(self.k_min, self.k_max, self.dk)

   def lnM(self):

           return np.arange(self.lnM_min, self.lnM_max, self.lnM)

   def T(self):

           return transfnc_eh(self.k)
           
   def unnormal_p_lin(self):
      
       	   return self.k**self.n * self.T

   def normal_p_lin(self):
	   
           # TODO : make sure rho_mean is delivered by the Universe

           sig = Sigma(self.rho_mean, self.delta_c, self.k, self.unnormal_p_lin)
	   norm = sig.normalize_power(self.sigma_8)
          
           return  norm * self.unnormal_p_lin	

   def p_lin(self):

           return self.normal_p_lin * self.growth 
