import numpy as np
from scipy.optimize import minimize
from scipy.special import sici
from universe import Universe
from sigma import Sigma
from eisenstein_hu import transfnc_eh

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
                 k=None,
                 k_min=np.exp(-5), k_max=np.exp(5), kbins=1000,
                 lnM=None,
                 lnM_min=np.log(1e11), lnM_max=np.log(1e15),
                 lnMbins=100,
                 transfer_fit="EH",
                 hmf_fit="Tinker",
                 bias_fit="Tinker",
                 incl_baryons=True,
                 **kwargs):

        #initializing cosmology inherited from the Universe class
        super(Matter, self).__init__(**kwargs)

        #initializing other stuff
        if k != None:
            self.k = k
        else:
            self.k_min = k_min
            self.k_max = k_max
       	    self.kbins = kbins
            self.gen_k()

        if lnM != None:
            self.lnM = lnM
        else:
            self.lnM_min = lnM_min
            self.lnM_max = lnM_max
            binfrac = 1. / lnMbins
            self.dlnM = binfrac * (lnM_max - lnM_min)
            self.gen_lnM()

   	self.transfer_fit = transfer_fit
        self.hmf_fit = hmf_fit
        self.bias_fit = bias_fit
        self.incl_baryons = incl_baryons

        self.trackz = -1

        self.master_sig = Sigma(1, 1.686, self.k, self.unnormal_p0_lin())
        self.norm_P0 = self.master_sig.normalize_power(self.sig_8)


    def gen_k(self):
        """
        array of wave numbers
        """
        self.k = np.logspace(np.log(self.k_min), np.log(self.k_max),
                             num=self.kbins)


    def gen_lnM(self):
        """
        array of log halo masses
        """
        self.lnM = np.arange(self.lnM_min, self.lnM_max, self.dlnM)


    def T(self):
        """
        transfer function
        """
        return transfnc_eh(self.k, Om=self.Om, h=self.h, T_cmb=self.T_cmb,
                           incl_baryons=self.incl_baryons)


    def delta_c(self, z):
        """
        critical matter over density for spherical
        collapse at redshift z, taken from van der Bosch (2013)
        """
 	return 1.686470199841145 * (self.omegamz(z) ** 0.0055)


    def unnormal_p0_lin(self):
        """
        basic power spectrum
        """
        return (self.k ** self.ns) * (self.T() ** 2)


    def normal_p0_lin(self):
        """
        normalized linear power spectrum at redshift zero
        """
        return  self.norm_P0 * self.unnormal_p0_lin()


    def normal_pz_lin(self, z):
        """
        linear power spectrum at redshift z
        """
        return self.normal_p0_lin() * self.gf(z)


    def reset_mastersig(self, z):

        self.master_sig = Sigma(self.rho_mean(z), self.delta_c(z), self.k,
                                self.normal_pz_lin(z))


    def check_z(self, z):

        if z != self.trackz:
            self.reset_mastersig(z)
            self.trackz = z

    def _hmf_fit(self, z):
        """
        Tinker, 2008 fitting function for halo mass function.
        """
        A_0 = 0.186
        a_0 = 1.47
        b_0 = 2.57
        c_0 = 1.19

        alph = np.exp( - (0.75 / np.log(200. / 75)) ** 1.2)

        A_z = A_0 * (1 + z) ** (-0.14)
        a_z = a_0 * (1 + z) ** (-0.06)
        b_z = b_0 * (1 + z) ** ( - alph)

        self.check_z(z)
        sig = np.sqrt(self.master_sig.sigma_squared_m(np.exp(self.lnM)))

        return A_z * ((sig / b_z) ** (- a_z) + 1) * np.exp(- c_0 / sig ** 2)


    def hmf(self, z):
        """
        Compute the halo mass function for a variety of models
        given the cosmological parameters
        """
        if self.hmf_fit  == "Tinker":

            f = self._hmf_fit(z)

            n = - f * (self.rho_mean(z) / np.exp(self.lnM) ** 2) * \
                      self.master_sig.dlnsigma_dlnm(np.exp(self.lnM))

        return n


    def c_minfunc(self, z, Rstar):

        s = self.sigma(Rstar)

        return np.abs(s - self.delta_c(z))


    def c200(self, m, z):
        """
        returns halo concentration
        """

        # From v. d. Bosch, Miyatake...
        F = 0.001
        K = 2.7

        Mstar = F * m
        Rstar = self.master_sig.mass_to_radius(Mstar)

        zc = minimize(self.c_minfunc, 0.5, args=(Rstar,))['x']

        return K * ((1 + zc) / (1 + z))


    def u_g(self, z):
        """
        returns a matrix of u_g(k|m, z) in k and m for each redshift
        """
        self.check_z(z)

        umat = np.zeros((self.k.shape[0] , self.lnM.shape[0]))

        for i, m in enumerate(np.exp(self.lnM)):

            c = c200(m, z)
            d200 = (200. / 3) * (c ** 3 / (np.log(1 + c) - c / (1 + c)))
            mu = self.k * (self.master_sig.mass_to_radius(np.exp(self.lnM)) / c)   # mu is a k vector
            umat[:, i] = ((3 * d200) / (200 * c ** 3)) * \
                         (np.cos(mu) * (sici(mu + mu * c)[1] - sici(mu)[1]) + \
                          np.sin(mu) * (sici(mu + mu * c)[0] - sici(mu)[0]) - \
                          np.sin(mu * c) / (mu + mu * c) )

        return umat


    def bias(self, z):
        """
        returns a 1d array of halo bias b_h(m,z) in m for each redshift
        (Tinker 2010 model)
        """
	m = np.exp(self.lnM)

        self.check_z(z)

        R = R_m(m, z)
        s = self.sigma(R)
        nu = self.delta_c(z) / s

        # parameters for halo overdensity of 200
        y = np.log10(200)
        ey4 = np.exp(- (4. / y) ** 4)
        A = 1 + 0.24 * y * ey4
        a = 0.44 * y - 0.88
        B = 0.183
        b = 1.5
        C = 0.019 + 0.107 * y + 0.19 * ey4
        c = 2.4

        bias = 1. - (nu ** a)/(nu ** a + self.delta_c(z) ** a) + \
	             B * nu ** b + C * nu ** c

        return bias
