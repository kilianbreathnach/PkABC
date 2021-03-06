import numpy as np
import collections
import scipy.integrate as intg
class Sigma(object):

    _defaults = {}

    def __init__(self, rho_mean, delta_c, k, p):

        self.rho_mean = rho_mean  #mean denisty of the universe
        self.delta_c = delta_c    #critical density of spherical collapse, default :1.68
        self.k = k                #wave number at which linear power spectrum is computed
        self.dk = k[1] - k[0]
        self.p = p                #linear power spectrum.


    def mass_to_radius(self, m):
        """
        The units of ``m`` should be consistent with
                ``rho_mean``.
        """
        return 0.6203504908994001 * (m / self.rho_mean) ** (1. / 3)


    def radius_to_mass(self, r):
        """
        The units of ``r`` should be consistent with
                ``rho_mean``.
        """

        return 4.1887902047863905 * r ** 3 * self.rho_mean


    def sigma_squared_r(self, r):
        """
        Calculate the squared mass variance : \sigma^{2}(r)
        """

        if not isinstance(r, collections.Iterable):
            r = [r]

        sigma = np.zeros_like(r)
        rest = self.p * self.k * self.k

        for i, rr in enumerate(r):
            integ = rest * self.tophat_k(rr * self.k) ** 2
            sigma[i] = (0.5 / np.pi ** 2) * intg.simps(integ, self.k)

        return sigma


    def sigma_squared_m(self, m):

        r = self.mass_to_radius(m)

        return self.sigma_squared_r(r)


    def dlnsigma_squared_dlnr(self, r):

        if not isinstance(r, collections.Iterable):
            r = [r]

        out = np.zeros_like(r)
        s = self.sigma_squared_r(r)
        rest = self.p * self.k ** 3.
        for i, rr in enumerate(r):
            y = rr * self.k
            w = self.tophat_k(y)
            dw = self.dwdx(y)
            integ = w * dw * rest
            out[i] = intg.simps(integ, self.k) / (np.pi ** 2 * s[i] ** 2)

        return out


    def dlnsigma_dlnm(self, m):
        """
        derivarive of the log of the mass variance
        (sigma not sigma^2!) w.r.t log mass.
        this is needed for calculating HMF.
        formulae:
        \frac{d\ln \sigma}{d\ln m}=.5*\frac{d\ln \sigma^2}{d\ln m}
                                  = .5* \frac{d\ln \sigma^2}{d\ln r} * \frac{d\ln r}{d \ln m}
                                  = \frac{1}{6} \frac{d\ln \sigma^2}{d\ln r}
        where \frac{d\ln r}{d \ln m} is equal to one third
        """
        r = self.mass_to_radius(m)

        return (1. / 6) * self.dlnsigma_squared_dlnr(r)


    def nu(self, m):
        """
        parameter \nu = \frac{\delta_c}{\sigma(r)}
        for standard parametrization of halo bias
        """
        r= self.mass_to_radius(m)
        return self.delta_c / np.sqrt(self.sigma_squared_r(r))

    def nu2(self, m):
        """
        parameter \nu^2 = (\frac{\delta_c}{\sigma(r)})^2
        for standard parametrization of halo bias
        """
        r = self.mass_to_radius(m)
        return self.delta_c**2. / self.sigma_squared_r(r)


    def tophat_k(self, kr):

        """
        vectorized form of top-hat window function in k-space,
        note the adhoc truncation of the analytic formula at small KR
        """
        W = np.ones(len(kr))
        KR = kr[kr>1.4e-6]
        W[kr>1.4e-6] = (3./KR**3.)*(np.sin(KR) - KR*np.cos(KR))

        return W

    def dwdx(self, kr):
        """
        dw(kr)/d(kr) : derivative of the tophat window
        function in kspace W(kr)
        """
        y = np.zeros(len(kr))
        KR = kr[kr>1.e-3]
        y[kr>1.e-3] = (9 * KR * np.cos(KR) + 3.*(KR**2.- 3.)*np.sin(KR))/(KR**4.)

        return y


    def normalize_power(self, sigma_8):
        """
        input: cosmological parameter sigma_8
        returns normalization of nunormalized powerspectrum(P0(k)=k^ns*T(k)^2)
        math : (\frac{\hat{\sigma}_8}{\sigma_8})^2
        """
        sigma_8_hat_squared = self.sigma_squared_r(8.0)
        norm = sigma_8**2./sigma_8_hat_squared

        return norm
