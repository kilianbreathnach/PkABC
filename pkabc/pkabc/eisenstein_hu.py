import numpy as np


def transfnc_eh(k, Om=0.3, Om_b=0.0455, h=0.7, T_cmb=2.725, incl_baryons=True):
    """
    Compute the Eisenstein & Hu (1997) Transfer function for wave number k
    given a set of cosmological parameters.
    """

    if incl_baryons:
        """
        Shamelessly copying Tinker's implementation of the Eisenstein-Hu
        Transfer fitting function including baryonic effects.
        """

        # convert k to Mpc^-1 rather than hMpc^-1
        rk = k * h
        Omh2 = Om * h ** 2

        theta = T_cmb / 2.7

        # Eqn 4 - redshift of drag epoch
        b1 = 0.313 * (Omh2 ** (-0.419)) * (1 + 0.607 * (Omh2 ** (0.674)))
        b2 = 0.238 * (Omh2 ** 0.223)
        zd = 1291 * (Omh2 ** 0.251 / (1 + 0.659 * Omh2 ** 0.828)) * \
                    (1 + b1 * Omh2 **b2)

        # Equation 2 - redshift of matter-radiation equality
        z_eq = 2.5e4 * (Omh2 / theta ** 4)

        # value of R=(ratio of baryon-photon momentum density) at drag epoch
        R_drag = (31.5e3 * Om_b * h ** 2) / (theta ** 4 * zd)

        # value of R=(ratio of baryon-photon momentum density) at epoch of
        # matter-radiation equality
        R_eq = (31.5e3 * Om_b * h ** 2) / (theta ** 4 * z_eq)

        # Equation 3 - scale of particle horizon at matter-radiation equality
        k_eq = (7.46e-2 * Omh2) / (theta ** 2)

        # Equation 6 - sound horizon at drag epoch
        s = (2./ (3 * k_eq)) * np.sqrt(6. / R_eq) * \
            np.log((np.sqrt(1.+ R_drag) + np.sqrt(R_drag + R_eq)) / (1 + np.sqrt(R_eq)))

        # Equation 7 - silk damping scale
        k_silk = 1.6 * (Om_b * h ** 2) ** 0.52 * (Omh2 ** 0.73) * \
                       (1 + (10.4 * Omh2) ** (-0.95))

        # Equation 10  - define q
        q = rk / (13.41 * k_eq)

        # Equations 11 - CDM transfer function fits
        a1 = (46.9 * Omh2) ** 0.67 * (1 + (32.1 * Omh2) ** (-0.532))
        a2 = (12 * Omh2) ** 0.424 * (1 + (45 * Omh2) ** (-0.582))
        ac = a1 ** (- Om_b / Om) * a2 ** (- (Om_b/Om) ** 3)

        # Equations 12
        b1 = 0.944 / (1 + (458 * Omh2) ** (-0.708))
        b2 = (0.395 * Omh2) ** (-0.0266)
        Bc = 1. / (1 + b1 * ((1 - Om_b/Om) ** b2 - 1))

        # Equation 18
        f = 1. / (1 + ((rk * s) / 5.4) ** 4)

        # Equation 20
        c1 = 14.2 + 386. / (1 + 69.9 * q ** 1.08)
        c2 = 14.2 / ac + 386. / (1 + 69.9 * q ** 1.08)

        # Equation 17 - CDM transfer function with T_0 explicitly included
        T_c = f * np.log(np.exp(1) + 1.8 * Bc * q) / \
                  (np.log(np.exp(1) + 1.8 * Bc * q) + c1 * q * q) + \
           (1 - f) * np.log(np.exp(1) + 1.8 * Bc * q) / \
                     (np.log(np.exp(1) + 1.8 * Bc * q) + c2 * q * q)

        # Equation 15
        y = (1. + z_eq) / (1. + zd)
        G = y * (-6 * np.sqrt(1 + y) + \
                 (2 + 3 * y) * np.log((np.sqrt(1 + y) + 1) / \
                                      (np.sqrt(1 + y) - 1)))

        # Equation 14
        ab = G * 2.07 * k_eq * s * (1 + R_drag) ** (- 3. / 4)

        # Equation 23
        B_node = 8.41 * Omh2 ** 0.435

        # Equation 22
        stild = s * (1 + (B_node / (rk * s)) ** 3) ** (- 1. / 3)

        # Equation 24
        B_b = .5 + (Om_b / Om) + (3 - 2 * (Om_b / Om)) * \
                                 np.sqrt((17.2 * Omh2) ** 2 + 1)

        # Equations 19 & 21
        T_b = (np.log(np.exp(1) + 1.8 * q) / \
               (np.log(np.exp(1) + 1.8 * q)+ c1 * q * q)) / \
                (1 + ((rk * s) / 5.2) ** 2)

        T_b = (T_b + (ab * np.exp(- (rk / k_silk) ** 1.4)) / \
                     (1 + (B_b / (rk * s)) ** 3)) * \
               (np.sin(rk * stild) / (rk * stild))

        # Equation 8
        return (Om_b / Om) * T_b + (1 - Om_b / Om) * T_c

    elif not incl_baryons:
        """
        returns fitting formula for transfer
        function in the zero baryon limit
        """
        theta = T_cmb / 2.7
        gamma = Om*h

	# equation 28
	q = (k * (theta) ** 2) / (gamma*h) #equation 29
	Cq = 14.2 + 731. / (1. + 62.5 * q)
        Lq = np.log(2. * np.exp(1.) + 1.8 * q)

        return Lq / (Lq + Cq * q ** 2)
