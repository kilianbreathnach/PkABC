import numpy as np
cimport numpy as np


cdef double DELTA_HALO = 200
cdef double RHO_CRIT = 2.775e11
cdef double OMEGA_M = 0.31
cdef double REDSHIFT = 0.515


def density(double r, double Rs):

    return Rs / (r * (1 + r / Rs) * (1 + r / Rs))


def galpos(double M):

    cdef np.ndarray[np.float64_t, ndim=1] pos = np.zeros(3, dtype=np.float64)
    cdef double cvir, Rvir, Rs, rho_max, r, rho_r
    cdef double costheta, sintheta, signs, phi1
    cvir = 12 * (M / 10 ** 12) ** (-0.11)
    Rvir = (3 * (M / (4 * DELTA_HALO * np.pi * RHO_CRIT * OMEGA_M))) ** (1./3)
    Rs = Rvir / cvir
    rho_max = density(Rvir, Rs) * Rvir * Rvir * 4.0 * np.pi

    while True:

        r = np.random.rand() * Rvir
        rho_r = density(r, Rs) * r * r * 4.0 * (np.pi / rho_max)

        if (np.random.rand() <= rho_r):

            costheta = 2 * (np.random.rand() - 0.5)
            sintheta = np.sqrt(1 - costheta * costheta)
            signs = 2 * (np.random.rand() - 0.5)
            costheta = signs * costheta / np.abs(signs)
            phi1 = 2 * np.pi * np.random.rand()

            pos[0] = r * sintheta * np.cos(phi1)
            pos[1] = r * sintheta * np.sin(phi1)
            pos[2] = r * costheta

            return pos


def galvel(double M):

    cdef np.ndarray[np.float64_t, ndim=1] vel = np.zeros(3, dtype=np.float64)

    cdef double fac, sigv

    fac = np.sqrt(4.499e-48) * (4 * DELTA_HALO * np.pi * OMEGA_M * \
          RHO_CRIT / 3.) ** (1. / 6) * 3.09e19 * np.sqrt(1+REDSHIFT)

    sigv = fac * (M ** (1./3) / np.sqrt(2.0))

    vel = np.random.normal(size=3) * sigv

    return vel
