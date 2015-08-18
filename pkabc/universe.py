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
