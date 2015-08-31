import numpy as np

def set_universe(Om=0.3, Ol=0.7, ns=0.96,
                 sig_8=0.82, h=0.673, T_cmb=2.725,
                 k_min=1.e-3, k_max=2.e3, dk=0.05,
                 lnM_min=np.log(1e11), lnM_max=np.log(1e15), dlnM=np.log(5e9),
                 transfer_fit="EH",
                 hmf_fit="Tinker",
                 bias_fit="Tinker")
