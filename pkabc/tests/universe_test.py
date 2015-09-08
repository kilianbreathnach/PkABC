import numpy as np
from pkabc.universe import Universe


testuniv = Universe()

zs = np.arange(0, 5, 0.05)

gf_vals = testuniv.gf(zs)

for i in range(len(zs)):
    print zs[i], gf_vals[i]
