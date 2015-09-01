import numpy as np
import matplotlib.pyplot as plt
from pkabc.matter import Matter


testmat = Matter()

transfnc = testmat.T()


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(testmat.k, transfnc)
fig.show()
