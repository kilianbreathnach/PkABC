import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from Pk_gal import Pwspec

zs = np.array([0.52])

# miyatake
hod = np.array([13.15, np.sqrt(0.16), 1.26, 14.16, np.log(0.71) + 13.15])

hods = [np.array([13.15, np.sqrt(0.16), 1.26, 14.16, np.log(0.71) + 13.15]),
        np.array([13.25, np.sqrt(0.17), 1.36, 14.56, np.log(0.71) + 13.35]),
        np.array([12.15, np.sqrt(0.15), 1.16, 13.56, np.log(0.81) + 13.35]),
        np.array([12.35, np.sqrt(0.17), 1.56, 13.56, np.log(0.61) + 12.95]),
        np.array([12.85, np.sqrt(0.15), 1.27, 14.16, np.log(0.91) + 13.45])]


psp = Pwspec(zs)
psp.precompute()

# get changs psp
chang = np.loadtxt("../observed/power_cmass-dr12v4-N-Reid.dat.fiducial.grid960.P020000.box3600", usecols=(0, 1))

fig = plt.figure()
ax = fig.add_subplot(111)
ls = []

mask = np.where((psp.universe.k >= np.min(chang[:, 0])) * (psp.universe.k <= np.max(chang[:, 0])))

cmap = plt.get_cmap("summer")
cNorm = colors.Normalize(vmin=0, vmax=len(hods))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

print "Done with precomputing bitch"

for i, hod in enumerate(hods):

    psp.set_HOD(hod)
    pg = psp.P_g(zs[0], zs[0])
    colorVal = scalarMap.to_rgba(i)
    ls.append(ax.loglog(psp.universe.k[mask], pg[mask], alpha=0.5, color=colorVal, lw=2))
    print "hod done"

ls.append(ax.loglog(chang[:, 0], chang[:, 1], alpha=0.5, color='r', lw=5))

fig.savefig("miyatake_chang2.pdf")
