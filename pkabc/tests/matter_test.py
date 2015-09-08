import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from matter import Matter


testmat = Matter()

zs = np.arange(0, 2, 0.1)
cmap = plt.get_cmap("summer")
cNorm = colors.Normalize(vmin=np.min(zs), vmax=np.max(zs))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

# # test transfer function
# testmat.incl_baryons = False
# transfnc_nobaryons = testmat.T()
# testmat.incl_baryons = True
# transfnc = testmat.T()
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.loglog(testmat.k, transfnc, 'b', testmat.k, transfnc_nobaryons, 'r')
# fig.savefig("transfncs.pdf")
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.loglog(testmat.k, transfnc - transfnc_nobaryons)
# fig.savefig("transfnc_diff.pdf")
#
# # test power spectrum normalisation
# unnorm = testmat.unnormal_p0_lin()
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.loglog(testmat.k, unnorm)
# fig.savefig("unnormed_p.pdf")
#
# # test normalised against unnormalised
# normed = testmat.normal_p0_lin()
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.loglog(testmat.k, unnorm, testmat.k, normed)
# fig.savefig("normed_p.pdf")
#
# # test power specturm as function of z
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ls = []
#
# for z in zs:
#     p_z = testmat.normal_pz_lin(z)
#     colorVal = scalarMap.to_rgba(z)
#     ls.append(ax.loglog(testmat.k, p_z, color=colorVal, alpha=0.5))
#
# cax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
# cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=cNorm,
#                                spacing='proportional')
#
# fig.savefig("normed_p_z.pdf")
#
# # testing Tinker's fitting function
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ls = []
#
# for z in zs:
#     fit = testmat._hmf_fit(z)
#     sig = np.sqrt(testmat.master_sig.sigma_squared_m(np.exp(testmat.lnM)))
#     loginvsig = np.log10(1. / sig)
#     colorVal = scalarMap.to_rgba(z)
#     ls.append(ax.plot(loginvsig, fit, color=colorVal, alpha=0.5))
#
# cax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
# cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=cNorm,
#                                spacing='proportional')
# fig.savefig("Tinker_fit_invsig.pdf")
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ls = []
#
# for z in zs:
#     fit = testmat._hmf_fit(z)
#     colorVal = scalarMap.to_rgba(z)
#     ls.append(ax.semilogx(np.exp(testmat.lnM), fit, color=colorVal, alpha=0.5))
#
# cax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
# cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=cNorm,
#                                spacing='proportional')
# fig.savefig("Tinker_fit_mass.pdf")
#
# # testing hmf
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ls = []
#
# for z in zs:
#     hmf = testmat.hmf(z)
#     colorVal = scalarMap.to_rgba(z)
#     ls.append(ax.loglog(np.exp(testmat.lnM), hmf, color=colorVal, alpha=0.5))
#
# cax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
# cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=cNorm,
#                                spacing='proportional')
# fig.savefig("hmf.pdf")
#
# # testing bias
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ls = []
#
# for z in zs:
#     bias = testmat.bias(z)
#     nu = testmat.delta_c(z) / np.sqrt(testmat.master_sig.sigma_squared_m(np.exp(testmat.lnM)))
#     colorVal = scalarMap.to_rgba(z)
#     ls.append(ax.plot(np.log10(nu), np.log10(bias), color=colorVal, alpha=0.5))
#
# cax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
# cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=cNorm,
#                                spacing='proportional')
# fig.savefig("bias.pdf")

# testing NFW power spectrum
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ls = []

cNorm_mass = colors.Normalize(vmin=0, vmax=len(testmat.lnM))
scalarMap_mass = cmx.ScalarMappable(norm=cNorm_mass, cmap=cmap)

for z in zs:

    nfwpsp = testmat.u_g(z)

    for i, m in enumerate(reversed(testmat.lnM)):

        u_g = nfwpsp[:, i]
        colorVal = scalarMap_mass.to_rgba(i)
        ls.append(ax.plot(np.log10(testmat.k), np.log(u_g), zs=m, color=colorVal, alpha=0.5))
#        time.sleep(3)
    break

cax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=cNorm_mass,
                               spacing='proportional')
fig.savefig("u_g_0.pdf")



