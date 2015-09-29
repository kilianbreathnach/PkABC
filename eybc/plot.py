import numpy as np
import corner
import matplotlib.pyplot as plt
import seaborn as sns
plt.switch_backend('agg')

T = 20

import matplotlib.cm as cm


for i in xrange(T):

   w = np.loadtxt("theta_w_t"+str(i)+".dat")
   print w.shape
   #print w.shape
   theta = w[:,0:2]
   weights = w[:,2:]

   sns.jointplot(x=theta[:,0], y=theta[:,1], kind="kde" , style="white" , weights = weights , xlim = [11 , 13.5]  ,ylim = [-.2 , 1.])
   print w.shape

   #sns.jointplot(x=theta[:,0], y=theta[:,1], 
   #               kind="kde", xlim=[11.,14.], ylim=[0.1,0.3],
   #               joint_kws={"gridsize":30, "cmap":cm.BuPu,
   #                     "extent":[11.,14., 0.1, 0.3]}, 
   #               marginal_kws={"bins":20, "hist":{"range":[11.,14., .1,.3]}})
   """
   figure = corner.corner(theta, range=[[11,14],[0.1 , 0.3]], weights=weights, color="k",
           smooth=False,smooth1d = False,
           labels=None, label_kwargs=None,
           show_titles=False, title_fmt=".2f", title_kwargs=None,
           truths=[12.02 , 0.26], truth_color="#4682b4")"""
   plt.savefig("scatter"+str(i)+".png")
