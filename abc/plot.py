import numpy as np
import triangle
import matplotlib.pyplot as plt

plt.switch_backend('agg')

T = 17

for i in xrange(T):

   w = np.loadtxt("theta_w_t"+str(i)+".dat")
   print w.shape
   data = w[:,0:2]
   weights = w[:,2:]
   figure = triangle.corner(data, labels=[r"$\mu$", r"$\sigma$"], range = [[-2,2],[0,1]],

                         show_titles=True, title_args={"fontsize": 12})
   figure.gca().annotate(str(i), xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")
   figure.savefig("scatter"+str(i)+".png")
~                                                          
