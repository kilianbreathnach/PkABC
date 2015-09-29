from simulator import Simul
from prior import Prior
from pmc_abc import PmcAbc



mydata = np.arange()

mysimul = Simul()

prior_dict = {}

myprior = Params(prior_dict)

my_abc = PmcAbc(mydata, mysimul, myprior)
