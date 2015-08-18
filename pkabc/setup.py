from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name = 'vpf counting',
    ext_modules = cythonize(["vpf.pyx", "NFW.pyx", "HOD.pyx"]),
    include_dirs=[np.get_include()]
)
