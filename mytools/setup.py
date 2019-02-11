# Setup file for cython
import sys
from distutils.core import setup
from Cython.Build import cythonize
'''
Compile the cython modules. Corresponding module shall be passed as argument.

'''
module_name = "met_tools"
setup(
    name = module_name,
    ext_modules = cythonize(module_name+".pyx"),
)



