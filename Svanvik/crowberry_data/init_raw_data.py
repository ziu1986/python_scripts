import os, sys, glob
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D # Register 3d projection
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.ozone_tools import flunder

# Data source
src = os.environ['DATA']+'/astra_data/observations/krekling_svanvik/raadata_felt/'

try:
    data_list
except NameError:
    data_list = []
    exp_list = []
    for file in sorted(glob.glob(src+'*2019*')):
        print file
        data_list.append(pd.read_csv(file, header=11, delimiter='\\t'))
        print 
        exp_list.append(os.path.basename(file)[:10])

execfile("fits_raw_data.py")
execfile("plot_raw_data.py")    

# Units are explained in handbook section 16
# [Pressure] = kPa
# [Temperature] = degC
# [Tleaf] = degC
# [H2O] = mmol/mol
# [CO2] = mumol/mol
# [Flux] = mumol/s
# [PAR] = mumol/m2/s
# [VPD] = kPa
# [gsto] = mol H2O/m2/s
# [An] = mol CO2/m2/s
