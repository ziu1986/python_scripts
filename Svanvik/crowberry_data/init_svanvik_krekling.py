import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mytools.met_tools import print_all
from mytools.netcdf_tools import *

# The source data
src = os.environ['DATA']+'/astra_data/observations/krekling_svanvik/Krekling juni aug sept 2019 Photo og Cond.xlsx'

# Open the file
data_krekling = pd.read_excel(src)

k_O3 = 1.67
