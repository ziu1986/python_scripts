import os, sys
import pandas as pd
import numpay as np
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.station_info import *

# Clean up
plt.close('all')

# source
src_temp_svanvik = os.environ['DATA']+'/astra_data/observations/svanvik_temp_deg_2013-2019.xls'

src_temp_cru =  os.environ['DATA']+''
