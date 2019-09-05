import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Clean up
plt.close('all')
nc_src = os.environ['DATA']+'/astra_data/emep_ctm_results/Base_fullrun.nc'
try:
    data
except NameError:
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src)):
        print(file)
        data = xr.open_dataset(file)
        data_list.append(data)
    # Concat
    data = xr.concat(data_list, dim='time')


# Plot it
levels = np.arange(0,49,4)
fig1 = plt.figure(1)
fig1.canvas.set_window_title("EMEP_POD1")
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)
data['POD1_DF'].where(data['POD1_DF']>0).plot(ax=ax11, levels=levels)
data['POD1_CF'].where(data['POD1_CF']>0).plot(ax=ax12, levels=levels)

#for ax in fig1.axes:
#    ax.set_xlim(1525,2707)
#    ax.set_ylim(-40,-2525)
#ax11.set_xlabel("Time (months)")
#ax11.set_ylabel("$F_{st}^{O_3}$ (mmol m$^{-2}$)")
# Show it
plt.show(block=False)




    
