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
nc_src = os.environ['DATA']+'/astra_data/ctm_results/C3RUN_emep_full_2005/scavenging_monthly/sum_scavenging_2005??_2d.nc'
try:
    data
except NameError:
    data = read_data(nc_src,var='FstO3_avg', datatype='osloctm')

ave2sum = np.array([seconds_in_month(imonth,2005) for imonth in np.arange(1,13)])
ave2sum_stacked = np.stack(np.repeat(ave2sum,80*160)).reshape(12,80,160)

# Plot it
fig1 = plt.figure(1)
fig1.canvas.set_window_title("Local_FstO3_2005")
ax11 = plt.subplot()
(data*ave2sum_stacked*1e-6).sel(lat=49, lon=8 ,method='nearest').groupby('time').sum().plot()
(data*ave2sum_stacked*1e-6).sel(lat=45, lon=7.6 ,method='nearest').groupby('time').sum().plot()
ax11.set_xlabel("Time (months)")
ax11.set_ylabel("$F_{st}^{O_3}$ (mmol m$^{-2}$)")
# Show it
plt.show(block=False)




    
