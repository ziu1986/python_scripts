import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Data source
nc_src = os.environ['DATA']+'/abel/C3RUN_nitrate_corr_2005/avgsav/avgsav_*.nc'

try:
    data
except NameError:
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src)):
        print(file)
        data = xr.open_dataset(file)
        # Define new time coordinates and drop the old one
        #data['time'].reset_coords(drop=True)
        data.coords['time'] = dt.datetime(data['START_TIME'][0],data['START_TIME'][1],data['START_TIME'][2])
        data_list.append(data['temperature'])
# Concatinating the listq
data = xr.concat(data_list, dim='time')

#Plotting
# Clean-up
plt.close('all')

(data.isel(lev=0).where(data.isel(lev=0)-273.15>5)-273.15).squeeze().plot(x='lon', y='lat', col='time', col_wrap=4)
# Display it
plt.show(block=False)


