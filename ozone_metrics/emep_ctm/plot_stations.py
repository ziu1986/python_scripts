import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
from pyproj import Proj
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Clean up
plt.close('all')
nc_src = os.environ['DATA']+'/abel/emepctm_nested/sites_2015.nc'
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
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot()

for istation in data['O3'].isel(station=slice(30,35)):
    istation.plot(ax=ax11)

# Show it
plt.show(block=False)

