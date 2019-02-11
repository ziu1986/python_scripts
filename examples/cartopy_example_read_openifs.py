####
# Script reads single file from OpenIFS output on Cicero directories
# and plots it.
####

import os, glob, sys
import numpy as np
import pandas as pd               # Analysis of timeseries data
import xarray as xr               # Reading netcdf files and dealing with global fields
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp              # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
import datetime as dt             # Time-objects and maipulation

# Closing plots from previous runs
plt.close('all')

# Path to the file
src_dir = os.environ['DATA'] + '/CTM3_input_data' + '/metdataOpenIFS/cy38r1nc4_1990/T159N80L60/01/ECopenIFSc38r1_y1990m01d01h00_T159N80L60.nc'

# Read the data only once in interactive mode
try:
    data
except NameError:
    for each in sorted(glob.glob(src_dir)):
        print("Reading file %s" % (each))
        data = xr.open_dataset(each)


# Selecting surface pressure
mslp = data['MSLP']
# Correcting the longitudinal coordinates
mslp.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)

# Plotting the data
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(projection=cp.crs.PlateCarree())

(mslp/hecto).plot(ax=ax11, transform=cp.crs.PlateCarree(), cmap=plt.cm.RdYlBu_r)

# Decorade the plot
for ax in fig1.axes[:-1]:
    ax.set_global()
    ax.set_aspect('auto')
    ax.set_title('')
    ax.set_xticks(np.arange(-180,181,60), crs=cp.crs.PlateCarree())
    ax.set_yticks(np.arange(-90,91,30), crs=cp.crs.PlateCarree())
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.coastlines()
# Set the labels
ax11.set_xlabel("Latitude (deg)")
ax11.set_ylabel("Longitude (deg)")
fig1.axes[-1].set_ylabel("%s (h%s)" % (mslp.name, mslp.units))
# Show it
plt.show(block=False)

