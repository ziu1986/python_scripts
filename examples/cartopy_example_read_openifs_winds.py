####
# Script reads single file from OpenIFS output on Cicero directories
# and plots it.
####

import os, glob, sys
import numpy as np
import pandas as pd               # Analysis of timeseries data
import xarray as xr               # Reading netcdf files and dealing with global fields
import matplotlib.pyplot as plt   # Plotting data
from scipy.constants import *     # Get physics constants
import datetime as dt             # Time-objects and maipulation
import cartopy as cp              # Globe projections
import cartopy.crs as ccrs
import cartopy.util as ccrs_util  # Add cyclic
import cartopy.io.img_tiles as cimgt

# Closing plots from previous runs
plt.close('all')

# Path to the file
src_dir = os.environ['DATA'] + '/CTM3_input_data' + '/metdataOpenIFS/cy38r1nc4_2012/T159N80L60/05/ECopenIFSc38r1_y2012m05d06h00_T159N80L60.nc'

# Read the data only once in interactive mode
try:
    data
except NameError:
    for each in sorted(glob.glob(src_dir)):
        print("Reading file %s" % (each))
        data = xr.open_dataset(each)

# Selecting 10 m winds      
u10 = data['U10M']
v10 = data['V10M']
# Correcting the longitudinal coordinates
u10.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)
v10.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)
# Select an area
u10_fenno = u10.sel(lat=slice(67.6,71.4),lon=slice(19.4,31.4))
v10_fenno = v10.sel(lat=slice(67.6,71.4),lon=slice(19.4,31.4))

# Plotting the data pressure
fig1 = plt.figure(1, figsize=(16,9))
stamen_terrain = cimgt.Stamen('terrain-background')
ax11 = fig1.add_subplot(1,1,1, projection=stamen_terrain.crs) #ccrs.PlateCarree()
ax11.set_extent([19.4,31.4,67.6,71.4], crs=ccrs.PlateCarree())
ax11.add_image(stamen_terrain, 8)
#ax11.barbs(u10_fenno.lon.data, u10_fenno.lat.data, u10_fenno.data, v10_fenno.data, transform=cp.crs.PlateCarree())
ax11.quiver(u10_fenno.lon.data, u10_fenno.lat.data, u10_fenno.data, v10_fenno.data, transform=cp.crs.PlateCarree())

# Plotting the data wind

plt.show(block=False)
