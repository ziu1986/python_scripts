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
from cartopy.feature import NaturalEarthFeature
import datetime as dt             # Time-objects and maipulation

# Closing plots from previous runs
plt.close('all')

# Path to the file
#src_dir = os.environ['DATA'] + '/CTM3_input_data' + '/metdataOpenIFS/cy38r1nc4_1990/T159N80L60/01/ECopenIFSc38r1_y1990m01d01h00_T159N80L60.nc'
#src_dir = '/uio/hume/student-u41/hanneeby/abel_work/cy38r1nc4_2006/T159N80L60/03/ECopenIFSc38r1_y2006m03d19h12_T159N80L60.nc'
src_dir = '/uio/hume/student-u41/hanneeby/05/ECopenIFSc38r1_y2012m05d06h00_T159N80L60.nc'

# Read the data only once in interactive mode
try:
    data
except NameError:
    for each in sorted(glob.glob(src_dir)):
        print("Reading file %s" % (each))
        data = xr.open_dataset(each)

#data = xr.open_dataset(src_dir)

lats = data['lat'][:]
lons = data['lon'][:]
# Selecting surface pressure
mslp = data['MSLP'][:]
# Correcting the longitudinal coordinates
mslp.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)

# Selecting LS precipitation
lsprec = data['LSPREC'][:]
# Correcting the longitudinal coordinates
lsprec.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)

# Selecting 10 m winds
u10 = data['U10M'][:] 
# Correcting the longitudinal coordinates
u10.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)
v10 = data['V10M'][:] 
# Correcting the longitudinal coordinates
v10.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)

#Selecting 2 m temperature
t2 = data['T2M']
# Correcting the longitudinal coordinates
t2.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)

#Selecting wind level 60 top
u = data['U'] 
# Correcting the longitudinal coordinates
u.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)
v = data['V'] 
# Correcting the longitudinal coordinates
v.coords['lon'] = np.linspace(data.lon[0], data.lon[-1], data.lon.size)

U,V = np.meshgrid(u10, v10) 

# Plotting the data pressure
#fig1 = plt.figure(1, figsize=(16,9))
#ax11 = plt.subplot(projection=cp.crs.PlateCarree())
#(mslp/hecto).plot(ax=ax11, transform=cp.crs.PlateCarree(), cmap=plt.cm.RdYlBu_r)

# Plotting the data wind


fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(projection=cp.crs.PlateCarree())
ax11.barbs(u.lon, u.lat, U, V, transform=cp.crs.PlateCarree(), cmap=plt.cm.RdYlBu_r, length=5, linewidth=1.5)
#(v10).plot(ax=ax11, transform=cp.crs.PlateCarree(), cmap=plt.cm.RdYlBu_r)


ax11.set_global()
ax11.set_aspect('auto')
ax11.set_title('')
ax11.set_xticks(np.arange(-180,181,60), crs=cp.crs.PlateCarree())
ax11.set_yticks(np.arange(-90,91,30), crs=cp.crs.PlateCarree())
ax11.set_xlabel("")
ax11.set_ylabel("")
ax11.coastlines()
plt.show()

"""

# Decorade the plot
for ax in fig1.axes[:-1]:
    ax.set_global()
    ax.set_aspect('auto')
    ax.set_title('')
    ax.set_xticks(np.arange(-180,181,60), crs=cp.crs.PlateCarree())
    ax.set_yticks(np.arange(-90,91,30), crs=cp.crs.PlateCarree())
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.coastlines(zorder=1)

# Set the labels pressure
#ax11.set_xlabel("Latitude (deg)")
#ax11.set_ylabel("Longitude (deg)")
#fig1.axes[-1].set_ylabel("%s (h%s)" % (mslp.name, mslp.units))

# Set the labels wind
ax11.set_xlabel("Latitude (deg)")
ax11.set_ylabel("Longitude (deg)")
fig1.axes[-1].set_ylabel((v10.name, v10.units))




# Show it
#plt.show(block=False)
plt.show()
"""
