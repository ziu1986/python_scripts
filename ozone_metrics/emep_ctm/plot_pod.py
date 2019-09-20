import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
from pyproj import Proj
import cartopy as crs        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Clean up
plt.close('all')
#nc_src = os.environ['DATA']+'/astra_data/emep_ctm_results/Base_fullrun.nc'
nc_src = os.environ['DATA']+'/abel/2015-emepctm_16cpu/Base_fullrun.nc'
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

# Convert EMEP grid
# earth radius is 6370/50=127.4
p = Proj('+ellps=sphere +a=127.4 +e=0 +proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110')
#, 
nx = [370, 420]#[56,81]
ny = [270, 320]#[0,40]
ni, nj = np.meshgrid(nx, ny)
#nlons, nlats = p(ni, nj, inverse=True)
#upLon,  upLat = np.meshgrid(nlons,nlats)

# Plot it
levels = np.arange(0,49,4)
fig1 = plt.figure(1)
fig1.canvas.set_window_title("EMEP_POD1")
ax11 = plt.subplot(211, projection=crs.crs.PlateCarree())
ax12 = plt.subplot(212, projection=crs.crs.PlateCarree())
data['POD1_DF'].where(data['POD1_DF']>0).plot(ax=ax11, x='lon', y='lat', transform=crs.crs.PlateCarree(), levels=levels, extend='max', cmap=plt.cm.hot_r)
data['POD1_CF'].where(data['POD1_CF']>0).plot(ax=ax12, x='lon', y='lat', transform=crs.crs.PlateCarree(), levels=levels[1:], extend='max', cmap=plt.cm.hot_r)

for ax in fig1.axes[:-2]:
    ax.set_extent([0,40,56,81], crs=crs.crs.PlateCarree()) #[18,31,66,70]
    ax.coastlines('10m')
    ax.add_feature(crs.feature.OCEAN, zorder=1, edgecolor='None')
    
#ax11.set_xlabel("Time (months)")
#ax11.set_ylabel("$F_{st}^{O_3}$ (mmol m$^{-2}$)")

#fig2 = plt.figure(2)
#fig2.canvas.set_window_title("EMEP_POD1_xy")
#ax21 = plt.subplot(211)
#ax22 = plt.subplot(212)
#data['POD1_DF'].where(data['POD1_DF']>0).plot(ax=ax21, levels=levels[1:], extend='max', cmap=plt.cm.hot_r)
#data['POD1_CF'].where(data['POD1_CF']>0).plot(ax=ax22, levels=levels[2:], extend='max', cmap=plt.cm.hot_r)

#for ax in fig2.axes[:2]:
#    ax.set_xlim(1500,2500)
#    ax.set_ylim(-2888,-500)

# Show it
plt.show(block=False)




    
