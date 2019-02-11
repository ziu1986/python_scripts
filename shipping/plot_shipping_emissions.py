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

# Data source
nc_src = os.environ['DATA']+'/CTM3_input_data/EMIS/ECLIPSE_V5a/ship_CLE_emis_*.nc'

try:
    data
except NameError:
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src)):
        print(file)
        data = xr.open_dataset(file)
        data_list.append(data)
# Concatinating the listq
data = xr.concat(data_list, dim='time')

#Plotting
# Clean-up
plt.close('all')
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("total_shipping_emissions")
for i, key in enumerate(data.data_vars):
    plt.subplot(3,3,i+1)
    data.sum(dim='lon')[key].transpose().plot(cmap=plt.cm.afmhot_r)

for ax in fig1.axes[::2]:
    ax.set_xlabel('')
    ax.set_ylabel('')
fig1.axes[6].set_ylabel('Latitude (deg)')
fig1.axes[14].set_xlabel('Time (year)')

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("relative1990_shipping_emissions")
for i, key in enumerate(data.data_vars):
    plt.subplot(3,3,i+1)
    ((data.sum(dim='lon')[key]-data.isel(time=0).sum(dim='lon')[key])/data.isel(time=0).sum(dim='lon')[key]).transpose().plot(cmap=plt.cm.seismic)

for ax in fig2.axes[::2]:
    ax.set_xlabel('')
    ax.set_ylabel('')
fig2.axes[6].set_ylabel('Latitude (deg)')
fig2.axes[14].set_xlabel('Time (year)')

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("zonalave_relative1990_shipping_emissions")
for i, key in enumerate(data.data_vars):
    plt.subplot(3,3,i+1)
    ((data.mean(dim='lon')[key]-data.isel(time=0).mean(dim='lon')[key])/data.isel(time=0).mean(dim='lon')[key]).sel(lat=68,method='nearest').plot(label='%d$^\circ$N'%(data['lat'].sel(lat=68,method='nearest').data))
    ((data.mean(dim='lon')[key]-data.isel(time=0).mean(dim='lon')[key])/data.isel(time=0).mean(dim='lon')[key]).sel(lat=77,method='nearest').plot(label='%d$^\circ$N'%(data['lat'].sel(lat=77,method='nearest').data))
    ((data.mean(dim='lon')[key]-data.isel(time=0).mean(dim='lon')[key])/data.isel(time=0).mean(dim='lon')[key]).sel(lat=52,method='nearest').plot(label='%d$^\circ$N'%(data['lat'].sel(lat=52,method='nearest').data))
    
    ax = plt.gca()
    
    ax.set_title(key,y=0.85,x=0.1)

for ax in fig3.axes:
    ax.set_xlabel('')
    ax.set_ylim(-10,50)
    ax.set_ylabel('')
fig3.axes[3].set_ylabel('(X(t)-X(1990))/X(1990)')
fig3.axes[7].set_xlabel('Time (year)')
fig3.axes[0].legend()

fig6 = plt.figure(6)
fig6.canvas.set_window_title("map_norsca")
plat = (data.lat.sel(lat=68,method='nearest').data, data.lat.sel(lat=77,method='nearest').data)


for i, key in enumerate(data.data_vars):
    minimum = (data.isel(time=10)-data.isel(time=0))[key].min()
    plt.subplot(3,3,i+1,projection=cp.crs.NorthPolarStereo())
    if minimum >= 0:
        (data.isel(time=10)-data.isel(time=0))[key].plot(transform=cp.crs.PlateCarree(), cmap=plt.cm.OrRd, vmin=0,
                                                         cbar_kwargs={'label':' %s (kt)' % (key),
                                                                      'orientation':'vertical',
                                                                      #'fraction':0.046,
                                                                      #'pad':0.04,
                                                                      #'aspect':30
                                                         })
    else:
        (data.isel(time=10)-data.isel(time=0))[key].plot(transform=cp.crs.PlateCarree(),
                                                         cbar_kwargs={'label':' %s (kt)' % (key),
                                                                      'orientation':'vertical',
                                                                      #'fraction':0.046,
                                                                      #'pad':0.04,
                                                                      #'aspect':30
                                                         })
    
    
for ax in fig6.axes[::2]:
    ax.set_extent([0, 22, 56, 80], cp.crs.PlateCarree())
    # Adding lines
    draw_parallels(ax, np.arange(60,91,5))
    draw_meridians(ax, np.arange(-180,181,45))    
    ax.coastlines(resolution='50m')
    ax.plot(np.linspace(0,180,100),np.repeat(plat[0], 100),transform=cp.crs.PlateCarree(),color='blue')
    ax.plot(np.linspace(0,180,100),np.repeat(plat[1], 100),transform=cp.crs.PlateCarree(),color='red')
        
# Show it
plt.show(block=False)
