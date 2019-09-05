import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
from pyproj import Proj
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
#import cartopy.geodesic as cgeo
from matplotlib.transforms import offset_copy
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *

# Data source
nc_src = os.environ['DATA']+'/nird_data/postproc/OsloCTM3/shipping/vmr/*avg*.nc'

try:
    data_list
except NameError:
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src)):
        print(file)
        tmp_data = xr.open_dataset(file)
        data_list.append(tmp_data)
    # Concat
    tmp_data = [ xr.concat(data_list[i*12:i*12+12], dim='time') for i in range(0,4) ]
    # Reding datasets
    data_scandic_slice = [ each.sel(lat=slice(0,40),lon=slice(56,81),drop=True) for each in tmp_data ]
    data_finnmark_slice = [ each.sel(lat=slice(16,33),lon=slice(66,71),drop=True) for each in tmp_data ]

# Clean up
plt.close('all')
label = ("2005_ssp585ship2050", "2005_ssp585ship2100", "2006_ssp585ship2050", "2006_ssp585ship2100")
# Show it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("diff_ozone_future_shipping")
ax11 = plt.subplot()
#ax12 = plt.subplot(212)

for each, i in zip(data_scandic_slice[::2], range(0,4)):
    weights = np.cos(each.lat*np.pi/180)
    print((each.isel(lev=0).mean(dim='lon')*weights*1e9).mean(dim='lat')['O3'])
    (each.isel(lev=0).mean(dim='lon')*weights*1e9).mean(dim='lat')['O3'].plot(ax=ax11, label=label[::2][i])
for each, i in zip(data_finnmark_slice[::2], range(0,4)):
    weights = np.cos(each.lat*np.pi/180)
    print((each.isel(lev=0).mean(dim='lon')*weights*1e9).mean(dim='lat')['O3'])
    (each.isel(lev=0).mean(dim='lon')*weights*1e9).mean(dim='lat')['O3'].plot(ax=ax11, ls='--', label=label[::2][i]+"_NK")

ax11.set_title("")
ax11.set_xlabel("Time (months)")
ax11.set_ylabel("<$\Delta [O_3]$> (ppb)")
#ax12.set_title("Finnmark")
#ax12.set_xlabel("Time (months)")
#ax12.set_ylabel("<$\Delta [O_3]$> (ppb)")
ax11.legend()

for ax in fig1.axes:
    ax.set_xticks(np.arange(0,12))
    ax.set_xticklabels([ get_month_name(each, length=3) for each in np.arange(1,13)])



plt.show(block=False)
