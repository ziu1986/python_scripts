import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import xarray as xr
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from matplotlib.dates import date2num # Convert dates to matplotlib axis coords
from matplotlib import dates

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#import cartopy.feature as cfeature
#import cartopy.io.img_tiles as cimgt

from mytools.ozone_tools import *
from mytools.plot_tools import print_all


plt.close('all')

dataset = 'CAMSAQ'

src = {'CAMS': os.environ['DATA']+'/nird_data/reanalysis/ECMWF/CAMS_reanalysis/netcdf/VMR/vmr_cams_r_o3_ml60_climatology.nc',
       'MACC': os.environ['DATA']+'/nird_data/reanalysis/ECMWF/MACC_reanalysis/netcdf/VMR/vmr_macc_r_o3_ml60_3h_climatology.nc',
       'CAMSAQ': os.environ['DATA']+'/nird_data/reanalysis/Copernicus/ensemble_ozone/SCA_ENSa.yearlyrea.climatology.nc'}

var = {'CAMS': 'go3',
       'MACC': 'go3',
       'CAMSAQ': 'o3'}

data = xr.open_dataset(src[dataset])
#selection = data.sel(latitude=slice(71.5,65), longitude=slice(14,33)).isel(time=1460)*1e9
selection = data.sel(latitude=slice(71.5,65), longitude=slice(14,33))
selection_seasonal_means = selection.groupby('time.season').mean('time')*1e9
selection_annual_mean = selection.mean(dim='time')*1e9


# This is the map projection we want to plot *onto*
map_proj = ccrs.PlateCarree(central_longitude=18.5)

levels = np.arange(0,61,2.5)

lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()

fig2 = selection_seasonal_means[var[dataset]].plot(levels=levels, transform=ccrs.PlateCarree(), robust=True, col='season', col_wrap=2,figsize=(16,8), cbar_kwargs={'label':'$[O_3]$ (ppb)', 'shrink':.7, 'anchor':(0.5,0.5), 'drawedges':True}, subplot_kws={'projection': map_proj})

for ax in fig2.axes.flat:
    ax.coastlines()
    ax.set_extent([13.75, 33.5, 64.75, 71.75])
    ax.set_xticks(np.arange(14,34,2), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(65, 72,  2), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_xlabel("")
    ax.set_ylabel("")


fig2.axes.flat[-2].set_xlabel("Longitude ($^\circ E$)", x=1)
fig2.axes.flat[-2].set_ylabel("Latitude ($^\circ N$)", y=1)
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

fig2.fig.canvas.set_window_title("ozone_seasonal_average_%s" % dataset)
#plt.tight_layout()
plt.show(block=False)
