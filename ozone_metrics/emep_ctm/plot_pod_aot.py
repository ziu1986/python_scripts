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

def plot_data(data, **karg):
    var = karg.pop("var", 'POD1_DF')
    b_diff = karg.pop("diff", False)
    if b_diff:
        levels = karg.pop('level', np.arange(-4,4.1,0.5))
        data2 = data.isel(experiment=0)
        data1 = data.isel(experiment=slice(1,5))
        data = data1-data2
        fp = data[var].plot(transform=crs.crs.PlateCarree(), x='lon', y='lat',
                            levels=levels, extend='max', cmap=plt.cm.coolwarm,
                            col='experiment', col_wrap=4,
                            subplot_kws={'projection':crs.crs.PlateCarree(), 'aspect':6},
                            cbar_kwargs={'aspect':30, 'fraction':0.1, 'pad':0.04},
                            figsize = (16,9))
        fig1 = plt.gcf()
        fig1.canvas.set_window_title("EMEP_Delta%s" % (var))
    else:
        levels = karg.pop('level', np.arange(0,49,4))
        fp = data[var].plot(transform=crs.crs.PlateCarree(), x='lon', y='lat',
                            levels=levels, extend='max', cmap=plt.cm.hot_r,
                            col='experiment', col_wrap=5,
                            subplot_kws={'projection':crs.crs.PlateCarree(), 'aspect':6},
                            cbar_kwargs={'aspect':30, 'fraction':0.1, 'pad':0.04},
                            figsize = (16,9))
        fig1 = plt.gcf()
        fig1.canvas.set_window_title("EMEP_%s" % (var))

    for ax in fig1.axes[:-1]:
        ax.set_extent([0,40,56,71], crs=crs.crs.PlateCarree()) #[18,31,66,70]
        ax.coastlines('10m')
        ax.add_feature(crs.feature.OCEAN, zorder=1, edgecolor='None')
        #ax.set_aspect('auto')
       

# Clean up
plt.close('all')
#nc_src = os.environ['DATA']+'/astra_data/emep_ctm_results/Base_fullrun.nc'
nc_src = os.environ['DATA']+'/nird_data/results/EMEP_MSC-W/'
exp = ('2015-emepctm/Base_fullrun.nc', '2015-emepctm_t15/Base_fullrun.nc', '2015-emepctm_t2/Base_fullrun.nc', '2015-emepctm_t25/Base_fullrun.nc', '2015-emepctm_t4/Base_fullrun.nc')
try:
    data
except NameError:
    data_list = []
    for iexp in exp:
        #data = []
        # Open dataset
        for file in sorted(glob.glob(nc_src+iexp)):
            print(file)
            data = xr.open_dataset(file)
        # Concat
        data_list.append(data)
    # Concatenate the experiements as new dimension
    data_list = xr.concat(data_list,dim='experiment')
    # Assign dimension and discription
    data_list['experiment'] = [0,1.5,2,2.5,4]
    data_list['experiment'].attrs['units'] = 'K'
    data_list['experiment'].attrs['descibtion'] = "Delta Temperature - Compared to present day, representing 2100 global average temperature increase in RCP scenarios."

# Convert EMEP grid
# earth radius is 6370/50=127.4
#p = Proj('+ellps=sphere +a=127.4 +e=0 +proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110')
#, 
#nx = [370, 420]#[56,81]
#ny = [270, 320]#[0,40]
#ni, nj = np.meshgrid(nx, ny)
#nlons, nlats = p(ni, nj, inverse=True)
#upLon,  upLat = np.meshgrid(nlons,nlats)

# Plot it
plot_data(data_list, var='POD1_DF')
plot_data(data_list, var='POD1_CF')
plot_data(data_list, var='POD1_DF', diff=True)
plot_data(data_list, var='POD1_CF', diff=True)

plot_data(data_list, var='EUAOT40_Forests', level=np.arange(2000,10000.1, 1000))
plot_data(data_list, var='EUAOT40_Forests', level=np.arange(-1000,1000.1, 100), diff=True)

plot_data(data_list, var='SURF_ppb_O3', level=np.arange(0,40.1, 5))
plot_data(data_list, var='SURF_ppb_O3', level=np.arange(-1.5,1.6, 0.01), diff=True)

# Show it
#plt.tight_layout()
plt.show(block=False)





    
