import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
from matplotlib.dates import date2num
import datetime as dt
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from scipy.stats import linregress     # Get linearregression stats
from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.color_maps import *

b_sigma = True
# Close the previous plots
plt.close('all')

# Directories of data
nc_src = os.environ['DATA']+'/processed_data/CTM3_oivind/osloctm_ozone*.nc'
nc_landsea = os.environ['DATA']+'/landsea_T42.nc'
# Read CTM data
try:
    data
except NameError:
    data = read_data(nc_src,var='O3',lev=(1,1))
mask_landsea = read_data(nc_landsea)
# Select data and compute AOT40
selection_ctm = (data/1e-9).where(
    (data.time.dt.month>=5)&
    (data.time.dt.month<8)&
    ((data/1e-9)>=40),
    drop=True)

selection2_ctm = (data/1e-9).where(
    (data.time.dt.month>=5)&
    (data.time.dt.month<8)&
    (data.time.dt.hour>=8)&
    (data.time.dt.hour<=20)&
    ((data/1e-9)>=40),
    drop=True)
sum04 = ((selection_ctm-40)*6).resample('Y','time',np.nansum)
sum04_2 = ((selection2_ctm-40)*6).resample('Y','time',np.nansum)

data_cyclic, cyclic_lon = ccrs_util.add_cyclic_point(sum04.data, sum04.lon)
data_cyclic = xr.DataArray(data_cyclic, coords=[sum04.time,sum04.lev,sum04.lat,cyclic_lon], dims=['time','lev','lat','lon'])
data_cyclic_2, cyclic_2_lon = ccrs_util.add_cyclic_point(sum04_2.data, sum04_2.lon)
data_cyclic_2 = xr.DataArray(data_cyclic_2, coords=[sum04_2.time,sum04_2.lev,sum04_2.lat,cyclic_2_lon], dims=['time','lev','lat','lon'])
mask_landsea_cyclic, cyclic_lon = ccrs_util.add_cyclic_point(mask_landsea['landsea'].data, mask_landsea.lon)
mask_landsea_cyclic = xr.DataArray(mask_landsea_cyclic, coords=[mask_landsea.lat, cyclic_lon], dims=['lat','lon'])
zonal_mean_land = data_cyclic.where(mask_landsea_cyclic>0).sel(lat=slice(60,80)).mean()*1e-3
zonal_mean_normstd_land = data_cyclic.where(mask_landsea_cyclic>0).sel(lat=slice(60,80)).std()/data_cyclic.where(mask_landsea_cyclic>0).sel(lat=slice(60,80)).mean()
# 
# Plot data
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title('ozone_sum04_allyears')
for i in range(len(np.unique(data.time.dt.year))):
    ax = plt.subplot(5,3,i+1, projection=cp.crs.PlateCarree())
    (sum04*1e-3).isel(time=i).sel(lat=slice(60,80)).plot(ax=ax, transform=cp.crs.PlateCarree(), vmin=0, vmax=14, cmap=plt.cm.Oranges)
for ax, year in zip(fig1.axes[::2], np.unique(data.time.dt.year)):
    ax.set_title("%s" % (year))
    ax.set_extent([-180, 180, 61, 80], cp.crs.PlateCarree())
    ax.set_aspect('auto')
    ax.coastlines()
fig1.axes[11].set_ylabel("SUM04 (ppm)",y=-1)

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("ozone_sum04_mean")
ax21 = plt.subplot(211, projection=cp.crs.PlateCarree())
ax22 = plt.subplot(212, projection=cp.crs.PlateCarree())
(sum04.mean(dim='time')*1e-3).sel(lat=slice(60,80)).plot(ax=ax21, transform=cp.crs.PlateCarree(), vmin=0, vmax=14, cmap=plt.cm.Oranges)
(sum04.std(dim='time')/sum04.mean(dim='time')).sel(lat=slice(60,80)).plot(ax=ax22, transform=cp.crs.PlateCarree(), vmin=0, vmax=4, cmap=plt.cm.RdYlBu_r)
for ax in fig2.axes[:-2]:
    ax.set_title('')
    ax.set_extent([-180, 180, 61, 80], cp.crs.PlateCarree())
    ax.set_aspect('auto')
    ax.coastlines()
fig2.axes[-2].set_ylabel("<SUM04> (ppm)")
fig2.axes[-1].set_ylabel("$\sigma_{SUM04}/<SUM04>$")

fig3 = plt.figure(3, figsize=(16,9))
if b_sigma:
    fig3.canvas.set_window_title("ozone_sum04_sigmalevels")
else:
    fig3.canvas.set_window_title("ozone_sum04_var")
ax31 = plt.subplot(211, projection=cp.crs.PlateCarree())
ax32 = plt.subplot(212, projection=cp.crs.PlateCarree())
(data_cyclic.mean(dim='time')*1e-3).where(mask_landsea_cyclic>0).sel(lat=slice(60,80)).squeeze().plot.contourf(ax=ax31, transform=cp.crs.PlateCarree(),
                                                                                  levels=np.arange(0,14.5),
                                                                                  vmin=0, vmax=14,
                                                                                  cmap=plt.cm.Oranges,
                                                                                  cbar_kwargs={#'orientation':'horizontal',
                                                                                               'label':'<SUM04> (ppm)'})
                                                                                               #'fraction':0.05,
                                                                                               #'pad':0.05,
                                                                                               #'shrink':0.6,
                                                                                               #'aspect':30})

deviation = ((data_cyclic.std(dim='time')/data_cyclic.mean(dim='time')).where(mask_landsea_cyclic>0).sel(lat=slice(60,80)).squeeze()-zonal_mean_normstd_land)
np.fabs(deviation).plot.contourf(ax=ax31, transform=cp.crs.PlateCarree(), levels=np.arange(0,3), 
                                 hatches=[ None, '\\\\', '////'], 
                                 colors='none',
                                 add_colorbar=False)    
                        #cbar_kwargs={'label':'|$\sigma_{SUM04}$/<SUM04>|'})

(data_cyclic_2.mean(dim='time')*1e-3).where(mask_landsea_cyclic>0).sel(lat=slice(60,80)).squeeze().plot.contourf(ax=ax32, transform=cp.crs.PlateCarree(),
                                                                                  levels=np.arange(0,7.5),
                                                                                  vmin=0, vmax=7,
                                                                                  cmap=plt.cm.Oranges,
                                                                                  cbar_kwargs={#'orientation':'horizontal',
                                                                                               'label':'<SUM04> (ppm)'})
                                                                                               #'fraction':0.046,
                                                                                               #'pad':0.05,
                                                                                               #'shrink':0.6,
                                                                                               #'aspect':30})
deviation2 = ((data_cyclic_2.std(dim='time')/data_cyclic_2.mean(dim='time')).where(mask_landsea_cyclic>0).sel(lat=slice(60,80)).squeeze()-zonal_mean_normstd_land)
np.fabs(deviation2).plot.contourf(ax=ax32, transform=cp.crs.PlateCarree(), levels=np.arange(0,3),
                                  hatches=[ None, '\\\\', '////'], 
                                  colors='none',
                                  add_colorbar=False)
                         #cbar_kwargs={'label':'|$\sigma_{SUM04}$/<SUM04>|'})
for ax in fig3.axes[:-2]:
    ax.set_title('')
    ax.set_extent([-180, 180, 62, 80], cp.crs.PlateCarree())
    ax.set_aspect('auto')
    ax.coastlines()
ax31.set_title("0-24h")
ax32.set_title("8-20h")

fig4 = plt.figure(4, figsize=(16,9))
fig4.canvas.set_window_title("ozone_sum04_8-20_mean")
ax41 = plt.subplot(211, projection=cp.crs.PlateCarree())
ax42 = plt.subplot(212, projection=cp.crs.PlateCarree())
(sum04_2.mean(dim='time')*1e-3).sel(lat=slice(60,80)).plot(ax=ax41, transform=cp.crs.PlateCarree(), vmin=0, vmax=7, cmap=plt.cm.Oranges)
(sum04_2.std(dim='time')/sum04_2.mean(dim='time')).sel(lat=slice(60,80)).plot(ax=ax42, transform=cp.crs.PlateCarree(), vmin=0, vmax=4, cmap=plt.cm.RdYlBu_r)
for ax in fig4.axes[:-2]:
    ax.set_title('')
    ax.set_extent([-180, 180, 61, 80], cp.crs.PlateCarree())
    ax.set_aspect('auto')
    ax.coastlines()
fig4.axes[-2].set_ylabel("<SUM04> (ppm)")
fig4.axes[-1].set_ylabel("$\sigma_{SUM04}/<SUM04>$")

# Show it
plt.show(block=False)

#selection_ctm.sel(time='1997-05').groupby('time.day').apply(lambda x : np.sum((x-40)*6)).plot()
#stacked = selection_ctm.sel(time='1997-05').stack(z=('time.day','lat','lon'))
