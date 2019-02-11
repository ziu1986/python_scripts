import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.color_maps import *

def is_mjj(month):
    return(month >=5) & (month<=6)

# Clean up
plt.close('all')
# Data source
nc_src = os.environ['DATA']+'/processed_data/CTM3_oivind/osloctm_ozone*.nc'
nc_landsea = os.environ['DATA']+'/landalt_T42.ctm'
# Reading the data
try:
    data
except NameError:
    data = read_data(nc_src,var='O3',lev=(1,1))

lines_landsea = open(nc_landsea).readlines()
mask_landsea_raw = []
for each in lines_landsea[1:821]:
    # Split the lines
    col = each.split()
    # Seperate ozone data
    for each in col:
        mask_landsea_raw.append(float(each))
mask_landsea = xr.DataArray(np.array(mask_landsea_raw).reshape(64,128), coords=[data.lat.data, data.lon.data], dims=['lat','lon'])
mask_landsea.name = 'landsea'
# Selecting data
monthly_mean = (data*1e9).sel(lat=slice(60,80)).resample('D', 'time', how=np.max).resample('M','time', how=np.mean)
try:
    sum_mm_aug
except NameError:
    sum_mm_may = []
    sum_mm_june = []
    sum_mm_july = []
    sum_mm_august = []
    for iyear in np.unique(data.time.dt.year.data):
        sum_mm_may.append(monthly_mean.sel(time='%s-05' % (iyear)).data)
        sum_mm_june.append(monthly_mean.sel(time='%s-06' % (iyear)).data)
        sum_mm_july.append(monthly_mean.sel(time='%s-07' % (iyear)).data)
        sum_mm_august.append(monthly_mean.sel(time='%s-08' % (iyear)).data)

    clim_may = xr.DataArray(np.mean(sum_mm_may,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    clim_june = xr.DataArray(np.mean(sum_mm_june,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    clim_july = xr.DataArray(np.mean(sum_mm_july,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    clim_august = xr.DataArray(np.mean(sum_mm_august,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    
    clim_std_may = xr.DataArray(np.std(sum_mm_may,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    clim_std_june = xr.DataArray(np.std(sum_mm_june,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    clim_std_july = xr.DataArray(np.std(sum_mm_july,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    clim_std_august = xr.DataArray(np.std(sum_mm_august,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])

    max_may = xr.DataArray(np.max(sum_mm_may,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    max_june = xr.DataArray(np.max(sum_mm_june,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    max_july = xr.DataArray(np.max(sum_mm_july,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    max_august = xr.DataArray(np.max(sum_mm_august,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    min_may = xr.DataArray(np.min(sum_mm_may,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    min_june = xr.DataArray(np.min(sum_mm_june,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    min_july = xr.DataArray(np.min(sum_mm_july,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])
    min_august = xr.DataArray(np.min(sum_mm_august,axis=0).squeeze(), coords=[monthly_mean.lat.data, monthly_mean.lon.data], dims=['lat', 'lon'])

# Plotting
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_monthly_mean_noon")
ax11 = plt.subplot(411, projection=cp.crs.PlateCarree())
ax12 = plt.subplot(412, projection=cp.crs.PlateCarree())
ax13 = plt.subplot(413, projection=cp.crs.PlateCarree())
ax14 = plt.subplot(414, projection=cp.crs.PlateCarree())
clim_may.plot(ax=ax11, transform=cp.crs.PlateCarree(), vmin=0, vmax=60, cmap=plt.cm.RdYlBu_r)
clim_june.plot(ax=ax12, transform=cp.crs.PlateCarree(), vmin=0, vmax=60, cmap=plt.cm.RdYlBu_r)
clim_july.plot(ax=ax13, transform=cp.crs.PlateCarree(), vmin=0, vmax=60, cmap=plt.cm.RdYlBu_r)
clim_august.plot(ax=ax14, transform=cp.crs.PlateCarree(), vmin=0, vmax=60, cmap=plt.cm.RdYlBu_r)

for ax in fig1.axes[:-4]:
    ax.set_extent([-180, 180, 61, 80], cp.crs.PlateCarree())
    ax.set_aspect('auto')
    # Adding lines
    draw_parallels(ax, np.arange(60,91,5))
    draw_meridians(ax, np.arange(-180,181,45))    
    ax.coastlines()
ax11.set_title('May')
ax12.set_title('June')
ax13.set_title('July')
ax14.set_title('August')
fig1.axes[-2].set_ylabel('<$O_{3}$> (ppb)', y=1)


fig2 = plt.figure(2, figsize=(16,9))
#mask_landsea.where(mask_landsea > 0 ).plot()
fig2.canvas.set_window_title("ozone_monthly_mean_noon_std")
ax21 = plt.subplot(411, projection=cp.crs.PlateCarree())
ax22 = plt.subplot(412, projection=cp.crs.PlateCarree())
ax23 = plt.subplot(413, projection=cp.crs.PlateCarree())
ax24 = plt.subplot(414, projection=cp.crs.PlateCarree())
clim_std_may.plot(ax=ax21, transform=cp.crs.PlateCarree(), vmin=0, vmax=10, cmap=plt.cm.RdYlBu_r)
clim_std_june.plot(ax=ax22, transform=cp.crs.PlateCarree(), vmin=0, vmax=10, cmap=plt.cm.RdYlBu_r)
clim_std_july.plot(ax=ax23, transform=cp.crs.PlateCarree(), vmin=0, vmax=10, cmap=plt.cm.RdYlBu_r)
clim_std_august.plot(ax=ax24, transform=cp.crs.PlateCarree(), vmin=0, vmax=10, cmap=plt.cm.RdYlBu_r)
for ax in fig2.axes[:-4]:
    ax.set_extent([-180, 180, 61, 80], cp.crs.PlateCarree())
    ax.set_aspect('auto')
    # Adding lines
    draw_parallels(ax, np.arange(60,91,5))
    draw_meridians(ax, np.arange(-180,181,45))    
    ax.coastlines()
ax21.set_title('May')
ax22.set_title('June')
ax23.set_title('July')
ax24.set_title('August')
fig2.axes[-2].set_ylabel('$\sigma_{O_{3}}$ (ppb)', y=1)


fig3 = plt.figure(3, figsize=(16,9))
#mask_landsea.where(mask_landsea > 0 ).plot()
fig3.canvas.set_window_title("ozone_monthly_mean_noon_stdDIVmean")
ax31 = plt.subplot(411, projection=cp.crs.PlateCarree())
ax32 = plt.subplot(412, projection=cp.crs.PlateCarree())
ax33 = plt.subplot(413, projection=cp.crs.PlateCarree())
ax34 = plt.subplot(414, projection=cp.crs.PlateCarree())
(clim_std_may/clim_may).plot(ax=ax31, transform=cp.crs.PlateCarree(), vmin=0, vmax=0.4, cmap=plt.cm.RdYlBu_r)
(clim_std_june/clim_june).plot(ax=ax32, transform=cp.crs.PlateCarree(), vmin=0, vmax=0.4, cmap=plt.cm.RdYlBu_r)
(clim_std_july/clim_july).plot(ax=ax33, transform=cp.crs.PlateCarree(), vmin=0, vmax=0.4, cmap=plt.cm.RdYlBu_r)
(clim_std_august/clim_august).plot(ax=ax34, transform=cp.crs.PlateCarree(), vmin=0, vmax=0.4, cmap=plt.cm.RdYlBu_r)
for ax in fig3.axes[:-4]:
    ax.set_extent([-180, 180, 61, 80], cp.crs.PlateCarree())
    ax.set_aspect('auto')
    # Adding lines
    draw_parallels(ax, np.arange(60,91,5))
    draw_meridians(ax, np.arange(-180,181,45))    
    ax.coastlines()
ax31.set_title('May')
ax32.set_title('June')
ax33.set_title('July')
ax34.set_title('August')
fig3.axes[-2].set_ylabel('$\sigma_{O_{3}}$/<$O_{3}$>', y=1)


fig4 = plt.figure(4, figsize=(16,9))
fig4.canvas.set_window_title("ozone_monthly_mean_noon_max_2sigma")
ax41 = plt.subplot(411, projection=cp.crs.PlateCarree())
ax42 = plt.subplot(412, projection=cp.crs.PlateCarree())
ax43 = plt.subplot(413, projection=cp.crs.PlateCarree())
ax44 = plt.subplot(414, projection=cp.crs.PlateCarree())
max_may.plot(ax=ax41, transform=cp.crs.PlateCarree(), vmin=0, vmax=60, cmap=plt.cm.RdYlBu_r)
max_june.plot(ax=ax42, transform=cp.crs.PlateCarree(), vmin=0, vmax=60, cmap=plt.cm.RdYlBu_r)
max_july.plot(ax=ax43, transform=cp.crs.PlateCarree(), vmin=0, vmax=60, cmap=plt.cm.RdYlBu_r)
max_august.plot(ax=ax44, transform=cp.crs.PlateCarree(), vmin=0, vmax=60, cmap=plt.cm.RdYlBu_r)

((max_may-clim_may)/clim_std_may).plot.contourf(ax=ax41, transform=cp.crs.PlateCarree(),
                                                levels=np.arange(0,5), hatches=[ '....', '\\\\', None, None], 
                                                colors='none',
                                                add_colorbar=False)
((max_june-clim_june)/clim_std_june).plot.contourf(ax=ax42, transform=cp.crs.PlateCarree(),
                                                levels=np.arange(0,5), hatches=[ '....', '\\\\', None, None], 
                                                colors='none',
                                                add_colorbar=False)
((max_july-clim_july)/clim_std_july).plot.contourf(ax=ax43, transform=cp.crs.PlateCarree(),
                                                levels=np.arange(0,5), hatches=[ '....', '\\\\', None, None], 
                                                colors='none',
                                                add_colorbar=False)
((max_august-clim_august)/clim_std_august).plot.contourf(ax=ax44, transform=cp.crs.PlateCarree(),
                                                levels=np.arange(0,5), hatches=[ '....', '\\\\', None, None], 
                                                colors='none',
                                                add_colorbar=False)

for ax in fig4.axes[:-4]:
    ax.set_extent([-180, 180, 61, 80], cp.crs.PlateCarree())
    ax.set_aspect('auto')
    # Adding lines
    draw_parallels(ax, np.arange(60,91,5))
    draw_meridians(ax, np.arange(-180,181,45))    
    ax.coastlines()
ax41.set_title('May')
ax42.set_title('June')
ax43.set_title('July')
ax44.set_title('August')
fig4.axes[-2].set_ylabel('$O^{max}_{3}$ (ppb)', y=1)

fig5 = plt.figure(5, figsize=(16,9))
fig5.canvas.set_window_title("ozone_monthly_mean_noon_zonal")
ax51 = plt.subplot(411)
ax52 = plt.subplot(412)
ax53 = plt.subplot(413)
ax54 = plt.subplot(414)
#mask_landsea.where(mask_landsea>0).plot(ax=ax51)
min_may.mean(dim='lon').plot(ax=ax51, ls='--')
min_june.mean(dim='lon').plot(ax=ax52, ls='--')
min_july.mean(dim='lon').plot(ax=ax53, ls='--')
min_august.mean(dim='lon').plot(ax=ax54, ls='--')
max_may.mean(dim='lon').plot(ax=ax51, ls='--')
max_june.mean(dim='lon').plot(ax=ax52, ls='--')
max_july.mean(dim='lon').plot(ax=ax53, ls='--')
max_august.mean(dim='lon').plot(ax=ax54, ls='--')

clim_may.where(mask_landsea>0).mean(dim='lon').plot(ax=ax51, label='land')
clim_june.where(mask_landsea>0).mean(dim='lon').plot(ax=ax52, label='land')
clim_july.where(mask_landsea>0).mean(dim='lon').plot(ax=ax53, label='land')
clim_august.where(mask_landsea>0).mean(dim='lon').plot(ax=ax54, label='land')
plot_error_bands(ax51, clim_may.lat.data, clim_may.where(mask_landsea>0).mean(dim='lon').data, clim_std_may.where(mask_landsea>0).mean(dim='lon').data, color='black')
plot_error_bands(ax52, clim_june.lat.data, clim_june.where(mask_landsea>0).mean(dim='lon').data, clim_std_june.where(mask_landsea>0).mean(dim='lon').data, color='black')
plot_error_bands(ax53, clim_july.lat.data, clim_july.where(mask_landsea>0).mean(dim='lon').data, clim_std_july.where(mask_landsea>0).mean(dim='lon').data, color='black')
plot_error_bands(ax54, clim_august.lat.data, clim_august.where(mask_landsea>0).mean(dim='lon').data, clim_std_august.where(mask_landsea>0).mean(dim='lon').data, color='black')

clim_may.where(mask_landsea==0).mean(dim='lon').plot(ax=ax51, label='sea')
clim_june.where(mask_landsea==0).mean(dim='lon').plot(ax=ax52, label='sea')
clim_july.where(mask_landsea==0).mean(dim='lon').plot(ax=ax53, label='sea')
clim_august.where(mask_landsea==0).mean(dim='lon').plot(ax=ax54, label='sea')
plot_error_bands(ax51, clim_may.lat.data, clim_may.where(mask_landsea==0).mean(dim='lon').data, clim_std_may.where(mask_landsea==0).mean(dim='lon').data, color='c')
plot_error_bands(ax52, clim_june.lat.data, clim_june.where(mask_landsea==0).mean(dim='lon').data, clim_std_june.where(mask_landsea==0).mean(dim='lon').data, color='c')
plot_error_bands(ax53, clim_july.lat.data, clim_july.where(mask_landsea==0).mean(dim='lon').data, clim_std_july.where(mask_landsea==0).mean(dim='lon').data, color='c')
plot_error_bands(ax54, clim_august.lat.data, clim_august.where(mask_landsea==0).mean(dim='lon').data, clim_std_august.where(mask_landsea==0).mean(dim='lon').data, color='c')

for ax in fig5.axes:
    ax.set_ylim(0,60)
    ax.text(80,30,'max', color='red')
    ax.text(80,20,'min', color='blue')
ax51.set_title('May', x=0.95, y=0.85)
ax52.set_title('June', x=0.95, y=0.85)
ax53.set_title('July', x=0.95, y=0.85)
ax54.set_title('August', x=0.95, y=0.85)
ax52.set_ylabel("<O$_{3}$> (ppb)", y=-0.25)
ax54.set_xlabel("Latitude (deg)")
ax51.legend(loc='lower left')
# Show it
plt.show(block=False)
