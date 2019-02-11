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


b_plot = True
# Clean up
plt.close('all')
#nc_src = os.environ['DATA']+"/CTM3_Ozone_data/CIXPAG_????/ctm3_icbc_glo/gas*.nc"
nc_src = os.environ['DATA']+'/processed_data/CTM3_oivind/osloctm_ozone*.nc'
try:
    data
except NameError:
    data = read_data(nc_src,var='O3',lev=(1,1))

# Mean ozone in northern Scandinavia
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_norsca_timeseries")
ave_ozone = data.sel(lev=1).sel(lat=68, method='nearest').sel(lon=slice(19,28)).mean(dim=('lon'))/1e-9
ave_ozone.plot()
ave_ozone.where(ave_ozone>=40).plot(ls='', marker='+')
ax11 = plt.gca()
ax11.set_xlabel("Time (years)")
ax11.set_ylabel("$O_3$ (ppb)")

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("ozone_norsca_resampling")
ax21 = plt.subplot(121)
ax22 = plt.subplot(122)
for iyear in np.arange(1997,1998):
    ts_ozone = ave_ozone.sel(time=slice('%s-05-01'% (iyear),'%s-08-11'% (iyear)))
    ts_ozone_rs = ts_ozone.resample(time='1H').interpolate()
    ts_ozone_rs_aot40_full = ts_ozone_rs.where((ts_ozone_rs>=40))
    ts_ozone_rs_aot40 = ts_ozone_rs.where((ts_ozone_rs>=40)&(ts_ozone_rs['time.hour'] >= 8)&(ts_ozone_rs['time.hour']<=20))
    ts_ozone.plot(ax=ax21, label='6-hourly')
    ts_ozone_rs.where(ts_ozone_rs>=40).plot(ax=ax21, ls='', marker='+', label='upsampled 1-hourly above 40 ppb [0-24]')
    ts_ozone_rs_aot40.plot(ax=ax21, ls='', marker='x', label='upsampled 1-hourly above 40 ppb [8-20]')
    (ts_ozone_rs_aot40_full-40).groupby('time.dayofyear').reduce(np.nansum).plot(ax=ax22, color='red', label='[0-24]h')
    (ts_ozone_rs_aot40-40).groupby('time.dayofyear').reduce(np.nansum).plot(ax=ax22, color='black', label='[8-20]h')
ax21.legend()
ax22.legend(loc='center right')
ax21.set_xlabel("")
ax21.set_ylabel("$O_3$ (ppb)")
ax22.set_xlabel("Day of Year")
ax22.set_ylabel("$\Sigma O_3$ (ppb)")
ax22.text(120,205, get_month_name((dt.date(1997,1,1)+dt.timedelta(120)).month), size='large')
ax22.axvspan(120,151,alpha=0.25,color='orange')
ax22.text(151,205, get_month_name((dt.date(1997,1,1)+dt.timedelta(151)).month), size='large')
ax22.text(181,205, get_month_name((dt.date(1997,1,1)+dt.timedelta(181)).month), size='large')
ax22.axvspan(181,212,alpha=0.25,color='orange')
ax22.text(212,205, get_month_name((dt.date(1997,1,1)+dt.timedelta(212)).month), size='large')

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("ozone_norsca_aot40")
ax31 = plt.subplot(211)
ax32 = plt.subplot(212)
aot40_list_full = []
aot40_list = []
for iyear in np.arange(1997,2011):
    ts_ozone = ave_ozone.sel(time=slice('%s-05-01'% (iyear),'%s-08-11'% (iyear)))
    ts_ozone_rs = ts_ozone.resample(time='1H').interpolate()
    ts_ozone_rs_cumsum_full = ts_ozone_rs.where((ts_ozone_rs>=40))
    ts_ozone_rs_aot40_full = (ts_ozone_rs_cumsum_full-40).groupby('time.dayofyear').reduce(np.nansum)
    aot40_list_full.append(ts_ozone_rs_aot40_full)
    ts_ozone_rs_cumsum = ts_ozone_rs.where((ts_ozone_rs>=40)&(ts_ozone_rs['time.hour']>=8)&(ts_ozone_rs['time.hour']<=20))
    ts_ozone_rs_aot40 = (ts_ozone_rs_cumsum-40).groupby('time.dayofyear').reduce(np.nansum)
    aot40_list.append(ts_ozone_rs_aot40)
    #ts_ozone_rs_aot40.plot(ax=ax31,ls='', marker='x')
aot40_test_full = xr.concat(aot40_list_full,dim='time')
mask = ~np.isnan(aot40_test_full.data)
filtered_data = [d[m] for d, m in zip(aot40_test_full.data.T, mask.T)]
ax31.violinplot(filtered_data,
                positions=aot40_test_full.dayofyear.data,
                points=len(aot40_test_full.time),
                showmeans=True,
                showmedians=True)
aot40_test_full.reduce(np.nanmean, dim='time').plot(ax=ax31,ls='', marker='x', label='mean')
aot40_test_full.reduce(np.nanmedian, dim='time').plot(ax=ax31,ls='', marker='+', label='median')

aot40_test = xr.concat(aot40_list,dim='time')
mask = ~np.isnan(aot40_test.data)
filtered_data = [d[m] for d, m in zip(aot40_test.data.T, mask.T)]
ax32.violinplot(filtered_data,
                positions=aot40_test.dayofyear.data,
                points=len(aot40_test.time),
                showmeans=True,
                showmedians=True)
aot40_test.mean(dim='time').plot(ax=ax32,ls='', marker='x', label='mean')
aot40_test.median(dim='time').plot(ax=ax32,ls='', marker='+', label='median')
for ax in fig3.axes:
    ax.legend()
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_ylim(-10,270)
    ax.text(120,255, get_month_name((dt.date(1997,1,1)+dt.timedelta(120)).month), size='large')
    ax.axvspan(120,151,alpha=0.25,color='orange')
    ax.text(151,255, get_month_name((dt.date(1997,1,1)+dt.timedelta(151)).month), size='large')
    ax.text(181,255, get_month_name((dt.date(1997,1,1)+dt.timedelta(181)).month), size='large')
    ax.axvspan(181,212,alpha=0.25,color='orange')
    ax.text(212,255, get_month_name((dt.date(1997,1,1)+dt.timedelta(212)).month), size='large')
ax31.set_title("[0-24]h")
ax32.set_title("[8-20]h")
ax32.set_xlabel("Day of Year")
ax32.set_ylabel("$\Sigma O_3$ (ppb)", y=1)

fig4 = plt.figure(4,figsize=(16,9))
fig4.canvas.set_window_title("ozone_norsca_sumaot40")
aot40_test_full.reduce(np.nansum, dim='dayofyear').plot(color='red', label='[0-24]h')
aot40_test.reduce(np.nansum, dim='dayofyear').plot(color='black', label='[8-20]h')
ax41 = plt.gca()
ax41.set_title("May 1 - August 11")
ax41.set_ylabel("SumAOT40 $O_3$ (ppb)")
ax41.set_xlabel("Time (years)")
ax41.set_xticklabels(np.arange(1995,2012,2))
ax41.legend()


# Plot it
if b_plot:
    plat, plon = data.sel(lev=1).sel(lat=68, method='nearest').sel(lon=slice(19,28)).lat.data, data.sel(lev=1).sel(lat=68, method='nearest').sel(lon=slice(19,28)).lon.data
    delta_lat = np.abs(data.lat.data[0]-data.lat.data[1])
    delta_lon = np.abs(data.lon.data[0]-data.lon.data[1])
    plat = np.repeat(plat, len(plon))
    
    fig5 = plt.figure(5, figsize=(16,9))
    fig5.canvas.set_window_title("ozone_map")
    for itime in range(4):
        ax = plt.subplot(2,2,itime+1,projection=cp.crs.PlateCarree())
        (data.sel(lev=1).isel(time=itime+720)/1e-9).plot(ax=ax,transform=cp.crs.PlateCarree(), y='lat', x='lon', vmin=0, vmax=100,
                                                     cbar_kwargs={'label':'$O_3$ (ppb)',
                                                                  'orientation':'vertical',
                                                                  'fraction':0.046,
                                                                  'pad':0.04,
                                                                  'aspect':30})
        
    for ax in fig5.axes[::2]:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title('')
        fig5.axes[4].set_ylabel('Latitude (deg)', y=1)
        fig5.axes[4].set_xlabel('Longitude (deg)', x=1.25)
        ax.coastlines()
        ax.set_global()
        ax.set_xticks(np.arange(-180,181,60), crs=cp.crs.PlateCarree())
        ax.set_yticks(np.arange(-90,91,45), crs=cp.crs.PlateCarree())
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)

    fig6 = plt.figure(6)
    fig6.canvas.set_window_title("ozone_map_norsca")
    #ax21 = plt.subplot(projection=cp.crs.NorthPolarStereo())
    cyclic_data, cyclic_lon = ccrs_util.add_cyclic_point(data.data, data.lon)
    data_cyclic = xr.DataArray(cyclic_data, coords=[data.time,data.lev,data.lat,cyclic_lon], dims=['time','lev','lat','lon'])
    
    for itime in range(4):
        ax = plt.subplot(2,2,itime+1,projection=cp.crs.NorthPolarStereo())
        (data_cyclic.isel(time=itime+720).sel(lev=1)/1e-9).plot(ax=ax,transform=cp.crs.PlateCarree(), vmin=0, vmax=70,
                                                            cbar_kwargs={'label':'$O_3$ (ppb)',
                                                                         'orientation':'vertical',
                                                                         #'fraction':0.046,
                                                                         #'pad':0.04,
                                                                         #'aspect':30
                                                            })
        ax.plot((plon[0]-0.5*delta_lon,
                 plon[0]-0.5*delta_lon,
                 plon[-1]+0.5*delta_lon,
                 plon[-1]+0.5*delta_lon,
                 plon[0]-0.5*delta_lon),
                (plat[0]-0.5*delta_lat,
                 plat[-1]+0.5*delta_lat,
                 plat[-1]+0.5*delta_lat,
                 plat[0]-0.5*delta_lat,
                 plat[0]-0.5*delta_lat),transform=cp.crs.PlateCarree(),color='orange')
        ax.plot((24.6,25.217),(69.45,69.467),transform=cp.crs.PlateCarree(),color='orange',ls='',marker='.')
        #ax.plot((0,22,52.75,0),(56,56,74.25,80),transform=cp.crs.PlateCarree(),color='orange',ls='',marker='.')
    for ax in fig6.axes[::2]:
        ax.set_extent([0, 22, 56, 80], cp.crs.PlateCarree())
        # Adding lines
        draw_parallels(ax, np.arange(60,91,5))
        draw_meridians(ax, np.arange(-180,181,45))    
        ax.coastlines(resolution='50m')
        
    
# Show it
plt.show(block=False)
