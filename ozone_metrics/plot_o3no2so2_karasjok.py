'''
Plot observed and modeled O3, NO2, and SO2 for Karasjok.
Compute climatology.
'''
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
from scipy.optimize import curve_fit
from station_info import station_location

def plot_month_span(ax):
    ax.axvspan(1,31, color='linen')
    ax.axvspan(31+28, 31+28+31, color='linen')
    ax.axvspan(31+28+31+30, 31+28+31+30+31, color='linen')
    ax.axvspan(31+28+31+30+31+30, 31+28+31+30+31+30+31, color='linen')
    ax.axvspan(31+28+31+30+31+30+31+31, 31+28+31+30+31+30+31+31+30, color='linen')
    ax.axvspan(31+28+31+30+31+30+31+31+30+31, 31+28+31+30+31+30+31+31+30+31+30, color='linen')
def plot_month_name(ax, ypos):
    xpos = (1, 31, 31+28, 31+28+31, 31+28+31+30, 31+28+31+30+31,
            31+28+31+30+31+30, 31+28+31+30+31+30+31, 31+28+31+30+31+30+31+31,
            31+28+31+30+31+30+31+31+30, 31+28+31+30+31+30+31+31+30+31, 31+28+31+30+31+30+31+31+30+31+30)
    for i in range(1,13):
        ax.text(xpos[i-1], ypos, get_month_name(i, length=3))
    

def ozone_rescale(**kwarg):
    '''
    Taken from sensitivity studies EMEP_ppgssh and EMEP_ppgssh_ice. For each month.
    '''
    abs_scaling = kwarg.pop('scale', 1)
    scaling_factors = np.array((0.09141522, 0.11777156, 0.12770195, 0.1175085, 0.11039527, 0.06264322,
                                0.0323245, 0.02716563, 0.05162368, 0.06723087, 0.10173934, 0.15230455))
    return(scaling_factors*abs_scaling)

# Close the previous plots
plt.close('all')

# Directories of data
#nc_src_o3 = os.environ['DATA']+'/astra_data/ctm_results/CTM3_oivind/osloctm_ozone*.nc'
#nc_src_so2 = os.environ['DATA']+'/astra_data/ctm_results/CTM3_oivind/osloctm_so2*.nc'
#nc_src_no2 = os.environ['DATA']+'/astra_data/ctm_results/CTM3_oivind/osloctm_no2*.nc'

nc_src_o3 = os.environ['DATA']+'/nird_data/results/OsloCTM3/ozone25/vmr/vmr_ozone*.nc'
nc_src_so2 = os.environ['DATA']+'/nird_data/results/OsloCTM3/ozone25/vmr/vmr_sul*.nc'
nc_src_no2 = os.environ['DATA']+'/nird_data/results/OsloCTM3/ozone25/vmr/vmr_no2*.nc'
src_karasjok_o3 = os.environ['DATA']+'/astra_data/observations/ozone/Karasjok/NO0055R.*ozone*.nas'
src_karasjok_so2 = os.environ['DATA']+'/astra_data/observations/sulfur_dioxide/Karasjok/NO0055R.*pack*.nas'
src_karasjok_no2 = os.environ['DATA']+'/astra_data/observations/nitrogen_dioxide/Karasjok/NO0055R.*nitrogen_dioxide*.nas'


try:
    data_karasjok_o3
except NameError:
    data_karasjok_o3 = []
    data_karasjok_so2 = []
    data_karasjok_no2 = []
    
    for file in sorted(glob.glob(src_karasjok_o3)):
        tmp = read_station_data_ebas(file)
        data_karasjok_o3.append(tmp)
    data_karasjok_o3 = pd.concat(data_karasjok_o3)
    for file in sorted(glob.glob(src_karasjok_so2)):
        tmp = read_station_data_ebas(file, tracer=('SO2',))
        data_karasjok_so2.append(tmp)
    data_karasjok_so2 = pd.concat(data_karasjok_so2)
    for file in sorted(glob.glob(src_karasjok_no2)):
        tmp = read_station_data_ebas(file, tracer=('NO2',))
        data_karasjok_no2.append(tmp)
    data_karasjok_no2 = pd.concat(data_karasjok_no2)

    # Round time index to full hour
    data_karasjok_o3.index = data_karasjok_o3.index.round("h")
    data_karasjok_so2.index = data_karasjok_so2.index.round("h")
    data_karasjok_no2.index = data_karasjok_no2.index.round("h")
try:
    data_o3
except NameError:
    data_o3 = []
    data_so2 = []
    data_no2 = []
    for ifile in sorted(glob.glob(nc_src_o3)):
        print("Reading %s" % (ifile))
        data = xr.open_dataset(ifile).isel(lev=0)
        data_o3.append(data.sel(lat=station_location['Karasjok'].lat,method='nearest').sel(lon=station_location['Karasjok'].lon,method='nearest')*1e9)
    for ifile in sorted(glob.glob(nc_src_so2)):
        print("Reading %s" % (ifile))
        data = xr.open_dataset(ifile).isel(lev=0)
        data_so2.append(data.sel(lat=station_location['Karasjok'].lat,method='nearest').sel(lon=station_location['Karasjok'].lon,method='nearest')*1e9)
    for ifile in sorted(glob.glob(nc_src_no2)):
        print("Reading %s" % (ifile))
        data = xr.open_dataset(ifile).isel(lev=0)
        data_no2.append(data.sel(lat=station_location['Karasjok'].lat,method='nearest').sel(lon=station_location['Karasjok'].lon,method='nearest')*1e9)
    #data_o3 = read_data(nc_src_o3,var='O3',lev=(1,1))
    #data_so2 = read_data(nc_src_so2,var='SO2',lev=(1,1))
    #data_no2 = read_data(nc_src_no2,var='NO2',lev=(1,1))

    # Data selection
    #sel_data_o3 = (data_o3.sel(lat=station_location['Karasjok'].lat,method='nearest').sel(lon=station_location['Karasjok'].lon,method='nearest')*1e9)
    #sel_data_so2 =  (data_so2.sel(lat=station_location['Karasjok'].lat,method='nearest').sel(lon=station_location['Karasjok'].lon,method='nearest')*1e9)
    #sel_data_no2 = (data_no2.sel(lat=station_location['Karasjok'].lat,method='nearest').sel(lon=station_location['Karasjok'].lon,method='nearest')*1e9)
    sel_data_o3 = xr.concat(data_o3, dim='time')
    sel_data_so2 = xr.concat(data_so2, dim='time')
    sel_data_no2 = xr.concat(data_no2, dim='time')
    
# Resampling the scaling factors to daily
#test = pd.Series(1-ozone_rescale(), index=pd.date_range('1988-01-01', periods=12, freq='1M'))
#dates = pd.date_range('1988-01-01', '1988-12-31', freq='D')
#scale = test.reindex(dates, method='bfill')
scale = []
for iyear in np.unique(sel_data_o3.time.dt.year):
    test = pd.Series(1-ozone_rescale(), index=pd.date_range('%s-01-01' % (iyear), periods=12, freq=pd.offsets.MonthBegin(1)))
    dates = pd.date_range('%s-01-01' % (iyear), '%s-12-31 21' % (iyear), freq='3H')
    scale.append(test.reindex(dates, method='ffill'))
scale = pd.concat(scale)

clim_o3 = (sel_data_o3*scale.values).groupby('time.dayofyear').mean()
clim_so2 = sel_data_so2.resample(time='1D').mean().groupby('time.dayofyear').mean()
clim_no2 = sel_data_no2.groupby('time.dayofyear').mean()
climerr_o3 = (sel_data_o3*scale.values).groupby('time.dayofyear').std()/np.sqrt((sel_data_o3*scale.values).groupby('time.dayofyear').sum()/clim_o3)
climerr_so2 = sel_data_so2.resample(time='1D').mean().groupby('time.dayofyear').std()/np.sqrt(sel_data_so2.groupby('time.dayofyear').sum()/clim_so2)
climerr_no2 = sel_data_no2.groupby('time.dayofyear').std()/np.sqrt(sel_data_no2.groupby('time.dayofyear').sum()/clim_no2)
anom_o3 = (sel_data_o3*scale.values).groupby('time.dayofyear')-clim_o3
anom_so2 = sel_data_so2.resample(time='1D').mean().groupby('time.dayofyear')-clim_so2
anom_no2 = sel_data_no2.groupby('time.dayofyear')-clim_no2

clim_karasjok_o3 = data_karasjok_o3.groupby(data_karasjok_o3.index.dayofyear).mean()
clim_karasjok_so2 = data_karasjok_so2.groupby(data_karasjok_so2.index.dayofyear).mean()
clim_karasjok_no2 = data_karasjok_no2.groupby(data_karasjok_no2.index.dayofyear).mean()
climerr_karasjok_o3 = data_karasjok_o3.groupby(data_karasjok_o3.index.dayofyear).std()/np.sqrt(data_karasjok_o3.groupby(data_karasjok_o3.index.dayofyear).sum()/clim_karasjok_o3)
climerr_karasjok_so2 = data_karasjok_so2.groupby(data_karasjok_so2.index.dayofyear).std()/np.sqrt(data_karasjok_so2.groupby(data_karasjok_so2.index.dayofyear).sum()/clim_karasjok_so2)
climerr_karasjok_no2 = data_karasjok_no2.groupby(data_karasjok_no2.index.dayofyear).std()/np.sqrt(data_karasjok_no2.groupby(data_karasjok_no2.index.dayofyear).sum()/clim_karasjok_no2)
#anom_karasjok_o3 = data_karasjok_o3.groupby(data_karasjok_o3.index.dayofyear).apply(lambda x: x - clim_karasjok_o3)
#anom_karasjok_so2 = data_karasjok_so2.apply(lambda x: x - clim_karasjok_so2)
#anom_karasjok_no2 = data_karasjok_no2.apply(lambda x: x - clim_karasjok_no2)

# Plotting
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("ebas_timeseries")
ax11 = plt.subplot(311)
ax12 = plt.subplot(312, sharex=ax11)
ax13 = plt.subplot(313, sharex=ax11)
(sel_data_o3*scale.values)['O3'].plot(ax=ax11, alpha=0.15, color='blue', label='OsloCTM3 v1.0')
sel_data_so2['SO2'].plot(ax=ax12, alpha=0.15, color='blue', label='OsloCTM3 v1.0')
sel_data_no2['NO2'].plot(ax=ax13, alpha=0.15, color='blue', label='OsloCTM3 v1.0')
data_karasjok_o3['O3'].plot(ax=ax11, marker='+', ls='none', color='grey', alpha=0.15, label='Karasjok') #['1998-01-01':]
data_karasjok_so2['SO2'].plot(ax=ax12, marker='+', ls='none', color='grey', alpha=0.15, label='Karasjok')
data_karasjok_no2['NO2'].plot(ax=ax13, marker='+', ls='none', color='grey', alpha=0.15, label='Karasjok')

for ax in fig1.axes[:-1]:
    ax.set_xlabel('')
    ax.set_xticklabels('')
for ax in fig1.axes:   
    ax.set_title('')
    ax.set_ylabel('%s (ppb)' % ax.get_ylabel())
    ax.legend()
# There seems to be an outlier (last data point)?
ax13.set_ylim(ax13.get_ylim()[0], 11)
    
fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("ebas_climatology_profile")
ax21 = plt.subplot(311)
ax22 = plt.subplot(312, sharex=ax21)
ax23 = plt.subplot(313, sharex=ax21)

clim_o3['O3'].plot(ax=ax21, color='blue', label='OsloCTM3 v1.0')
clim_so2['SO2'].plot(ax=ax22, color='blue', label='OsloCTM3 v1.0')
clim_no2['NO2'].plot(ax=ax23, color='blue', label='OsloCTM3 v1.0')
ax21.fill_between(clim_o3.dayofyear.data, clim_o3['O3']-climerr_o3['O3'], clim_o3['O3']+climerr_o3['O3'], alpha=0.25)
ax22.fill_between(clim_so2.dayofyear.data, clim_so2['SO2']-climerr_so2['SO2'], clim_so2['SO2']+climerr_so2['SO2'], alpha=0.25)
ax23.fill_between(clim_no2.dayofyear.data, clim_no2['NO2']-climerr_no2['NO2'], clim_no2['NO2']+climerr_no2['NO2'], alpha=0.25)
#plot_error_bands(ax21,
#                 clim_o3.dayofyear.data,
#                 clim_o3.data_vars,
#                 climerr_o3['O3'].data)
#plot_error_bands(ax22,
#                 clim_so2.dayofyear,
#                 clim_so2.data_vars,
#                 climerr_so2['SO2'].data)
#plot_error_bands(ax23,
#                 clim_no2.dayofyear,
#                 clim_no2.data_vars,
#                 climerr_no2['NO2'].data)

clim_karasjok_o3.plot(ax=ax21,
                      yerr=climerr_karasjok_o3,
                      marker='.', ls='none', color='grey', label='Karasjok')
clim_karasjok_so2.plot(ax=ax22,
                       yerr=climerr_karasjok_so2,
                       marker='.', ls='none', color='grey', label='Karasjok')
clim_karasjok_no2.plot(ax=ax23,
                       yerr=climerr_karasjok_no2,
                       marker='.', ls='none', color='grey', label='Karasjok')

for ax in fig2.axes[:-1]:
    ax.set_xlabel('')
    ax.set_xticklabels('')
for ax in fig2.axes:   
    ax.set_title('')
    ax.set_ylabel('%s (ppb)' % ax.get_ylabel())
    ax.set_xlim(0,367)
    plot_month_span(ax)
    ax.legend()
plot_month_name(ax21, 55)
ax21.set_ylim(0,60)
ax22.set_ylim(0,3)
ax23.set_ylim(0,3)

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("ebas_anomalies")
ax31 = plt.subplot(311)
ax32 = plt.subplot(312, sharex=ax31)
ax33 = plt.subplot(313, sharex=ax31)
anom_o3['O3'].plot(ax=ax31, color='blue', label='OsloCTM3')
anom_so2['SO2'].plot(ax=ax32, color='blue', label='OsloCTM3')
anom_no2['NO2'].plot(ax=ax33, color='blue', label='OsloCTM3')
#anom_karasjok_o3.plot(ax=ax31, marker='.', ls='none', color='grey', label='Karasjok')
#anom_karasjok_so2.plot(ax=ax32, marker='.', ls='none', color='grey', label='Karasjok')
#anom_karasjok_no2.plot(ax=ax33, marker='.', ls='none', color='grey', label='Karasjok')

for ax in fig3.axes[:-1]:
    ax.set_xlabel('')
    ax.set_xticklabels('')
for ax in fig3.axes:   
    ax.set_title('')
    ax.set_ylabel('%s (ppb)' % ax.get_ylabel())
    ax.legend()


fig4 = plt.figure(4, figsize=(9,16))
fig4.canvas.set_window_title("ebas_corr_dens")

ax41 = plt.subplot(311)
ax42 = plt.subplot(312)
ax43 = plt.subplot(313)
hist_o3 = ax41.hist2d(sel_data_o3.sel(time=data_karasjok_o3[2::3].dropna().index)['O3'], data_karasjok_o3[2::3].dropna()['O3'], bins=np.arange(0,81), cmap=plt.cm.hot_r)

hist_so2 = ax42.hist2d(sel_data_so2.resample(time='1D').mean().sel(time=data_karasjok_so2.dropna().index.date)['SO2'], data_karasjok_so2.dropna()['SO2'], bins=(np.arange(0,0.205,0.01), np.arange(0,0.205,0.01)), cmap=plt.cm.hot_r)

hist_no2 = ax43.hist2d(sel_data_no2.resample(time='1D').mean().sel(time=data_karasjok_no2.dropna().index.date)['NO2'], data_karasjok_no2.dropna()['NO2'], bins=(np.arange(0,2.01,0.1), np.arange(0,2.01,0.1)), cmap=plt.cm.hot_r)


ax41.plot(np.arange(0,81),np.arange(0,81), color='grey', ls=':')
cb = fig4.colorbar(hist_o3[3], ax=ax41)
cb.set_label("counts")
ax42.plot(np.arange(0,2),np.arange(0,2), color='grey', ls=':')
cb = fig4.colorbar(hist_so2[3], ax=ax42)
cb.set_label("counts")
cb = fig4.colorbar(hist_no2[3], ax=ax43)
cb.set_label("counts")
ax41.set_ylabel("$[O_3]_{obs}$ (ppb)")
ax41.set_xlabel("$[O_3]_{model} (ppb)$")
ax42.set_ylabel("$[SO_2]_{obs}$ (ppb)")
ax42.set_xlabel("$[SO_2]_{model} (ppb)$")
ax43.set_ylabel("$[NO_2]_{obs}$ (ppb)")
ax43.set_xlabel("$[NO_2]_{model} (ppb)$")

# Show it
plt.show(block=False)
