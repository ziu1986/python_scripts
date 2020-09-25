import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from matplotlib.dates import date2num # Convert dates to matplotlib axis coords
from matplotlib import dates
from scipy import fftpack
from scipy import stats
#from read_sunspots import *
from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.ozone_tools import *
#from mytools.cartopy_tools import scale_bar
from mytools.station_info import station_location

#---------------------------------------------------------------------------------------------------------------------------------
# Read data
#src = os.environ['DATA']+'/astra_data/observations/ozone/'
#src_stations = ('Esrange', 'Pallas', 'Prestebakke')
#src_rra = os.environ['DATA']+'/nird_data/reanalysis/Copernicus/ensemble_ozone/SCA_ENSa.2018.O3.yearlyrea.nc'
src_ref = os.environ['DATA']+'/nird_data/models/results/OsloCTM3/drydepdevel/version2/C3RUN_mOSaic/trop_tracer/*.nc'
src_ice = os.environ['DATA']+'/nird_data/models/results/OsloCTM3/drydepdevel/version2/C3RUN_mOSaic_ice/trop_tracer/*.nc'
src_gs = "GROWING_SEASON_2005.nc"

src_obs = os.environ['DATA']+'/astra_data/observations/ozone/Pallas/*2005*'
station = station_location['Pallas']

# Clean up
plt.close('all')

try:
    data
except NameError:
    data_list = []
    data_list_ref = []

    for each in sorted(glob.glob(src_ref))[100:161]:
        print(each)
        data_list_ref.append(xr.open_dataset(each)['O3'].sel(lat=station.lat, lon=station.lon, lev=1000, method='nearest'))

    for each in sorted(glob.glob(src_ice))[100:161]:
        print(each)
        data_list.append(xr.open_dataset(each)['O3'].sel(lat=station.lat, lon=station.lon, lev=1000, method='nearest'))

    data_ref = xr.concat(data_list_ref, dim='time')
    data = xr.concat(data_list, dim='time')
    # Reading aproximate SGS
    gs = xr.open_dataset(src_gs).sel(lat=station.lat, lon=station.lon, method='nearest')
    # Reading ozone obs
    obs = read_station_data_ebas(glob.glob(src_obs)[0])
    
# Scaling factor g/m^3 -> vmr (ppb)
s_factor = 0.5*1e-3*1e9
# Plot it
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)

#(data_ref*s_factor).groupby("time.dayofyear").mean().plot(ax=ax11, label='ref')
#(data*s_factor).groupby("time.dayofyear").mean().plot(ax=ax11, label='ice')
ra_ref = (data_ref*s_factor).rolling(time=3*24, center=True)
ra_ice = (data*s_factor).rolling(time=3*24, center=True)
ra_ref.mean().plot(ax=ax11, label='ref', color='blue')
ra_ice.mean().plot(ax=ax11, label='ice', color='blue')
obs.rolling(24, center=True).mean().plot(ax=ax11)

#((data-data_ref)*s_factor).groupby("time.dayofyear").mean().plot(ax=ax12)
ax12.plot(ra_ref.mean().dropna('time').time.values, np.gradient(ra_ref.mean().dropna('time').values), label='ref', color='blue')
ax12.plot(ra_ice.mean().dropna('time').time.values, np.gradient(ra_ice.mean().dropna('time').values), label='ice', color='blue')

SGS = (gs.where(gs['GDAY']==1, drop=True))['GDAY'].time
for ax in fig1.axes:
    ax.axvline(pd.to_datetime(SGS.values), color='red')
    #ax.axvline(SGS.time.dt.dayofyear.values, color='red')
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_title('')
ax12.set_ylabel("$O_3$ (ppb)", y=1.1)
ax12.set_xlabel("Time (day of year)")
ax12.set_ylim(-0.5,0.5)
ax11.set_xticklabels("")
ax11.set_ylim(0,80)
ax11.legend()

# Show it
plt.show(block=False)
