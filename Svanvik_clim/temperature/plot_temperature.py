import os, sys
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.station_info import *


# Clean up
plt.close('all')

# Source
src_temp_svanvik = os.environ['DATA']+'/astra_data/observations/svanvik_temp_deg_2013-2019.xls'
src_temp_cru_clim =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/SCA_Tair_WFDE5_CRU_climatology.nc'
src_temp_cru_clim_std =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/SCA_Tair_WFDE5_CRU_climatology_std.nc'
src_temp_cru =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/2017-2018/SCA_Tair_WFDE5_CRU_*.nc'
# Load data
try:
    data_svanvik
except NameError:
    
    data_svanvik = pd.read_excel(src_temp_svanvik)
    data_svanvik.index = data_svanvik[u'Fra-tid']

    data_cru_clim = xr.open_dataset(src_temp_cru_clim)
    data_cru_clim_std = xr.open_dataset(src_temp_cru_clim_std)

svanvik_temp = (273.15+data_svanvik[u"Svanvik | Ambient Temperature | degC"].where((data_svanvik[u"Svanvik | Ambient Temperature | degC"]>-50) & (data_svanvik[u"Svanvik | Ambient Temperature | degC"]<50)).dropna())
svanvik_temp_cru = data_cru_clim.sel(lat=station_location['Svanvik'].lat, lon=station_location['Svanvik'].lon, method='nearest')['Tair']
svanvik_temp_cru_std = data_cru_clim_std.sel(lat=station_location['Svanvik'].lat, lon=station_location['Svanvik'].lon, method='nearest')['Tair']

svanvik_temp_cru_test_list = []
svanvik_temp_cru_test_std_list = []
for iyear in range(2014,2020):
    if iyear != 2016:
        #print(iyear)
        svanvik_temp_cru_test = svanvik_temp_cru.drop(pd.date_range('2016-02-29',freq='H', periods=24), dim='time').copy()
        svanvik_temp_cru_test['time'] = pd.date_range('%d-01-01' % iyear, periods=365*24, freq='H')
        svanvik_temp_cru_test_std = svanvik_temp_cru_std.drop(pd.date_range('2016-02-29',freq='H', periods=24), dim='time').copy()
        svanvik_temp_cru_test_std['time'] = pd.date_range('%d-01-01' % iyear, periods=365*24, freq='H')
    else:
        #print(iyear, "ja")
        svanvik_temp_cru_test = svanvik_temp_cru.copy()
        svanvik_temp_cru_test['time'] = pd.date_range('%d-01-01' % iyear, periods=366*24, freq='H')
        svanvik_temp_cru_test_std = svanvik_temp_cru_std.copy()
        svanvik_temp_cru_test_std['time'] = pd.date_range('%d-01-01' % iyear, periods=366*24, freq='H')
        
    svanvik_temp_cru_test_list.append(svanvik_temp_cru_test)
    svanvik_temp_cru_test_std_list.append(svanvik_temp_cru_test_std)

svanvik_temp_cru_test_list = xr.concat(svanvik_temp_cru_test_list, dim='time')
svanvik_temp_cru_test_std_list = xr.concat(svanvik_temp_cru_test_std_list, dim='time')
svanvik_pull = ((svanvik_temp['2014':'2019']-svanvik_temp_cru_test_list.sel(time=svanvik_temp['2014':'2019'].index.values))/svanvik_temp_cru_test_std_list.sel(time=svanvik_temp['2014':'2019'].index.values))

# Plot it
fig1 = plt.figure(1)
fig1.canvas.set_window_title("plot_temperature_anomalies_svanvik")
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)

svanvik_temp.plot(ax=ax11, color='blueviolet', label='Svanvik')
svanvik_temp_cru_test_list.plot(ax=ax11, color='black', label="CRU reanalysis")
(svanvik_temp_cru_test_list+svanvik_temp_cru_test_std_list).plot(ax=ax11, color='black', alpha=0.45)
(svanvik_temp_cru_test_list-svanvik_temp_cru_test_std_list).plot(ax=ax11, color='black', alpha=0.45)

ax11.legend(ncol=2, loc='upper left')

svanvik_pull.plot(ax=ax12)

ax12.axhline(4, color='red', ls=':')
ax12.axhline(2, color='orange', ls=':')
ax12.axhline(-4, color='red', ls=':')
ax12.axhline(-2, color='orange', ls=':')

for iyear in range(2014,2020):
    print(iyear, svanvik_pull.where(svanvik_pull['%d' % iyear]>4).dropna().count(), svanvik_pull.where(svanvik_pull['%d' % iyear]>2).dropna().count())

ax11.set_xlabel("Time (years)")
ax11.set_ylabel("$T_{2m}$ (K)")
ax12.set_xlabel("Time (months)")
ax12.set_ylabel("$\Delta_{clim} T_{2m} / \sigma T_{2m, clim}$")

# Show it
plt.show(block=False)
