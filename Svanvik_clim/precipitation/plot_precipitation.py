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
src_precip_svanvik = os.environ['DATA']+'/astra_data/observations/svanvik_accu_precip_2009-2020.xlsx'
src_precip_cru_clim =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/SCA_precip_flux_WFDE5_CRU_climatology.nc'
src_precip_cru_clim_std =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/SCA_precip_flux_WFDE5_CRU_climatology_std.nc'
src_precip_cru =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/2017-2018/SCA_precip_flux_WFDE5_CRU_*.nc'
# Load data
try:
    data_svanvik
except NameError:
    
    data_svanvik = pd.read_excel(src_precip_svanvik).dropna()
    data_svanvik[u"sum(precipitation_amount P1D)"] = [float(each.replace(',','.')) for each in data_svanvik[u"sum(precipitation_amount P1D)"].values]
    data_svanvik.index = [np.datetime64("%s-%s-%s" % (itime[-4:], itime[-7:-5], itime[:2])) for itime in data_svanvik['time'].dropna().values ]

    data_cru_clim = xr.open_dataset(src_precip_cru_clim)
    data_cru_clim_std = xr.open_dataset(src_precip_cru_clim_std)

svanvik_precip = data_svanvik[u"sum(precipitation_amount P1D)"]
svanvik_precip_cru = data_cru_clim.sel(lat=station_location['Svanvik'].lat, lon=station_location['Svanvik'].lon, method='nearest')['Rainf']
svanvik_precip_cru_std = data_cru_clim_std.sel(lat=station_location['Svanvik'].lat, lon=station_location['Svanvik'].lon, method='nearest')['Rainf']

svanvik_ap_cru = (svanvik_precip_cru*60**2).resample(time='1D').sum()
svanvik_ap_cru_std = np.sqrt(((svanvik_precip_cru_std*60**2)**2).resample(time='1D').sum())

svanvik_ap_cru_test_list = []
svanvik_ap_cru_test_std_list = []
for iyear in range(2014,2020):
    if iyear != 2016:
        #print(iyear)
        svanvik_ap_cru_test = svanvik_ap_cru.drop(pd.date_range('2016-02-29',freq='D', periods=1), dim='time').copy()
        svanvik_ap_cru_test['time'] = pd.date_range('%d-01-01' % iyear, periods=365, freq='D')
        svanvik_ap_cru_test_std = svanvik_ap_cru_std.drop(pd.date_range('2016-02-29',freq='D', periods=1), dim='time').copy()
        svanvik_ap_cru_test_std['time'] = pd.date_range('%d-01-01' % iyear, periods=365, freq='D')
    else:
        #print(iyear, "ja")
        svanvik_ap_cru_test = svanvik_ap_cru.copy()
        svanvik_ap_cru_test['time'] = pd.date_range('%d-01-01' % iyear, periods=366, freq='D')
        svanvik_ap_cru_test_std = svanvik_ap_cru_std.copy()
        svanvik_ap_cru_test_std['time'] = pd.date_range('%d-01-01' % iyear, periods=366, freq='D')
        
    svanvik_ap_cru_test_list.append(svanvik_ap_cru_test)
    svanvik_ap_cru_test_std_list.append(svanvik_ap_cru_test_std)

svanvik_ap_cru_test_list = xr.concat(svanvik_ap_cru_test_list, dim='time')
svanvik_ap_cru_test_std_list = xr.concat(svanvik_ap_cru_test_std_list, dim='time')
svanvik_pull = ((svanvik_precip['2014':'2019']-svanvik_ap_cru_test_list.sel(time=svanvik_precip['2014':'2019'].index.values))/svanvik_ap_cru_test_std_list.sel(time=svanvik_precip['2014':'2019'].index.values))

# Plot it
fig1 = plt.figure(1)
fig1.canvas.set_window_title("plot_precipitation_anomalies_svanvik")
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)


svanvik_ap_cru_test_list.plot(ax=ax11, color='black', label="CRU reanalysis", ls='--')
svanvik_precip['2014':].plot(ax=ax11, color='blueviolet', label='Svanvik', ls='-')
(svanvik_ap_cru_test_list+svanvik_ap_cru_test_std_list).plot(ax=ax11, color='black', alpha=0.45)
(svanvik_ap_cru_test_list-svanvik_ap_cru_test_std_list).plot(ax=ax11, color='black', alpha=0.45)

ax11.legend(ncol=2, loc='upper left')
ax11.set_xlabel("")
ax11.set_ylabel("$Precip$ (mm day)")

pull_plot = svanvik_pull.where((svanvik_pull.index.month>6)&(svanvik_pull.index.month<10)).dropna()
pull_plot.plot(ax=ax12, ls='None', marker='x')

ax12.axhline(2, color='red', ls=':')
ax12.axhline(1, color='orange', ls=':')
ax12.axhline(-2, color='red', ls=':')
ax12.axhline(-1, color='orange', ls=':')

for iyear in range(2014,2020):
    print(iyear, pull_plot.where(pull_plot['%d' % iyear]<-2).dropna().count(), pull_plot.where(pull_plot['%d' % iyear]<-1).dropna().count())

ax12.set_xlabel("Time (years)")
ax12.set_ylabel("$\Delta_{clim} Precip / \sigma Precip_{clim}$")
ax12.set_ylim(-2.5, 2.5)

#fig2 = plt.figure("precipitation_hist")
#ax21 = plt.subplot()

#svanvik_precip.hist(ax=ax21, bins=np.arange(0,50,0.1), density=True, histtype='step', color='blueviolet')
#ax21.hist((svanvik_precip_cru*60**2).resample(time='1D').sum(), density=True, bins=np.arange(0,50,0.1), histtype='step')


# Show it
plt.show(block=False)
