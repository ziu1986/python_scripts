import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import *     # Get physics constants
import matplotlib.pyplot as plt
import datetime as dt
#from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.station_info import station_location
from mytools.plot_tools import *

def read_date(input_files):
    data_list = []
    print("Reading")
    for file_name in sorted(glob.glob(input_files)):
        print(file_name)
        data = xr.open_dataset(file_name)
        data_list.append(data)
        

    # Concat list
    data = xr.concat(data_list, dim='time')
    return(data)

def select(ag_data, lat,lon):
    ag_data_sel = {}
    for each in ag_data:
        ag_data_sel[each] = ag_data[each].sel(latitude=lat,longitude=lon, method='nearest')
    return(ag_data_sel)

try:
    data_macc
except NameError:
    nc_src = os.environ['DATA']+"/astra_data/ECMWF/MACC_reanalysis/netcdf/VMR/*climatology.nc" 
    data_macc = read_date(nc_src)
    nc_src = os.environ['DATA']+"/astra_data/ECMWF/CAMS_reanalysis/netcdf/VMR/*climatology.nc" 
    data_cams = read_date(nc_src)


# Load climatologies from observations
import pickle
with open( os.environ['PY_SCRIPTS']+'/ozone_metrics/ozone_anomalies/obs_climatologies.pkl', 'rb') as input:
    climatology_nord = pickle.load(input)
    climatology_svanvik = pickle.load(input)
    climatology_prestebakke = pickle.load(input)
    yerr_mean_nord = pickle.load(input)
    yerr_mean_svanvik = pickle.load(input)
    yerr_mean_prestebakke = pickle.load(input)

doy_macc = {}
doy_cams = {}   
sample_spl_nord =  climatology_nord(np.unique(data_macc.time.dt.dayofyear))
sample_spl_prestebakke =  climatology_prestebakke(np.unique(data_macc.time.dt.dayofyear))

for istation in ('Esrange', 'Pallas', 'Jergul', 'Svanvik', 'Prestebakke'):
    doy_macc[istation] = (data_macc['go3']*1e9).sel(latitude=station_location[istation].lat, longitude=station_location[istation].lon, method='nearest').groupby(data_macc.time.dt.dayofyear).mean()
    doy_cams[istation] = (data_cams['go3']*1e9).sel(latitude=station_location[istation].lat, longitude=station_location[istation].lon, method='nearest').groupby(data_cams.time.dt.dayofyear).mean()
    


# Plot it
plt.close('all')
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_climatologies")
ax11 = plt.subplot(121)
ax12 = plt.subplot(122)


for istation, color, marker in zip(('Esrange', 'Pallas', 'Jergul', 'Svanvik'), ('orange','black','blue', 'blueviolet'), ('o','^','v','d')):
    doy_macc[istation].plot(ax=ax11, color=color, label=istation, marker=marker, ls='None')
    doy_cams[istation].plot(ax=ax12, color=color, label=istation, marker=marker, fillstyle='none', ls='None')

#
doy_macc['Prestebakke'].plot(ax=ax11, color='red', label='Prestebakke', alpha=0.5)
doy_cams['Prestebakke'].plot(ax=ax12, color='red', label='Prestebakke', alpha=0.5)

for ax in fig1.axes:
    
    # Climatology from obs
    ax.plot(np.linspace(1,367), climatology_nord(np.linspace(1,367)), ls='--', label='obs_clim: Nordkalotten')
    ax.plot(np.linspace(1,367), climatology_prestebakke(np.linspace(1,367)), ls='-.', label='obs_clim: Prestebakke')
    plot_error_bands(ax, np.arange(1,367), climatology_nord(np.arange(1,367)), yerr_mean_nord.values, color='black',ls='None')
    plot_error_bands(ax, np.arange(1,367), climatology_prestebakke(np.arange(1,367)), yerr_mean_prestebakke.values, color='grey',ls='None')
    ax.set_ylim(0,60)
    ax.legend(loc='lower left', ncol=2)
    ax.set_xlabel("Time (day of year)")
    ax.set_ylabel("[$O_3$] (ppb)")
    plot_month_span(ax)
    plot_month_name(ax, ypos=58)

ax11.set_title("MACC reanalysis")
ax12.set_title("CAMS reanalysis")

fig2 = plt.figure(2, figsize=(10,12))
fig2.canvas.set_window_title("ozone_climatology_fenoscandic_obs_residuals")

ax21 = plt.subplot(321)
ax23 = plt.subplot(323, sharex=ax21)
ax25 = plt.subplot(325, sharex=ax21)

ax22 = plt.subplot(322)
ax24 = plt.subplot(324, sharex=ax22)
ax26 = plt.subplot(326, sharex=ax22)

ax21.set_title("Esrange", y=0.85, x=0.15)
ax23.set_title("Pallas", y=0.85, x=0.15)
ax25.set_title("Prestebakke", y=0.85, x=0.15)

ax22.set_title("Esrange", y=0.85, x=0.15)
ax24.set_title("Pallas", y=0.85, x=0.15)
ax26.set_title("Prestebakke", y=0.85, x=0.15)


hist21_2018 = ax21.hist((doy_macc['Esrange']-sample_spl_nord).values, density=True, histtype='step', color='blue', label='Esrange 2018')
hist23_2018 = ax23.hist((doy_macc['Pallas']-sample_spl_nord).values, density=True, histtype='step', color='black', label='Pallas 2018')
hist25_2018 = ax25.hist((doy_macc['Prestebakke']-sample_spl_prestebakke).values, density=True, histtype='step', color='red', label='Prestebakke 2018')

hist22_2018 = ax22.hist((doy_cams['Esrange']-sample_spl_nord).values, density=True, histtype='step', color='blue', label='Esrange 2018')
hist24_2018 = ax24.hist((doy_cams['Pallas']-sample_spl_nord).values, density=True, histtype='step', color='black', label='Pallas 2018')
hist26_2018 = ax26.hist((doy_cams['Prestebakke']-sample_spl_prestebakke).values, density=True, histtype='step', color='red', label='Prestebakke 2018')
#ax21.plot(x_sample_esrange, pdf_esrange, color='black', label='Skew normal fit')
#stats_text(ax21, stat_esrange, fit_esrange, ypos=0.4)

#ax22.plot(x_sample_pallas, pdf_pallas, color='black', label='Skew normal fit')
#stats_text(ax22, stat_pallas, fit_pallas, ypos=0.4)

#ax23.plot(x_sample_prestebakke, pdf_prestebakke, color='black', label='Skew normal fit')
#stats_text(ax23, stat_prestebakke, fit_prestebakke, ypos=0.4)

ax26.set_xlabel("$<O_3>_{daily}-O_3^{clim}$ (ppb)", x=-0.25)

# Show it
plt.show(block=False)
