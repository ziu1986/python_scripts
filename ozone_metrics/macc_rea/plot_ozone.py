import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.station_info import station_location

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
    data
except NameError:
    nc_src = os.environ['DATA']+"/astra_data/ECMWF/MACC_reanalysis/netcdf/VMR/*climatology.nc" 
    data = read_date(nc_src)


# Load climatologies from observations
import pickle
with open('../obs_climatologies.pkl', 'rb') as input:
    climatology_nord = pickle.load(input)
    climatology_svanvik = pickle.load(input)
    climatology_prestebakke = pickle.load(input)

sample_macc = {}
sample_x = {}
sample_spl_nord =  climatology_nord(data.time.dt.dayofyear)
sample_spl_prestebakke =  climatology_prestebakke(data.time.dt.dayofyear)

for istation in ('Esrange', 'Pallas', 'Jergul', 'Svanvik', 'Prestebakke'):
    sample_macc[istation] = (data['go3']*1e9).sel(latitude=station_location[istation].lat, longitude=station_location[istation].lon, method='nearest').groupby(data.time.dt.dayofyear).mean()
    #sample_x[istation] = 


# Plot it
plt.close('all')
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_climatology_macc")
sample_macc['Prestebakke'].plot(color='red', label='Prestebakke')

for istation, color, marker in zip(('Esrange', 'Pallas', 'Jergul', 'Svanvik'), ('orange','black','blue', 'blueviolet'), ('x','^','v','+')):
    sample_macc[istation].plot(color=color, label=istation, marker=marker, ls='None')

# Get current axis object
ax = plt.gca()
# Climatology from obs
ax.plot(np.linspace(1,367), climatology_nord(np.linspace(1,367)), ls='--', label='climatology: nord')
ax.plot(np.linspace(1,367), climatology_prestebakke(np.linspace(1,367)), ls='--', label='climatology: Prestebakke')

ax.set_xlabel("Time (day of year)")
ax.set_ylabel("[$O_3$] (ppb)")
for ax in fig1.axes:
    ax.set_ylim(0,60)
    ax.legend(loc='lower left')

plot_month_span(ax)
plot_month_name(ax, ypos=58)

fig2 = plt.figure(2, figsize=(10,12))
fig2.canvas.set_window_title("ozone_climatology_fenoscandic_obs_residuals")

ax21 = plt.subplot(311)
ax22 = plt.subplot(312, sharex=ax21)
ax23 = plt.subplot(313, sharex=ax21)

ax21.set_title("Esrange", y=0.85, x=0.05)
ax22.set_title("Pallas", y=0.85, x=0.05)
ax23.set_title("Prestebakke", y=0.85, x=0.05)


hist21_2018 = ax21.hist((sample_macc['Esrange']-sample_spl_nord).values, density=True, histtype='step', color='blue', label='Esrange 2018')
hist22_2018 = ax22.hist((sample_macc['Pallas']-sample_spl_nord).values, density=True, histtype='step', color='black', label='Pallas 2018')
hist23_2018 = ax23.hist((sample_macc['Prestebakke']-sample_spl_prestebakke).values, density=True, histtype='step', color='red', label='Prestebakke 2018')

#ax21.plot(x_sample_esrange, pdf_esrange, color='black', label='Skew normal fit')
#stats_text(ax21, stat_esrange, fit_esrange, ypos=0.4)

#ax22.plot(x_sample_pallas, pdf_pallas, color='black', label='Skew normal fit')
#stats_text(ax22, stat_pallas, fit_pallas, ypos=0.4)

#ax23.plot(x_sample_prestebakke, pdf_prestebakke, color='black', label='Skew normal fit')
#stats_text(ax23, stat_prestebakke, fit_prestebakke, ypos=0.4)

ax23.set_xlabel("$<O_3>_{daily}-O_3^{clim}$ (ppb)")

# Show it
plt.show(block=False)
