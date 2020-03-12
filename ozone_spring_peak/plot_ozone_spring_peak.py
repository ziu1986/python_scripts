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
src = os.environ['DATA']+'/astra_data/observations/ozone/'
src_stations = ('Esrange', 'Pallas', 'Prestebakke')
src_rra = os.environ['DATA']+'/nird_data/reanalysis/Copernicus/ensemble_ozone/SCA_ENSa.2018.O3.yearlyrea.nc'

station_colors = {'Esrange':'blue','Pallas':'black','Prestebakke':'red'}
# Clean up
plt.close('all')

try:
    data
except NameError:
    data = {}
    for station in src_stations:
        if station=='Barrow':
            data.update({station:load_data(src+station+'/*', type="Barrow")})
        else:
            data.update({station:load_data(src+station+'/*.nas')})

    # Load regional model reanalysis 2018 and set time axis
    data_rra = xr.open_dataset(src_rra)
    data_rra['time'] = pd.date_range("2018-01-01", periods=365*24, freq='H')


# Plot it
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)

for istation in src_stations:
    selection = data[istation].where(((data[istation].index.month>=2)&(data[istation].index.month<=5))).dropna()
    ax11.hist(selection, histtype='step', color=station_colors[istation])
    ax12.plot(selection.groupby(selection.index.hour).mean(), color=station_colors[istation])

ax11.axvline(30, ls=':', color='grey', linewidth=5)

fig2 = plt.figure(2, figsize=(16,9))
ax21 = plt.subplot()
for istation in src_stations:
    selection = data[istation].where(((data[istation].index.month>=2)&(data[istation].index.month<=5)&(data[istation].index.hour>10)&(data[istation].index.hour<18)&(data[istation]<30))).dropna()
    ax21.plot(selection.groupby(selection.index.hour).mean(), color=station_colors[istation])
# Show it
plt.show(block=False)
