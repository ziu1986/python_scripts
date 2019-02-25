import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from scipy import fftpack
from read_sunspots import *
from mytools.met_tools import *
from mytools.netcdf_tools import *

def load_data(src):
    data = []
    for file in sorted(glob.glob(src)):
        print("Reading file %s" % (file))
        tmp = read_station_data_ebas(file)
        data.append(tmp['O3'])   
    # Concatenate the lists
    print('Concatenating data...')
    data = pd.concat(data)
    # Round to full hours
    data.index = data.index.round("h")
    return(data)
    
# Clean up
plt.close('all')

# Load functions to read NOAA data
#execfile("../BrXplo/read_station_data.py")
#from read_station_data import read_station_data_noaa
src = os.environ['DATA']+'/astra_data/observations/ozone/'
src_stations = ('Barrow', 'Esrange', 'Janiskoski', 'Jergul', 'Karasjok', 'Pallas', 'Prestebakke', 'Svanvik')

try:
    data_barrow
except NameError:
    data_barrow = []
    for file in sorted(glob.glob(src+src_stations[0]+'/*')):
        if int(file[-4:]) < 2003:
            print(file)
            tmp = (read_station_data_noaa(file, utc=-9, start_data=28))
        elif int(file[-4:]) < 2012:
            print(file)
            tmp = (read_station_data_noaa(file, utc=-9, station='other', column=5))
        data_barrow.append(tmp)
    # Concatenate the lists
    data_barrow = pd.concat(data_barrow)
    # Round time index to full hour
    data_barrow.index = data_barrow.index.round("h")

try:
    data
except NameError:
    data = []
    for station in src_stations[1:]:
        data.append(load_data(src+station+'/*.nas'))
        
data_prestebakke = data[5]
data_jergkara = pd.concat(data[2:4])
    
# Spectral analysis
from scipy import fftpack
fft_barrow = fftpack.fft(data_barrow.resample('1M').mean().fillna(method='ffill'))
freqs_barrow = fftpack.fftfreq(len(fft_barrow))

fft_prestebakke = fftpack.fft(data_prestebakke.resample('1M').mean().fillna(method='ffill'))
freqs_prestebakke = fftpack.fftfreq(len(fft_prestebakke))

fft_jergkara = fftpack.fft(data_jergkara.resample('1M').mean().fillna(method='ffill'))
freqs_jergkara = fftpack.fftfreq(len(fft_jergkara))

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_timeseries_ltobs")
ax11 = plt.subplot(311)
ax12 = plt.subplot(312, sharex=ax11)
ax13 = plt.subplot(313, sharex=ax11)
data_barrow.plot(ax=ax11, ls='None', marker='x', label='Utqiagvik (USA)')
data_prestebakke.plot(ax=ax12, ls='None', marker='.', label='Prestebakke (NOR)', color='red')
data_jergkara.plot(ax=ax13, ls='None', marker='+', label='Jergul/Karasjok (NOR)', color='orange')

ax13.set_xlabel("Time (year)")
ax12.set_ylabel("[$O_3$] (ppb)")

for ax in fig1.axes:
    ax.set_ylim(0,100)
    ax.legend()

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("ozone_timeseries_fenoscandicobs")
ax21 = plt.subplot(411)
ax22 = plt.subplot(412, sharex=ax21)
ax23 = plt.subplot(413, sharex=ax21)
ax24 = plt.subplot(414, sharex=ax21)
ax21.plot(data[0].index, data[0], ls='None', marker='+', label='Esrange (SWE)', color='blue')
ax22.plot(data[4].index, data[4], ls='None', marker='+', label='Pallas (FIN)', color='black')
data_jergkara.plot(ax=ax23, ls='None', marker='+', label='Jergul/Karasjok (NOR)', color='orange')
ax24.plot(data[-1].index, data[-1], ls='None', marker='+', label='Svanvik (NOR)', color='blueviolet')

ax24.set_xlabel("Time (year)")
ax22.set_ylabel("[$O_3$] (ppb)", y=1)
for ax in fig2.axes:
    ax.set_ylim(0,100)
    ax.legend()

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("ozone_climatology_fenoscandicobs")
ax31 = plt.subplot()

data[0].groupby(data[0].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Esrange (SWE)', color='blue')
data[4].groupby(data[4].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Pallas (FIN)', color='black')
data_jergkara.groupby(data_jergkara.index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Jergul/Karasjok (NOR)', color='orange')
data_prestebakke.groupby(data_prestebakke.index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Prestebakke (NOR)', color='red')
data[-1].groupby(data[-1].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Svanvik (NOR)', color='blueviolet')
data[1].groupby(data[1].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Janiskoski (RUS)', color='grey')

ax31.set_xlabel("Time (day of year)")
ax31.set_ylabel("[$O_3$] (ppb)")
for ax in fig3.axes:
    ax.set_ylim(0,100)
    ax.legend()

if(False):
    fig2.canvas.set_window_title("fequency_spectrum")
    ax21 = plt.subplot()
    ax21.stem(1/freqs_barrow/12, np.abs(fft_barrow)/np.abs(fft_barrow).max(), label='Utqiagvik (USA)')
    markerline, stemlines, baseline = ax21.stem(1/freqs_prestebakke/12, np.abs(fft_prestebakke)/np.abs(fft_prestebakke).max(), label='Prestebakke (NOR)')
    plt.setp(stemlines, color='red')
    plt.setp(markerline, color='red')
    markerline, stemlines, baseline = ax21.stem(1/freqs_jergkara/12, np.abs(fft_jergkara)/np.abs(fft_jergkara).max(), label='Jergul/Karasjok (NOR)')
    plt.setp(stemlines, color='orange')
    plt.setp(markerline, color='orange')
    ax21.set_xlim(0,40)
    ax21.set_ylabel("Normalized Amplitude")
    ax21.set_xlabel("Frequency (years)")

    ax21.legend()
# Show it
plt.show(block=False)



