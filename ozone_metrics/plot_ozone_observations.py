import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from scipy import fftpack
from read_sunspots import *
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Clean up
plt.close('all')

# Load functions to read NOAA data
#execfile("../BrXplo/read_station_data.py")
#from read_station_data import read_station_data_noaa

src_barrow = os.environ['DATA']+'/astra_data/observations/ozone/Barrow/*'
src_prestebakke = os.environ['DATA']+'/astra_data/observations/ozone/Prestebakke/*'
src_jergul = os.environ['DATA']+'/astra_data/observations/ozone/Jergul/NO0030R.*ozone*.nas'
src_karasjok = os.environ['DATA']+'/astra_data/observations/ozone/Karasjok/NO0055R.*ozone*.nas'

try:
    data_barrow
except NameError:
    data_barrow = []
    
    for file in sorted(glob.glob(src_barrow)):
        if int(file[-4:]) < 2003:
            print(file)
            tmp = (read_station_data_noaa(file, utc=-9, start_data=28))
        elif int(file[-4:]) < 2012:
            print(file)
            tmp = (read_station_data_noaa(file, utc=-9, station='other', column=5))
        data_barrow.append(tmp)
    # Concatenate the lists
    data_barrow = pd.concat(data_barrow)

try:
    data_prestebakke
except NameError:
    data_prestebakke = []
    for file in sorted(glob.glob(src_prestebakke)):
        print("Reading file %s" % (file))
        tmp = read_station_data_ebas(file)
        data_prestebakke.append(tmp['O3'])   
   
    # Concatenate the lists
    data_prestebakke = pd.concat(data_prestebakke)
try:
    data_jergul
except NameError:
    data_jergul = []
    data_karasjok = []
    for file in sorted(glob.glob(src_jergul)):
        print("Reading file %s" % (file))
        tmp = read_station_data_ebas(file)
        data_jergul.append(tmp)   
    for file in sorted(glob.glob(src_karasjok)):
        print("Reading file %s" % (file))
        tmp = read_station_data_ebas(file)
        data_karasjok.append(tmp)
    data_jergkara = pd.concat((pd.concat(data_jergul),pd.concat(data_karasjok)))['O3']

    # Round time index to full hour
    data_jergkara.index = data_jergkara.index.round("h")
    data_barrow.index = data_barrow.index.round("h")
    data_prestebakke.index = data_prestebakke.index.round("h")
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



