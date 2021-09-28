import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
#import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
#from matplotlib.dates import date2num # Convert dates to matplotlib axis coords
#from matplotlib import dates
from scipy import fftpack
from scipy import stats
#from read_sunspots import *
#from mytools.met_tools import *
#from mytools.netcdf_tools import *
from mytools.ozone_tools import *
#from mytools.cartopy_tools import scale_bar
from mytools.station_info import station_location

# Clean up
plt.close('all')

#---------------------------------------------------------------------------------------------------------------------------------
# Read data
src = os.environ['DATA']+'/astra/observations/ozone/'
src_stations = ('Barrow', 'Esrange', 'Jergul', 'Karasjok', 'Pallas', 'Prestebakke')

try:
    data
except NameError:
    data = {}
    for station in src_stations:
        if station=='Barrow':
            data.update({station:load_data(src+station+'/*', type="Barrow")})
        else:
            data.update({station:load_data(src+station+'/*.nas')})

    # Concate Jergul and Karasjok data
    data_jergkara = pd.concat((data['Jergul'], data['Karasjok']))

# Analysis
try:
    fft
except NameError:
    fft = {}
    fft_freq = {}
    for each in data:
        tmp = fftpack.fft(data[each].resample('1D').mean().fillna(method='ffill').fillna(method='bfill'))
        #tmp = fftpack.fft(data[each].resample('1M').mean().fillna(method='ffill'))
        fft.update({each:tmp})
        num = len(tmp) 
        fft_freq.update({each:fftpack.fftfreq(num)})

# Plot
fig1 = plt.figure(1,figsize=(16,9))
for i, istat in zip(np.arange(1,6),src_stations):
    ax11 = plt.subplot(5,1,i)
    ax11.set_title("%s" % (istat))
    ax11.stem(1/fft_freq[istat][1:], np.abs(fft[istat][1:]))
    

for ax in fig1.axes:
    ax.set_xlim(0,30)
    ax.set_ylim(0,3600)
    ax.set_xlabel('Frequency (days)')
    ax.set_ylabel('Power')
    ax.axvspan(27, 28, color='red', alpha=0.5, linewidth=4)

plt.show(block=False)


