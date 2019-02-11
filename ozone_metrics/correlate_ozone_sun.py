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

# Load modules and EBAS reading routine
execfile('read_ebas.py')

# The data
nc_src = os.environ['DATA']+'/astra_data/observations/sunspots/SN_d_tot_V2.0.txt'
src_jergul = os.environ['DATA']+'/astra_data/observations/ozone/Jergul/NO0030R.*ozone*.nas'
src_karasjok = os.environ['DATA']+'/astra_data/observations/ozone/Karasjok/NO0055R.*ozone*.nas'

station_location = {"Jergul":(69.45,24.6),"Karasjok":(69.467,25.217),"Svanvik":(69.45,30.03)}

# Loop through EBAS data and transform them to pandas timeseries
try:
    data_jergul
except NameError:
    data_jergul = []
    data_karasjok = []
    for file in sorted(glob.glob(src_jergul)):
        print("Reading file %s" % (file))
        tmp = read_station_data(file)
        data_jergul.append(pd.Series(tmp['O3'],index=tmp['time']))   
    for file in sorted(glob.glob(src_karasjok)):
        print("Reading file %s" % (file))
        tmp = read_station_data(file)
        data_karasjok.append(pd.Series(tmp['O3'],index=tmp['time']))
    # Concatenate the lists
    data_jergul = pd.concat(data_jergul)
    data_karasjok = pd.concat(data_karasjok)
    data_jerg_kara = pd.concat((data_jergul, data_karasjok))
    # Sunspot data
    data_ss = read_sunspots(nc_src)


# Plot it
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot()

data_jerg_kara.plot(ax=ax11, color='orange', ls='None', marker='+', label='Jergul/Karasjok')
(data_ss['1988-04-18':'2010-02-28']['Ntot']/data_ss['1988-04-18':'2010-02-28'].max()['Ntot']*80).plot(ax=ax11, marker='x', alpha=0.25, color='black', label='sunspots')

ax11.legend()

# Show it
plt.show(block=False)
