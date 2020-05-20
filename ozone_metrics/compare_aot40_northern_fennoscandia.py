import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from matplotlib.dates import date2num # Convert dates to matplotlib axis coords
from matplotlib import dates
from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.ozone_tools import *
from mytools.station_info import station_location

# Clean up
plt.close('all')

# Read data
src = os.environ['DATA']+'/astra_data/observations/ozone/'
src_svanvik_OzoNorClim = os.environ['DATA']+'/astra_data/input_data/DO3SE_input/svanvik_ozone_20*.csv'
src_stations = ('Esrange', 'Jergul', 'Karasjok', 'Pallas', 'Prestebakke', 'Svanvik')

try:
    data
except NameError:
    data = {}
    for station in src_stations:
        data.update({station:load_data(src+station+'/*.nas')})

    # Concate Jergul and Karasjok data
    data_jergkara = pd.concat((data['Jergul'], data['Karasjok']))

    # Read and convert file data
    data_svanvik_OzoNorClim = []
    for file in sorted(glob.glob(src_svanvik_OzoNorClim)):
        tmp_data_svanvik = pd.read_csv(file)
        data_svanvik_OzoNorClim.append(tmp_data_svanvik)
    # Concat data Svanvik data
    data_svanvik_OzoNorClim = pd.concat(data_svanvik_OzoNorClim)
    data_svanvik_OzoNorClim = data_svanvik_OzoNorClim.set_index('Year_Mnth_Date_Time(NMT)')
    data_svanvik_OzoNorClim.index = pd.to_datetime(data_svanvik_OzoNorClim.index)

# Plot it
fig11 = plt.figure(11, figsize=(12,10))
fig11.canvas.set_window_title("ozone_fennoscandic_obs_aot40")

ax111 = plt.subplot(211)
ax111.set_title("(a)") # 3 months
aot40_prestebakke = compute_aot(data['Prestebakke'], time_start=1, time_end=23, month_start=5, month_end=7)
aot40_jerkara = compute_aot(data_jergkara, time_start=1, time_end=23, month_start=5, month_end=7)
aot40_pallas = compute_aot(data['Pallas'], time_start=1, time_end=23, month_start=5, month_end=7)
aot40_esrange = compute_aot(data['Esrange'], time_start=1, time_end=23, month_start=5, month_end=7)
aot40_svanvik = compute_aot(pd.concat((data['Svanvik'],data_svanvik_OzoNorClim.iloc[:,0])), time_start=1, time_end=23, month_start=5, month_end=9)

aot40 = pd.DataFrame({"Esrange":aot40_esrange, "Pallas":aot40_pallas, "Jergul/Karasjok":aot40_jerkara, "Svanvik":aot40_svanvik}) # "Prestebakke":aot40_prestebakke,

aot40.plot.bar(ax=ax111, width=1., color={'blue', 'orange', 'black', 'blueviolet'}) # 'red', 
ax111.set_ylim(0,8000)
ax111.legend(ncol=2)
ax111.axhline(3000, ls=':', color='red')

ax112 = plt.subplot(212)
ax112.set_title("(b)") # 6 months
aot40_prestebakke = compute_aot(data['Prestebakke'], time_start=1, time_end=23, month_start=4, month_end=9)
aot40_jerkara = compute_aot(data_jergkara, time_start=1, time_end=23, month_start=4, month_end=9)
aot40_pallas = compute_aot(data['Pallas'], time_start=1, time_end=23, month_start=4, month_end=9)
aot40_esrange = compute_aot(data['Esrange'], time_start=1, time_end=23, month_start=4, month_end=9)
aot40_svanvik = compute_aot(pd.concat((data['Svanvik'],data_svanvik_OzoNorClim.iloc[:,0])), time_start=1, time_end=23, month_start=4, month_end=9)

aot40 = pd.DataFrame({"Esrange":aot40_esrange, "Pallas":aot40_pallas, "Jergul/Karasjok":aot40_jerkara, "Svanvik":aot40_svanvik}) # "Prestebakke":aot40_prestebakke,

aot40.plot.bar(ax=ax112, width=1., color={'blue', 'orange', 'black', 'blueviolet'}) # 'red', 
ax112.set_ylim(0,16000)
ax112.legend(ncol=2)
ax112.axhline(5000, ls=':', color='red')

ax111.set_xlabel("")
ax111.set_ylabel("")
ax111.set_xticklabels("")
ax112.set_xlabel("Time (years)")
ax112.set_ylabel("AOT40 (ppb h)", y=1)



# Show it
plt.show(block=False)
