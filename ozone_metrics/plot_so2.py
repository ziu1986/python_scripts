import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
execfile('read_ebas.py')
# Clean up
plt.close('all')

# Directories of data
src_jergul = os.environ['DATA']+'/processed_data/observations/Jergul/NO0030R.*pack*.nas'
src_karasjok = os.environ['DATA']+'/processed_data/observations/Karasjok/NO0055R.*pack*.nas'
src_svanvik = os.environ['DATA']+'/processed_data/observations/Svanvik/NO0047R.*pack*.nas'
src = os.environ['DATA']+'/processed_data/CTM3_oivind/'
nc_src = 'osloctm_so2*.nc'
# Loop through EBAS data and transform them to pandas timeseries
data_jergul = []
data_karasjok = []
data_svanvik = []
try:
    data
except NameError:
    data = read_data(src+nc_src,var='SO2',lev=(1,1))

for file in sorted(glob.glob(src_svanvik)):
    tmp = read_station_data(file, tracer=('SO4','SO2'))
    data_svanvik.append(pd.Series(tmp['SO2'], index=tmp['time']))
    
for file in sorted(glob.glob(src_jergul)):
    tmp = read_station_data(file, tracer=('SO2',))
    data_jergul.append(pd.Series(tmp['SO2'],index=tmp['time']))
for file in sorted(glob.glob(src_karasjok)):
    tmp = read_station_data(file, tracer=('SO2',))
    data_karasjok.append(pd.Series(tmp['SO2'],index=tmp['time']))
# Concatenate the lists
data_jergul = pd.concat(data_jergul)
data_karasjok = pd.concat(data_karasjok)
data_svanvik = pd.concat(data_svanvik)

# Plotting
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("so2_ebas_timeseries")
ax11 = plt.subplot()
data_jergul.plot(marker='x', ls='none', color='grey', alpha=0.15, label='Jergul')
data_karasjok.plot(marker='+', ls='none', color='grey', alpha=0.15, label='Karasjok')
#data_svanvik.plot(marker='', ls='-', color='blueviolet', alpha=0.15, label='Svanvik')
(data.sel(lat=69,method='nearest').sel(lon=24,method='nearest')/1e-9).plot(color='blue', alpha=0.15, label='OsloCTM3')
ax11.set_xlabel('Time')
ax11.set_ylabel('$SO_2$ (ppb)')
ax11.legend()
# Show it
plt.show(block=False)
