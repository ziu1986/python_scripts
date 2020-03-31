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
src_rad_svanvik_clim = os.environ['DATA']+'/astra_data/observations/svanvik_glob_rad_2000_2017.csv'
src_rad_svanvik = os.environ['DATA']+'/astra_data/observations/svanvik_glob_rad_2018_2019.csv'
# Load data
try:
    data_svanvik
except NameError:
    
    data_svanvik = pd.read_csv(src_rad_svanvik)
    data_svanvik_clim = pd.read_csv(src_rad_svanvik_clim)
    data_svanvik.index = pd.date_range("%s" % data_svanvik['Time measured'][0][:10], "%s" % data_svanvik['Time measured'].iloc[-1][:-3], freq='H')
    data_svanvik_clim.index = pd.date_range("%s" % data_svanvik_clim['Time measured'][0][:10], "%s" % data_svanvik_clim['Time measured'].iloc[-1][:-3], freq='H')
    data_svanvik = data_svanvik.drop(columns=['Time measured'])
    data_svanvik_clim = data_svanvik_clim.drop(columns=['Time measured'])

    data_svanvik_clim.loc[:,'hour'] = data_svanvik_clim.index.hour.values
    data_svanvik_clim.loc[:,'day'] = data_svanvik_clim.index.day.values
    data_svanvik_clim.loc[:,'month'] = data_svanvik_clim.index.month.values

    data_svanvik.loc[:,'hour'] = data_svanvik.index.hour.values
    data_svanvik.loc[:,'day'] = data_svanvik.index.day.values
    data_svanvik.loc[:,'month'] = data_svanvik.index.month.values
    
# Plot it
fig1 = plt.figure(1)
ax11 = plt.subplot(211)
ax11.set_title("(a)")
xtime_hour = np.arange(1, data_svanvik_clim.groupby(['month','day','hour']).max().size+1)
ax11.plot(xtime_hour, data_svanvik_clim.groupby(['month','day','hour']).max(), color='red', label='max', alpha=0.5)
ax11.plot(xtime_hour, data_svanvik_clim.groupby(['month','day','hour']).mean(), color='black', label='mean', ls=':', alpha=0.75)
ax11.plot(xtime_hour, data_svanvik_clim.groupby(['month','day','hour']).min(), color='blue', label='min', ls='-', alpha=1)

ax12 = plt.subplot(223)
ax12.set_title("(b)")

ax12.plot(xtime_hour, (data_svanvik['2018'].groupby(['month','day','hour']).max()-data_svanvik_clim.groupby(['month','day','hour']).max()), color='red', ls='-', alpha=0.5)
ax12.plot(xtime_hour, (data_svanvik['2018'].groupby(['month','day','hour']).min()-data_svanvik_clim.groupby(['month','day','hour']).min()), color='blue', ls='-')
ax12.plot(xtime_hour, (data_svanvik['2018'].groupby(['month','day','hour']).mean()-data_svanvik_clim.groupby(['month','day','hour']).mean()), color='black', ls=':', alpha=0.75)

ax13 = plt.subplot(224)
ax13.set_title("(c)")

ax13.plot(xtime_hour, (data_svanvik['2019'].groupby(['month','day','hour']).max()-data_svanvik_clim.groupby(['month','day','hour']).max()), color='red', ls='-', alpha=0.5)
ax13.plot(xtime_hour, (data_svanvik['2019'].groupby(['month','day','hour']).min()-data_svanvik_clim.groupby(['month','day','hour']).min()), color='blue', ls='-')
ax13.plot(xtime_hour, (data_svanvik['2019'].groupby(['month','day','hour']).mean()-data_svanvik_clim.groupby(['month','day','hour']).mean()), color='black', ls=':', alpha=0.75)

for ax in fig1.axes:
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("")
    ybound = ax.get_ybound()[1]-50
    ax.text(1*24, ybound, "Jan")
    ax.text((31+1)*24, ybound,"Feb")
    ax.text((31+29+1)*24, ybound,"Mar")
    ax.text((31+29+31+1)*24, ybound,"Apr")
    ax.text((31+29+31+30+1)*24, ybound,"May")
    ax.text((31+29+31+30+31+1)*24, ybound, "Jun")
    ax.text((31+29+31+30+31+30+1)*24, ybound, "Jul")
    ax.text((31+29+31+30+31+30+31+1)*24, ybound, "Aug")
    ax.text((31+29+31+30+31+30+31+31+1)*24, ybound, "Sep")
    ax.text((31+29+31+30+31+30+31+31+30+1)*24, ybound, "Oct")
    ax.text((31+29+31+30+31+30+31+31+30+31+1)*24, ybound, "Nov")
    ax.text((31+29+31+30+31+30+31+31+30+31+30+1)*24, ybound, "Dec")
ax11.set_ylabel("Global radiation (W$\,m^{-2}s^{-1}$)")
#ax11.set_xlabel("")
ax12.set_ylabel("$\Delta$Global radiation (W$\,m^{-2}s^{-1}$)")
ax11.legend()

#print("clim %")
for iyear in (2018, 2019):
    print('%d %1.2f' % (iyear,(data_svanvik['%d' % iyear].groupby(['month','day','hour']).max()-data_svanvik_clim.groupby(['month','day','hour']).mean()).dropna().mean()))


# Show it
plt.show(block=False)

