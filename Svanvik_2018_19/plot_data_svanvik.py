import os, glob, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
from scipy.constants import *     # Get physics constants
#from station_info import station_location
from mytools.plot_tools import *

# Clean up
plt.close('all')

# Data sources
src_svanvik = os.environ['DATA']+'/astra_data/observations/ozone/Svanvik/NO0047R.*ozone*.xls'
src_svanvik_precip = os.environ['DATA']+'/astra_data/observations/svanvik_accu_precip_2009-2020.xlsx'

data_list = []
# Read data
try:
    data_svanvik
except NameError:
    for file in sorted(glob.glob(src_svanvik)):
        data_list.append(pd.read_excel(file, index_col=0, header=0))
    data_svanvik = pd.concat(data_list)
    data_svanvik_ozone = data_svanvik['O3_mugm-3'].where(data_svanvik['O3_mugm-3']>=0).dropna()/2.

    data_svanvik_temp = (data_svanvik['Temp_degC'].where((data_svanvik['Temp_degC']>-50) & (data_svanvik['Temp_degC']<50)).dropna())

    data_svanvik_precip = pd.read_excel(src_svanvik_precip).dropna()
    data_svanvik_precip[u"sum(precipitation_amount P1D)"] = [float(each.replace(',','.')) for each in data_svanvik_precip[u"sum(precipitation_amount P1D)"].values]
    data_svanvik_precip.index = [np.datetime64("%s-%s-%s" % (itime[-4:], itime[-7:-5], itime[:2])) for itime in data_svanvik_precip['time'].dropna().values ]
    data_svanvik_precip = data_svanvik_precip[u"sum(precipitation_amount P1D)"]



# Plot it
fig1 = plt.figure(1, figsize=(16,12))
fig1.canvas.set_window_title("data_svanvik_2018_2019")

ax11 = plt.subplot(311)
ax11.set_title('(a)')
data_svanvik_ozone.plot(ax=ax11, color='blueviolet', ls='None', marker='+')

ax11.set_ylabel("$[O_3]$ (ppb)")
ax11.set_xticklabels("")
ax11.axvspan(dt.date(2018, 1, 1), data_svanvik_ozone.index[0], color='grey', alpha=0.5)
ax11.axvspan(data_svanvik_ozone['2018'].index[-1], data_svanvik_ozone['2019'].index[0], color='grey', alpha=0.5)
ax11.axvspan(data_svanvik_ozone.index[-1], dt.date(2020, 1, 1), color='grey', alpha=0.5)
ax11.axvspan(dt.date(2018, 7, 10), dt.date(2018, 7, 22), color='grey', alpha=0.5)

ax12 = plt.subplot(312)
ax12.set_title('(b)')
data_svanvik_temp.plot(ax=ax12, color='red')
ax12.set_xticklabels("")
ax12.set_ylabel("T ($^\circ$C)")

ax13 = plt.subplot(313)
ax13.set_title('(c)')
data_svanvik_precip['2018':'2019'].plot(ax=ax13, color='blue')
ax13.set_ylabel('$\Sigma_{day}$Precip (mm)')


for ax in fig1.axes:
    ax.set_xlabel("")
    ax.set_xlim(dt.date(2018, 1, 1), dt.date(2020, 1, 1))
    ax.axvspan(dt.date(2018, 1, 1), data_svanvik_ozone.index[0], facecolor='None', edgecolor='black', hatch='//', alpha=0.5)
    ax.axvspan(data_svanvik_ozone['2018'].index[-1], data_svanvik_ozone['2019'].index[0], facecolor='None', edgecolor='black', hatch='//', alpha=0.5)
    ax.axvspan(data_svanvik_ozone.index[-1], dt.date(2020, 1, 1), facecolor='None', edgecolor='black',  hatch='//', alpha=0.5)
    ax.axvspan(dt.date(2018, 7, 10), dt.date(2018, 7, 22), facecolor='None', edgecolor='black',  hatch='//', alpha=0.5)

ax13.set_xlabel("Time (months)")
# Show it
plt.show(block=False)
