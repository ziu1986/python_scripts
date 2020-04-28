import os, glob, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.station_info import station_location

# Clean up
plt.close('all')

# Source
src = "growing_season_senorge_1960-2019.txt"

# Read data
data = pd.read_fwf(src)
data.index = data['year']
data = data.drop(columns=['year'])

# Fits
fit = np.polyfit(data.index.values, data['Glength'].values, 1)
fit_function = np.poly1d(fit)

fit_2 = np.polyfit(data.index.values, data['Gbegin'].values, 1)
fit_2_function = np.poly1d(fit_2)

fit_3 = np.polyfit(data.index.values, data['Gend'].values, 1)
fit_3_function = np.poly1d(fit_3)


# Plot it
fig1 = plt.figure(1, figsize=(12,10))
fig1.canvas.set_window_title("greening_season_change_%s" % "Svanvik")

ax11 = plt.subplot(311)

if np.sign(fit[1]) > 0:
    fit_label = "$f(x)=%1.2f\cdot x + %3.2f$" % (fit[0], fit[1])
else:
    fit_label = "$f(x)=%1.2f\cdot x %3.2f$" % (fit[0], fit[1])
data['Glength'].plot(ax=ax11, ls='none', marker='x', color='blue', label="SeNorge")
ax11.plot(data.index, fit_function(data.index.values), color='red', label=fit_label)

ylim = np.ceil(data['Glength'].std())*3.5

ax11.set_xlabel("")
ax11.set_ylabel("G$_{length}$ (days)")
ax11.set_ylim(data['Glength'].mean()-ylim, data['Glength'].mean()+ylim)
ax11.set_xticklabels("")
ax11.axhline(data['Glength'].mean(), ls="--", color='grey')
ax11.axhspan(data['Glength'].mean()-data['Glength'].std(), data['Glength'].mean()+data['Glength'].std(), color='grey', alpha=0.25)
ax11.set_title("")

ax12 = plt.subplot(312)
if np.sign(fit_2[1]) > 0:
    fit_label = "$f(x)=%1.2f\cdot x + %3.2f$" % (fit_2[0], fit_2[1])
else:
    fit_label = "$f(x)=%1.2f\cdot x %3.2f$" % (fit_2[0], fit_2[1])
ax12.plot(data['Glength'].index.values, data['Gbegin'].values, ls='none', marker='x', color='blue', label="SeNorge")
ax12.plot(data['Glength'].index.values, fit_2_function(data['Glength'].index.values), color='red', label=fit_label)

ylim = np.ceil(data['Gbegin'].std())*3.5

ax12.set_xlabel("")
ax12.set_ylabel("G$_{begin}$ (doy)")
ax12.set_ylim(np.ceil(data['Gbegin'].values.mean())-ylim, np.ceil(data['Gbegin'].values.mean())+ylim)

ax12.set_title("")
ax12.set_xticklabels("")
ax12.axhline(data['Gbegin'].mean(), ls="--", color='grey')
ax12.axhspan(data['Gbegin'].mean()-data['Gbegin'].std(), data['Gbegin'].mean()+data['Gbegin'].std(), color='grey', alpha=0.25)  

ax13 = plt.subplot(313)
if np.sign(fit_3[1]) > 0:
    fit_label = "$f(x)=%1.2f\cdot x + %3.2f$" % (fit_3[0], fit_3[1])
else:
    fit_label = "$f(x)=%1.2f\cdot x %3.2f$" % (fit_3[0], fit_3[1])
ax13.plot(data['Glength'].index.values, data['Gend'].values, ls='none', marker='x', color='blue', label="SeNorge")
ax13.plot(data['Glength'].index.values, fit_3_function(data['Glength'].index.values), color='red', label=fit_label)

ylim = np.ceil(data['Gend'].std())*3.5

ax13.set_xlabel("Time (years)")
ax13.set_ylabel("G$_{end}$ (doy)")
ax13.set_ylim(np.ceil(data['Gend'].values.mean())-ylim, np.ceil(data['Gend'].values.mean())+ylim)
ax13.axhline(data['Gend'].mean(), ls="--", color='grey')
ax13.axhspan(data['Gend'].mean()-data['Gend'].std(), data['Gend'].mean()+data['Gend'].std(), color='grey', alpha=0.25)

for ax in fig1.axes:
    ax.legend(loc='upper left')
# Show it
plt.show(block=False)

