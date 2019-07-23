'''
Look at change in greening season begin and length for verious stations.
'''
import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

# In situ data from dellbart et al. (2006)
aksarka_year = np.array([1982,] + range(1984,2003))
aksarka_gstart = (np.array((770.15,776.15,775.95,772.56,784.47,776.94,761.28,762.67,764.85,772.78,764.85,764.85,742.25,770.20,764.25,783.09,783.48,754.54, 766.43,760.49))*0.747-408.13).astype(int)

# Clean up
plt.close('all')
selection = -1
labels = ('Sodankyla','Karasjok','Ulborg','Castel Porziano','Auchencorth Moss','Hyytiala','Harvard Forest','Citrus Orchard','Blodgett Forest', 'Svanvik', 'Aksarka')
station_lat = (67,69,56,41,55,62,42,36,39,69.45,66.55)
station_lon = (27,24,8,12,360-3,24,360-77,360-120,360-120,30,67.79)
# Data source
nc_src = os.environ['DATA']+'/astra_data/input_data/ctm_input/GROWING_SEASON/*.nc'
glen_data = []
gday_data = []
try:
    data2
except NameError:
    for date in sorted(glob.glob(nc_src)):
        print(date)
        data = xr.open_dataset(date)
        glen_data.append(data["GLEN"].sel(lat=station_lat[selection],method='nearest').sel(lon=station_lon[selection],method='nearest'))
        gday_data.append(data["GDAY"].sel(lat=station_lat[selection],method='nearest').sel(lon=station_lon[selection],method='nearest'))

    data_glen = xr.concat(glen_data, dim='year')
    data_glen.coords['year'] = np.arange(1990,2018)
    data_gday = xr.concat(gday_data, dim='time')

gbegin = data_gday.where(data_gday==1, drop=True).time.dt.dayofyear
gend = data_gday.groupby(data_gday.time.dt.year).max().data+data_gday.where(data_gday==1, drop=True).time.dt.dayofyear

fit = np.polyfit(data_glen['year'], data_glen, 1)
fit_function = np.poly1d(fit)

fit_2 = np.polyfit(data_gday.where(data_gday==1, drop=True).time.dt.year, gbegin, 1)
fit_2_function = np.poly1d(fit_2)

fit_3 = np.polyfit(data_gday.where(data_gday==1, drop=True).time.dt.year, gend, 1)
fit_3_function = np.poly1d(fit_3)

# Plot it
fig1 = plt.figure(1,figsize=(12,10))
fig1.canvas.set_window_title("greening_season_change_%s" % labels[selection])

ax11 = plt.subplot(311)

if np.sign(fit[1]) > 0:
    fit_label = "$f(x)=%1.2f\cdot x + %3.2f$" % (fit[0], fit[1])
else:
    fit_label = "$f(x)=%1.2f\cdot x %3.2f$" % (fit[0], fit[1])
data_glen.plot(ax=ax11, ls='none', marker='x', label="OpenIFS")
ax11.plot(data_glen.coords['year'], fit_function(data_glen.coords['year']), label=fit_label)
ax11.set_xlabel("")
ax11.set_ylabel("G$_{length}$ (days)")
ax11.set_ylim(data_glen.mean()-60,data_glen.mean()+60)
ax11.set_xticklabels("")
ax11.axhline(data_glen.mean(),ls="--",color='grey')
ax11.axhspan(data_glen.mean()-data_glen.std(), data_glen.mean()+data_glen.std(),color='grey',alpha=0.25)
ax11.set_title(labels[selection])

ax12 = plt.subplot(312)
if np.sign(fit_2[1]) > 0:
    fit_label = "$f(x)=%1.2f\cdot x + %3.2f$" % (fit_2[0], fit_2[1])
else:
    fit_label = "$f(x)=%1.2f\cdot x %3.2f$" % (fit_2[0], fit_2[1])
ax12.plot(data_glen.coords['year'], gbegin.data, ls='none', marker='x', label="OpenIFS")
ax12.plot(data_glen.coords['year'], fit_2_function(data_glen.coords['year']), label=fit_label)
ax12.set_xlabel("")
ax12.set_ylabel("G$_{begin}$ (doy)")
ax12.set_ylim(round(gbegin.data.mean())-30,round(gbegin.data.mean())+30)

if(labels[selection]=="Aksarka"):
    ax12.plot(aksarka_year[7:], aksarka_gstart[7:], "d", color='black', label="in-situ")

ax12.set_title("")
ax12.set_xticklabels("")
ax12.axhline(gbegin.mean(),ls="--",color='grey')
ax12.axhspan(gbegin.mean()-gbegin.std(), gbegin.mean()+gbegin.std(),color='grey',alpha=0.25)  

ax13 = plt.subplot(313)
if np.sign(fit_3[1]) > 0:
    fit_label = "$f(x)=%1.2f\cdot x + %3.2f$" % (fit_3[0], fit_3[1])
else:
    fit_label = "$f(x)=%1.2f\cdot x %3.2f$" % (fit_3[0], fit_3[1])
ax13.plot(data_glen.coords['year'], gend.data, ls='none', marker='x', label="OpenIFS")
ax13.plot(data_glen.coords['year'], fit_3_function(data_glen.coords['year']), label=fit_label)
ax13.set_xlabel("Time (years)")
ax13.set_ylabel("G$_{end}$ (doy)")
ax13.set_ylim(round(gend.data.mean())-30,round(gend.data.mean())+30)
ax13.axhline(gend.mean(),ls="--",color='grey')
ax13.axhspan(gend.mean()-gend.std(), gend.mean()+gend.std(),color='grey',alpha=0.25)

for ax in fig1.axes:
    ax.legend(loc='upper left')
# Show it
plt.show(block=False)
