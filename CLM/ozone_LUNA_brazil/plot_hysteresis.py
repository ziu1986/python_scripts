import os, sys, glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all, plot_error_bands, get_month_name, seconds_in_month
from mytools.ozone_tools import flunder
from mytools.clm_tools import *

# Clean up
plt.close("all")

# CLM simulations
# Reference - no ozone - ozone_luna_100
src = os.environ['CESM_RUN'] + '/archive/'
exp = ('brazil_2000',
       'brazil_2000_5.0.34_ozone_luna_100_wos',
       'brazil_2000_5.0.34_ozone_luna_80_wos',
       'brazil_2000_5.0.34_ozone_luna_60_wos',
       'brazil_2000_5.0.34_ozone_luna_40_wos',
       'brazil_2000_5.0.34_ozone_luna_20_wos',
       'brazil_2000_ozone_luna_100',
       'brazil_2000_ozone_luna_80',
       'brazil_2000_ozone_luna_60',
       'brazil_2000_ozone_luna_40',
       'brazil_2000_ozone_luna_20')

data = []

for iexp in exp:
    data.append(load_data(src + iexp + '/lnd/hist/*.nc', var=['GPP']))

# Plot it
marker = ('x', '+','o', 's', 'd', '*')
fig1 = plt.figure(1, figsize=(16,8))
ax11 = plt.subplot(321)
ax12 = plt.subplot(322)
ax13 = plt.subplot(323)
ax14 = plt.subplot(324)
ax15 = plt.subplot(325)

for i, imarker, iexp in zip((1,6), marker, ('spin-up no O3', 'spin-up 100ppb')):
    ax11.plot(data[i]['GPP']*1e6, data[0]['GPP']*1e6, ls='None', marker=imarker, label=iexp)
    ax12.plot(data[i+1]['GPP']*1e6, data[0]['GPP']*1e6, ls='None', marker=imarker)
    ax13.plot(data[i+2]['GPP']*1e6, data[0]['GPP']*1e6, ls='None', marker=imarker)
    ax14.plot(data[i+3]['GPP']*1e6, data[0]['GPP']*1e6, ls='None', marker=imarker)
    ax15.plot(data[i+4]['GPP']*1e6, data[0]['GPP']*1e6, ls='None', marker=imarker)


for ax in fig1.axes:
    ax.set_xlim(0,110)
    ax.set_ylim(0,110)
    ax.set_ylabel("GPP$_{exp}$")
    ax.set_xlabel("GPP$_{ref}$")
    ax.plot((0,110), (0,110), ls=':', color='grey')
    ax.legend()

ax11.set_title("100 ppb", x=0.8, y=0.88)
ax12.set_title("80 ppb", x=0.8, y=0.88)
ax13.set_title("60 ppb", x=0.8, y=0.88)
ax14.set_title("40 ppb", x=0.8, y=0.88)
ax15.set_title("20 ppb", x=0.8, y=0.88)

fig2 = plt.figure(2, figsize=(16,8))

ax21 = plt.subplot(321)
ax22 = plt.subplot(322)
ax23 = plt.subplot(323)
ax24 = plt.subplot(324)
ax25 = plt.subplot(325)

for i, imarker, iexp in zip((1,6), marker, ('spin-up no O3', 'spin-up 100ppb')):
    ((data[i]['GPP']-data[0]['GPP'])/data[0]['GPP']*100).plot(ax=ax21, ls='None', marker=imarker, label=iexp)
    ((data[i+1]['GPP']-data[0]['GPP'])/data[0]['GPP']*100).plot(ax=ax22, ls='None', marker=imarker)
    ((data[i+2]['GPP']-data[0]['GPP'])/data[0]['GPP']*100).plot(ax=ax23, ls='None', marker=imarker)
    ((data[i+3]['GPP']-data[0]['GPP'])/data[0]['GPP']*100).plot(ax=ax24, ls='None', marker=imarker)
    ((data[i+4]['GPP']-data[0]['GPP'])/data[0]['GPP']*100).plot(ax=ax25, ls='None', marker=imarker)

for ax in fig2.axes:
    ax.set_ylim(-15, 15)
    ax.set_ylabel("$\Delta GPP_{exp-ref}$ (%)")
    ax.set_xlabel("Time (years)")
    ax.legend()

ax21.set_title("100 ppb", x=0.88, y=0.88)
ax22.set_title("80 ppb", x=0.88, y=0.88)
ax23.set_title("60 ppb", x=0.88, y=0.88)
ax24.set_title("40 ppb", x=0.88, y=0.88)
ax25.set_title("20 ppb", x=0.88, y=0.88)

def func(x, limit, start, tau):
    return(limit-(limit-start)*np.exp(-tau*x))
def timedelta_to_days(timedelta):
    days = timedelta.astype('timedelta64[D]')
    return(days / np.timedelta64(1,'D'))

x = np.arange(len(data[0].time))
par0 = [(-4.7, 0, 1/(365.*3)), (-2.7, -7.8, 1/(365.*5))]
#ax25.plot(data[0].time, func(x,*par0[0]))
#ax25.plot(data[0].time, func(x,*par0[1]))

from scipy.optimize import curve_fit

for i in np.arange(1,6):
    test = ((data[i]-data[0])/data[0]).dropna(dim='time')*100
    vals, covar = curve_fit(func, timedelta_to_days(test.time-test.time[0]).values, flunder(test['GPP'].values), p0=par0[0])
    fig2.axes[i-1].plot(data[0].time, func(x,*vals), ls='--', color='black', label='fit spin-up no O3')
    print(exp[i], vals, covar)

    test = ((data[i+5]-data[0])/data[0]).dropna(dim='time')*100
    vals, covar = curve_fit(func, timedelta_to_days(test.time-test.time[0]).values, flunder(test['GPP'].values), p0=par0[1])
    fig2.axes[i-1].plot(data[0].time, func(x,*vals), ls=':', color='black', label='fit spin-up 100ppb')
    print(exp[i+5], vals, covar)

# Show it
plt.show(block=False)
