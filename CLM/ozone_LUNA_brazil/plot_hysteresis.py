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
exp = ('brazil_2000_5.0.34',
       'brazil_2000_5.0.34_ozone_luna_100_wos',
       'brazil_2000_5.0.34_ozone_luna_80_wos',
       'brazil_2000_5.0.34_ozone_luna_60_wos',
       'brazil_2000_5.0.34_ozone_luna_40_wos',
       'brazil_2000_5.0.34_ozone_luna_20_wos',
       'brazil_2000_5.0.34_ozone_luna_100',
       'brazil_2000_5.0.34_ozone_luna_80',
       'brazil_2000_5.0.34_ozone_luna_60',
       'brazil_2000_5.0.34_ozone_luna_40',
       'brazil_2000_5.0.34_ozone_luna_20',
       'brazil_2000_5.0.34_ozone_luna_0')

data = []
var = 'TOTVEGN'

for iexp in exp:
    data.append(load_data(src + iexp + '/lnd/hist/*.nc', var=[var]))

# Plot it
marker = ('x', '+','o', 's', 'd', '*')

fig2 = plt.figure(2, figsize=(16,8))
fig2.canvas.set_window_title("hysteresis_%s_ref_exp_time" % var.lower())
ax21 = plt.subplot(321)
ax22 = plt.subplot(322)
ax23 = plt.subplot(323)
ax24 = plt.subplot(324)
ax25 = plt.subplot(325)
ax26 = plt.subplot(326)

for i, imarker, iexp in zip((1,6), marker, ('spin-up no O3', 'spin-up 100 ppb')):
    for j, ax in zip(range(5), fig2.axes):
        print(exp[i+j])
        ((data[i+j][var]-data[0][var])/data[0][var]*100).plot(ax=ax, ls='None', marker=imarker, label=iexp)

((data[-1][var]-data[0][var])/data[0][var]*100).plot(ax=ax26, ls='None', marker=marker[1], label='spin-up 100 ppb', color=ax25.lines[1].get_color())
    
for ax in fig2.axes:
    ax.set_ylim(-29, 29) # 29, 15, 4
    ax.set_ylabel("$\Delta %s_{exp-ref}$ (%%)" % var)
    ax.set_xlabel("Time (years)")
    ax.legend(loc='lower left')

ax21.set_title("100 ppb", x=0.88, y=0.88)
ax22.set_title("80 ppb", x=0.88, y=0.88)
ax23.set_title("60 ppb", x=0.88, y=0.88)
ax24.set_title("40 ppb", x=0.88, y=0.88)
ax25.set_title("20 ppb", x=0.88, y=0.88)
ax26.set_title("0 ppb", x=0.88, y=0.88)

b_fit = True

if b_fit:

    def funcExp(x, limit, start, tau):
        return(limit-(limit-start)*np.exp(-tau*x))
    def funcP1d(x, slope, intercept):
        return(slope*x+intercept)
    def timedelta_to_days(timedelta):
        days = timedelta.astype('timedelta64[D]')
        return(days / np.timedelta64(1,'D'))

    x = np.arange(len(data[0].time))
    
    
    from scipy.optimize import curve_fit

    vals_list = {}
    covar_list = {}
    func = funcP1d
    par0 = ((0,0),(0,0))#[(-4.7, 0, 1/(365.*3)), (-2.7, -7.8, 1/(365.*5))]
    # Output
    f = open('fit_results_%s.txt' % var, 'w')

    for i in np.arange(1,6):
        for j, iexp in zip((0,5), ('fit spin-up no O3', 'fit spin-up 100 ppb')):
            test = ((data[i+j]-data[0])/data[0]).dropna(dim='time')*100
            vals, covar = curve_fit(func, timedelta_to_days(test.time-test.time[0]).values, flunder(test[var].values), p0=par0[0])
            fig2.axes[i-1].plot(data[0].time, func(x,*vals), ls='--', color='black', label=iexp)
        
            vals_list.update({exp[i+j]:vals})
            covar_list.update({exp[i+j]:covar})
        
            f.write("%s %s %s" % (exp[i+j], vals, covar))

    test = ((data[-1]-data[0])/data[0]).dropna(dim='time')*100
    vals, covar = curve_fit(func, timedelta_to_days(test.time-test.time[0]).values, flunder(test[var].values), p0=par0[1])
    vals_list.update({exp[-1]:vals})
    covar_list.update({exp[-1]:covar})
    fig2.axes[-1].plot(data[0].time, func(x,*vals), ls=':', color='black', label='fit spin-up 100 ppb')
    f.write("%s %s %s" % (exp[-1], vals, covar))
    # Close file
    f.close()

fig3 = plt.figure(3, figsize=(16,8))
fig3.canvas.set_window_title("hysteresis_%s_ref_exp_time_100_0" % var.lower())
ax31 = plt.subplot()

((data[1][var]-data[0][var])/data[0][var]*100).plot(ax=ax31, ls='None', marker=marker[0], label='spin-up no O3 -> 100 ppb')
((data[-1][var]-data[0][var])/data[0][var]*100).plot(ax=ax31, ls='None', marker=marker[0], label='spin-up 100 ppb -> 0')
if b_fit:
    ax31.plot(data[0].time, func(x,*vals_list[exp[1]]), ls='--', color='black', label='fit spin-up no O3')
    ax31.plot(data[0].time, func(x,*vals_list[exp[-1]]), ls=':', color='black', label='fit spin-up 100 ppb')
ax31.set_ylim(-20,2)
ax31.set_xlabel("Time (years)")
ax31.set_ylabel("$\Delta %s_{exp-ref}$ (%%)" % var)
ax31.legend()

# Show it
plt.show(block=False)
