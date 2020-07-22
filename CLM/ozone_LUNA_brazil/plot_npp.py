import os, sys, glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all, plot_error_bands, get_month_name, seconds_in_month
from clm_tools import *

# Clean up
plt.close("all")

threshold = (0, 0.2, 0.5, 0.8, 1, 2, 3, 4, 5)
ozone = np.arange(20,100,20)

brazil_0_ref = load_data("brazil_0_ref.nc")
brazil = {}

for ioz in ozone:
    for ithr in threshold:
        brazil.update({(ithr,ioz):load_data("brazil_%s_%s.nc" % (ioz, ithr))})

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot()
fig2 = plt.figure(2, figsize=(16,9))
ax21 = plt.subplot()
fig3 = plt.figure(3, figsize=(16,9))
ax31 = plt.subplot(211)
ax32 = plt.subplot(212)


icolor = {20:'lightblue',40:'blue',60:'darkblue',80:'red'}
sec_in_months = [seconds_in_month(each,2019) for each in np.arange(1,13)]
sec_in_months = xr.DataArray(data=sec_in_months, coords={'month':np.arange(1,13)},dims='month')
area = 49089.5*1e6 # m**2

for ioz in ozone:
    color = icolor[ioz]
    alpha = 1
    for ithr in threshold:
        #print(brazil[(ithr,ioz)]['TOTVEGC'])
        if ithr < 1:
            ils = '-.'
        else:
            ils = '--'
            
        (brazil[(ithr,ioz)]['TOTVEGC']/brazil[(ithr,ioz)]['TOTVEGN']).plot(ax=ax11, color=color, ls=ils, alpha=alpha, label="%s, %s" % (ioz,ithr))
        (100*(brazil[(ithr,ioz)]['TOTVEGC']-brazil_0_ref['TOTVEGC'])/brazil_0_ref['TOTVEGC']).plot(ax=ax21, color=color, ls=ils, alpha=alpha, label="%s, %s" % (ioz,ithr))
        (brazil[(ithr,ioz)]['NPP']*sec_in_months*area).plot(ax=ax31, color=color, ls=ils, alpha=alpha, label="%s, %s" % (ioz,ithr))
        (100*(brazil[(ithr,ioz)]['NPP']-brazil_0_ref['NPP'])/brazil_0_ref['NPP']).plot(ax=ax32, color=color, ls=ils, alpha=alpha, label="%s, %s" % (ioz,ithr))
        alpha = alpha-0.05

(brazil_0_ref['TOTVEGC']/brazil_0_ref['TOTVEGN']).plot(ax=ax11, color='black', linewidth=3, label='ref')
(brazil_0_ref['NPP']*sec_in_months*area).plot(ax=ax31, color='black', linewidth=3, label='ref')

ax11.set_xticks(np.arange(1,13))
ax11.set_xticklabels([get_month_name(each, length=3) for each in np.arange(1,13)])

ax21.set_xticks(np.arange(1,13))
ax21.set_xticklabels([get_month_name(each, length=3) for each in np.arange(1,13)])

for ax in fig3.axes:
    ax.set_xticks(np.arange(1,13))
    ax.set_xticklabels("")
ax32.set_xticklabels([get_month_name(each, length=3) for each in np.arange(1,13)])

ax11.legend(ncol=5)
ax21.legend(ncol=5)
ax31.legend(ncol=5)


plt.show(block=False)
