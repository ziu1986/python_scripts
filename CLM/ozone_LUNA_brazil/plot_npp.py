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
ozone = np.arange(20,110,20)

brazil_0_ref = load_data("brazil_0_ref.nc")
brazil = {}

for ioz in ozone:
    for ithr in threshold:
        try:
            brazil.update({(ithr,ioz):load_data("brazil_%s_%s.nc" % (ioz, ithr))})
        except ValueError:
            print("(%s, %s) not found." % (ioz,ithr))
# Plot it
        

fig1, (ax11, ax12) = plt.subplots(1,2, gridspec_kw={'width_ratios': [8, 1]}, figsize=(16,9))
fig1.canvas.set_window_title("npp_C_N_ratio_abs")
ax12.remove()
fig2, (ax21, ax22) = plt.subplots(1,2, gridspec_kw={'width_ratios': [8, 1]}, figsize=(16,9))
fig2.canvas.set_window_title("npp_C_N_ratio_rel")
ax22.remove()
fig3, ((ax31, ax33), (ax32, ax34)) = plt.subplots(2,2, gridspec_kw={'width_ratios': [8, 1], 'height_ratios': [1, 1]}, figsize=(16,9))
fig3.canvas.set_window_title("npp_npp")
ax33.remove()
ax34.remove()

icolor = {20:'dodgerblue',40:'blue',60:'blueviolet', 80:'deeppink', 100:'red'}
sec_in_months = [seconds_in_month(each,2019) for each in np.arange(1,13)]
sec_in_months = xr.DataArray(data=sec_in_months, coords={'month':np.arange(1,13)},dims='month')
area = 49089.5*1e6 # m**2

for ioz in ozone:
    color = icolor[ioz]
    alpha = 1
    for ithr in threshold:
        #print(brazil[(ithr,ioz)]['TOTVEGC'])
        if ithr == 0:
            ils = ':'
        elif ithr <= 1:
            ils = '-.'
        else:
            ils = '--'
        try:
            (brazil[(ithr,ioz)]['TOTVEGC']/brazil[(ithr,ioz)]['TOTVEGN']).plot(ax=ax11, color=color, linewidth=2.5, ls=ils, alpha=alpha, label="(%s, %s)" % (ioz,ithr))
            (100*(brazil[(ithr,ioz)]['TOTVEGC']/brazil[(ithr,ioz)]['TOTVEGN']-brazil_0_ref['TOTVEGC']/brazil_0_ref['TOTVEGN'])/brazil_0_ref['TOTVEGC']*brazil_0_ref['TOTVEGN']).plot(ax=ax21, color=color, ls=ils, alpha=alpha, linewidth=2.5, label="%s, %s" % (ioz,ithr))
            (brazil[(ithr,ioz)]['NPP']*sec_in_months*area*1e-12).plot(ax=ax31, color=color, ls=ils, alpha=alpha, linewidth=2.5, label="(%s, %s)" % (ioz,ithr))
            (100*(brazil[(ithr,ioz)]['NPP']-brazil_0_ref['NPP'])/brazil_0_ref['NPP']).plot(ax=ax32, color=color, ls=ils, alpha=alpha, linewidth=2.5, label="(%s, %s)" % (ioz,ithr))
        
            alpha = alpha-0.1
        except KeyError:
            print("(%s, %s) not found." % (ioz,ithr))

(brazil_0_ref['TOTVEGC']/brazil_0_ref['TOTVEGN']).plot(ax=ax11, color='black', linewidth=3, label='ref')
(brazil_0_ref['NPP']*sec_in_months*area*1e-12).plot(ax=ax31, color='black', linewidth=3, label='ref')

ax11.set_xticks(np.arange(1,13))
ax11.set_xticklabels([get_month_name(each, length=3) for each in np.arange(1,13)])
ax11.set_xlabel("Time (month)")
ax11.set_ylabel("C:N $(gC m^{-2} / gN m^{-2})$")

ax21.set_xticks(np.arange(1,13))
ax21.set_xticklabels([get_month_name(each, length=3) for each in np.arange(1,13)])
ax21.set_xlabel("Time (month)")
ax21.set_ylabel("$\Delta$C:N (%)")

for ax in fig3.axes:
    ax.set_xticks(np.arange(1,13))
    ax.set_xticklabels("")
    ax.set_xlabel("")
ax32.set_xticklabels([get_month_name(each, length=3) for each in np.arange(1,13)])
ax32.set_xlabel("Time (month)")
ax31.set_ylabel("NPP (TgC)")
ax32.set_ylabel("$\Delta$NPP (%)")

ax11.legend(title="$([O_3],\Phi_{O_3}^{thresh})$", title_fontsize=12, bbox_to_anchor=(1.275, 1),loc='upper right', borderaxespad=0., ncol=2)
ax21.legend(title="$([O_3],\Phi_{O_3}^{thresh})$", title_fontsize=12, bbox_to_anchor=(1.275, 1),loc='upper right', borderaxespad=0., ncol=2)
ax31.legend(title="$([O_3],\Phi_{O_3}^{thresh})$", title_fontsize=12, bbox_to_anchor=(1.275, 1),loc='upper right', borderaxespad=0., ncol=2)


plt.show(block=False)
