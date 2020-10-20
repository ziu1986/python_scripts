import os, sys, glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all, plot_error_bands, get_month_name, seconds_in_month
from mytools.clm_tools import *

# Clean up
plt.close("all")
b_hydro = False

# Reference data from Almeida et al. (2018)
src_obs = os.environ['DATA'] + '/astra/observations/metadata_OzoneLUNA/2020-09-amazonas_gpp.xlsx'
obs = pd.read_excel(src_obs, header=2)

# CLM simulations
# Reference without hydraulic stress
if b_hydro:
    src_wo_hyd = os.environ['CESM_RUN'] + '/archive/brazil_2000_wohydr_ozone_luna_100/lnd/hist/*'
    wo_hydr_ozone_luna = load_data(src_wo_hyd, var=['NPP','GPP'])
    sum_wo_hydr_ozone_luna = wo_hydr_ozone_luna.apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size*60**2) 

    src_wo_hyd = os.environ['CESM_RUN'] + '/archive/brazil_2000_wohydr/lnd/hist/*'
    wo_hydr_ref = load_data(src_wo_hyd, var=['NPP','GPP'])
    sum_wo_hydr_ref = wo_hydr_ref.apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size*60**2) 

# Senisitvity simulations
src = os.environ['DATA'] + '/preprocessed_data/CLM50_ozone_luna_brazil/'
threshold = (0, 0.2, 0.5, 0.8, 1, 2, 3, 4, 5)
ozone = np.arange(20,110,20)

icolor = {20:'dodgerblue',40:'blue',60:'blueviolet', 80:'deeppink', 100:'red'}
sec_in_months = [seconds_in_month(each,2019) for each in np.arange(1,13)]
sec_in_months = xr.DataArray(data=sec_in_months, coords={'month':np.arange(1,13)},dims='month')
area = 49089.5*1e6 # m**2 taken from original history file

# Mean values
brazil_0_ref = load_data(src + "brazil_0_ref.nc", var=['TOTVEGC', 'TOTVEGN', 'NPP'])
brazil = {}
# Monthly sums
brazil_0_ref_sum = load_data(src + "brazil_npp_0_ref.nc", var=['NPP','GPP'])
brazil_sum = {}
# Annual sum
sums = []

pp = 'GPP'

for ioz in ozone:
    for ithr in threshold:
        try:
            brazil.update({(ithr,ioz):load_data(src + "brazil_%s_%s.nc" % (ioz, ithr), var=['TOTVEGC', 'TOTVEGN', 'NPP'])})
            brazil_sum.update({(ithr,ioz):load_data(src + "brazil_npp_%s_%s.nc" % (ioz, ithr), var=['NPP','GPP'])})
        except ValueError:
            print("(%s, %s) not found." % (ioz,ithr))

# Plot it
fig1, (ax11, ax12) = plt.subplots(1,2, gridspec_kw={'width_ratios': [8, 1]}, figsize=(16,9))
fig1.canvas.set_window_title("%s_C_N_ratio_abs" % pp)
ax12.remove()
fig2, (ax21, ax22) = plt.subplots(1,2, gridspec_kw={'width_ratios': [8, 1]}, figsize=(16,9))
fig2.canvas.set_window_title("%s_C_N_ratio_rel" % pp)
ax22.remove()
fig3, ((ax31, ax33), (ax32, ax34)) = plt.subplots(2,2, gridspec_kw={'width_ratios': [4, 1], 'height_ratios': [1, 1]}, figsize=(12,9))
fig3.canvas.set_window_title("npp_%s" % pp)
ax33.remove()
ax34.remove()

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
            #(brazil[(ithr,ioz)]['NPP']*sec_in_months*area*1e-12).plot(ax=ax31, color=color, ls=ils, alpha=alpha, linewidth=2.5, label="(%s, %s)" % (ioz,ithr))
            (brazil_sum[(ithr,ioz)][pp]*60**2).plot(ax=ax31, color=color, ls=ils, alpha=alpha, linewidth=2.5, label="(%s, %s)" % (ioz,ithr))
            #(100*(brazil[(ithr,ioz)]['NPP']-brazil_0_ref['NPP'])/brazil_0_ref['NPP']).plot(ax=ax32, color=color, ls=ils, alpha=alpha, linewidth=2.5, label="(%s, %s)" % (ioz,ithr))
            (100*(brazil_sum[(ithr,ioz)]-brazil_0_ref_sum)/brazil_0_ref_sum)[pp].plot(ax=ax32, color=color, ls=ils, alpha=alpha, linewidth=2.5, label="(%s, %s)" % (ioz,ithr))
            # Get annual sums of NPP/GPP
            sums.append(brazil_sum[(ithr,ioz)].sum()*60**2*area)
            alpha = alpha-0.1

        except KeyError:
            print("(%s, %s) not found." % (ioz,ithr))

(brazil_0_ref['TOTVEGC']/brazil_0_ref['TOTVEGN']).plot(ax=ax11, color='black', linewidth=3, label='ref')
#(brazil_0_ref['NPP']*sec_in_months*area*1e-12).plot(ax=ax31, color='black', linewidth=3, label='ref')
(brazil_0_ref_sum[pp]*60**2).plot(ax=ax31, color='black', linewidth=3, label='ref')

if b_hydro:
    sum_wo_hydr_ozone_luna[pp].plot(ax=ax31, color='orange', ls='--', label='ozoneluna wohydr')
    sum_wo_hydr_ref[pp].plot(ax=ax31, color='orange', label='ref wohydr')
# Annual sum
sums.append(brazil_0_ref_sum.sum()*60**2*area)

# Reference obs
if pp == 'GPP':
    obs_sel = obs.where((obs.Site=='K34') | (obs.Site=='K67') | (obs.Site=='RJA')).dropna()
    ax31.fill_between(np.arange(1,13), obs_sel.min().values[1:].astype(float), obs_sel.max().values[1:].astype(float), color='grey', alpha=0.5, label='obs')
    
ax11.set_xticks(np.arange(1,13))
ax11.set_xticklabels([get_month_name(each, length=3) for each in np.arange(1,13)])
ax11.set_xlabel("Time (month)")
ax11.set_ylabel("C:N $(gC m^{-2} / gN\,m^{-2})$")

ax21.set_xticks(np.arange(1,13))
ax21.set_xticklabels([get_month_name(each, length=3) for each in np.arange(1,13)])
ax21.set_xlabel("Time (month)")
ax21.set_ylabel("$\Delta$C:N (%)")

for ax in fig3.axes:
    ax.set_xticks(np.arange(1,13))
    ax.set_xticklabels("")
    ax.set_xlabel("")

ax31.set_ylim(0,15)
ax32.set_ylim(-10,0)
ax32.set_xticklabels([get_month_name(each, length=3) for each in np.arange(1,13)])
ax32.set_xlabel("Time (month)")
#ax31.set_ylabel("NPP (TgC)")
ax31.set_ylabel("%s ($gC m^{-2}$)" % pp)
ax32.set_ylabel("$\Delta$%s (%%)" % pp)
ax31.set_title("(a)")
ax32.set_title("(b)")

ax31.text(9,12.35, "$\Sigma NPP_{ref} =  %2.2f\,GgC\,yr^{-1}$\n$\Sigma NPP_{min} = %2.2f\,GgC\,yr^{-1}$" % (xr.concat(sums, dim='exp')[pp].max()*1e-12, 
                                                        xr.concat(sums, dim='exp')[pp].min()*1e-12), size=12)

#ax32.axhline(0, ls=':', color='grey', lw=3)

ax11.legend(title="$([O_3],\Phi_{O_3}^{TH})$", title_fontsize=12, bbox_to_anchor=(1.275, 1),loc='upper right', borderaxespad=0., ncol=2)
ax21.legend(title="$([O_3],\Phi_{O_3}^{TH})$", title_fontsize=12, bbox_to_anchor=(1.275, 1),loc='upper right', borderaxespad=0., ncol=2)
ax31.legend(title="$([O_3],\Phi_{O_3}^{TH})$", title_fontsize=12, bbox_to_anchor=(1.425, 1),loc='upper right', borderaxespad=0., ncol=2)


plt.show(block=False)
