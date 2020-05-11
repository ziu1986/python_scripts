import os, sys, glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all, plot_error_bands
from clm_tools import *

# Clean up
plt.close("all")

# Source
ref_data_src = os.environ['PY_SCRIPTS'] + '/plant_model/test.cvs'
run_archive = os.environ['CESM_RUN'] + '/archive/'
land_hist = '/lnd/hist/*.clm2.h0.*.nc'
case = ('brazil_2000_ozone_luna_100', 'brazil_2000_ozone_luna_100_thresh_', 'brazil_2000_ozone_luna_')
threshold = (0.8, 0.4, 0.6, 0.7, 0.85, 0.9, 1)
ozone = (100, 0, 40, 60, 80)

# Load data
ref_data = pd.read_csv(ref_data_src)

brazil_test = {}
brazil_test_ozone = {}
brazil_ref_ozone = {}
# Load reference simulation
brazil_src = run_archive + case[0] + land_hist
brazil_test.update({threshold[0]:load_data(brazil_src)})
brazil_test_ozone.update({ozone[0]:load_data(brazil_src)})
brazil_test_ozone.update({ozone[1]:load_data(brazil_src.replace('_ozone_luna_100', ''))})
brazil_ref_ozone.update({ozone[0]:load_data(brazil_src.replace('_luna_100', ''))})
brazil_ref_ozone.update({ozone[1]:load_data(brazil_src.replace('_ozone_luna_100', ''))})

# Load sensitivity tests
for ithresh in threshold[1:]:
    brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
    brazil_test.update({ithresh:load_data(brazil_src)})
for iozone in ozone[2:]:
    brazil_src = run_archive + case[2] + "%s" % iozone + land_hist
    brazil_test_ozone.update({iozone:load_data(brazil_src)})
    brazil_ref_ozone.update({iozone:load_data(brazil_src.replace('luna_',''))})


# Plot it
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("OzoneLunaMod_threshold_sensitivity")
ax11 = plt.subplot(221)
ax12 = plt.subplot(222)
ax13 = plt.subplot(223)
ax14 = plt.subplot(224)

for ithresh in threshold:
    probe = (brazil_test[ithresh]-brazil_test[threshold[0]])#/brazil_test[threshold[0]]*100

    ax11.errorbar(ithresh, probe['GSSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSHA'], color='blue', marker='o', fillstyle='none', label="shade")
    ax11.errorbar(ithresh, probe['GSSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSUN'], color='black', marker='o', label='sun')

    ax12.errorbar(ithresh, probe['PSNSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSHA'], color='blue', marker='o', fillstyle='none', label='shade')
    ax12.errorbar(ithresh, probe['PSNSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSUN'], color='black', marker='o', label='sun')

    ax13.errorbar(ithresh, probe['Jmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Jmx25Z'], color='blue', marker='o', fillstyle='none', label='avg')
    ax13.errorbar(ithresh, probe['JMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['JMX25T'], color='black', marker='o', label='top')

    ax14.errorbar(ithresh, probe['Vcmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Vcmx25Z'], color='blue', marker='o', fillstyle='none', label='avg')
    ax14.errorbar(ithresh, probe['VCMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['VCMX25T'], color='black', marker='o', label='top')

for ax in fig1.axes:
    ax.set_yscale('symlog', linthreshy=1e-10)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2])
    ax.axhline(1e-8, color='red', ls='--')
    ax.axhline(-1e-8, color='red', ls='--')
    ax.axhline(0, color='grey', ls=':')
    ax.axhspan(-1e-8, 1e-8, edgecolor='black', facecolor='None', hatch="//")


#ax11.set_ylabel("$\Delta G_{sto}$ (%)")
#ax12.set_ylabel("$\Delta A_{n}$ (%)")
#ax13.set_ylabel("$\Delta J_{max}$ (%)")
#ax14.set_ylabel("$\Delta V_{cmax}$ (%)")

ax11.set_ylabel("$\Delta G_{sto}$ ($\mu mol H_20 m^{-2}s^{-1}$)")
ax12.set_ylabel("$\Delta A_{n}$ ($\mu molCO_2 m^{-2}s^{-1}$)")
ax13.set_ylabel("$\Delta J_{max}$ ($\mu mol m^{-2}s^{-1}$)")
ax14.set_ylabel("$\Delta V_{cmax}$ ($\mu mol m^{-2}s^{-1}$)")


ax13.set_xlabel("$O_3^{threshold}$", x=1)

fig2 = plt.figure(2,figsize=(16,9))
fig2.canvas.set_window_title("ozone_sensitivity")
ax21 = plt.subplot(221)
ax22 = plt.subplot(222)
ax23 = plt.subplot(223)
ax24 = plt.subplot(224)

for ithresh in ozone:
    probe = (brazil_test_ozone[ithresh]-brazil_test_ozone[ozone[0]])#/brazil_test_ozone[ozone[0]]*100
    
    ax21.errorbar(ithresh, probe['GSSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSHA'], color='blue', marker='o', fillstyle='none', label="shade - OzoneLunaMod")
    ax21.errorbar(ithresh, probe['GSSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSUN'], color='blue', marker='o', label='sun - OzoneLunaMod')

    ax22.errorbar(ithresh, probe['PSNSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSHA'], color='blue', marker='o', fillstyle='none', label='shade - OzoneLunaMod')
    ax22.errorbar(ithresh, probe['PSNSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSUN'], color='blue', marker='o', label='sun - OzoneLunaMod')

    ax23.errorbar(ithresh, probe['Jmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Jmx25Z'], color='blue', marker='o', fillstyle='none', label='avg - OzoneLunaMod')
    ax23.errorbar(ithresh, probe['JMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['JMX25T'], color='blue', marker='o', label='top - OzoneLunaMod')

    ax24.errorbar(ithresh, probe['Vcmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Vcmx25Z'], color='blue', marker='o', fillstyle='none', label='avg - OzoneLunaMod')
    ax24.errorbar(ithresh, probe['VCMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['VCMX25T'], color='blue', marker='o', label='top - OzoneLunaMod')

    probe = (brazil_ref_ozone[ithresh]-brazil_ref_ozone[ozone[0]])#/brazil_ref_ozone[ozone[0]]*100
    
    ax21.errorbar(ithresh, probe['GSSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSHA'], color='black', marker='s', fillstyle='none', label="shade - OzoneMod")
    ax21.errorbar(ithresh, probe['GSSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSUN'], color='black', marker='s', label='sun - OzoneMod')

    ax22.errorbar(ithresh, probe['PSNSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSHA'], color='black', marker='s', fillstyle='none', label='shade - OzoneMod')
    ax22.errorbar(ithresh, probe['PSNSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSUN'], color='black', marker='s', label='sun - OzoneMod')

    ax23.errorbar(ithresh, probe['Jmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Jmx25Z'], color='black', marker='s', fillstyle='none', label='avg - OzoneMod')
    ax23.errorbar(ithresh, probe['JMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['JMX25T'], color='black', marker='s', label='top - OzoneMod')

    ax24.errorbar(ithresh, probe['Vcmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Vcmx25Z'], color='black', marker='s', fillstyle='none', label='avg - OzoneMod')
    ax24.errorbar(ithresh, probe['VCMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['VCMX25T'], color='black', marker='s', label='top - OzoneMod')


for ax in fig2.axes:
    ax.set_yscale('symlog', linthreshy=1e-10)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:4], labels[:4], ncol=2)
    ax.axhline(1e-8, color='red', ls='--')
    ax.axhline(-1e-8, color='red', ls='--')
    ax.axhline(0, color='grey', ls=':')
    ax.axhspan(-1e-8, 1e-8, edgecolor='black', facecolor='None', hatch="//")

ax21.set_ylabel("$\Delta G_{sto}$ ($\mu mol H_20 m^{-2}s^{-1}$)")
ax22.set_ylabel("$\Delta A_{n}$ ($\mu molCO_2 m^{-2}s^{-1}$)")
ax23.set_ylabel("$\Delta J_{max}$ ($\mu mol m^{-2}s^{-1}$)")
ax24.set_ylabel("$\Delta V_{cmax}$ ($\mu mol m^{-2}s^{-1}$)")

#ax21.set_ylabel("$\Delta G_{sto}$ (%)")
#ax22.set_ylabel("$\Delta A_{n}$ (%)")
#ax23.set_ylabel("$\Delta J_{max}$ (%)")
#ax24.set_ylabel("$\Delta V_{cmax}$ (%)")

ax23.set_xlabel("$[O_3]$ (ppb)", x=1)

"""
fig3 = plt.figure(3,figsize=(16,9))
fig3.canvas.set_window_title("OzoneMod_ozone_sensitivity")
ax31 = plt.subplot(221)
ax32 = plt.subplot(222)
ax33 = plt.subplot(223)
ax34 = plt.subplot(224)

for ithresh in ozone:
    probe = (brazil_ref_ozone[ithresh]-brazil_ref_ozone[ozone[0]])#/brazil_ref_ozone[ozone[0]]*100
    
    ax31.errorbar(ithresh, probe['GSSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSHA'], color='blue', marker='o', fillstyle='none', label="shade")
    ax31.errorbar(ithresh, probe['GSSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSUN'], color='black', marker='o', label='sun')

    ax32.errorbar(ithresh, probe['PSNSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSHA'], color='blue', marker='o', fillstyle='none', label='shade')
    ax32.errorbar(ithresh, probe['PSNSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSUN'], color='black', marker='o', label='sun')

    ax33.errorbar(ithresh, probe['Jmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Jmx25Z'], color='blue', marker='o', fillstyle='none', label='avg')
    ax33.errorbar(ithresh, probe['JMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['JMX25T'], color='black', marker='o', label='top')

    ax34.errorbar(ithresh, probe['Vcmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Vcmx25Z'], color='blue', marker='o', fillstyle='none', label='avg')
    ax34.errorbar(ithresh, probe['VCMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['VCMX25T'], color='black', marker='o', label='top')

for ax in fig3.axes:
    ax.set_yscale('symlog', linthreshy=1e-10)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2])
    ax.axhline(1e-8, color='red', ls='--')
    ax.axhline(-1e-8, color='red', ls='--')
    ax.axhline(0, color='grey', ls=':')
    ax.axhspan(-1e-8, 1e-8, edgecolor='black', facecolor='None', hatch="//")

ax31.set_ylabel("$\Delta G_{sto}$ ($\mu mol H_20 m^{-2}s^{-1}$)")
ax32.set_ylabel("$\Delta A_{n}$ ($\mu molCO_2 m^{-2}s^{-1}$)")
ax33.set_ylabel("$\Delta J_{max}$ ($\mu mol m^{-2}s^{-1}$)")
ax34.set_ylabel("$\Delta V_{cmax}$ ($\mu mol m^{-2}s^{-1}$)")

#ax31.set_ylabel("$\Delta G_{sto}$ (%)")
#ax32.set_ylabel("$\Delta A_{n}$ (%)")
#ax33.set_ylabel("$\Delta J_{max}$ (%)")
#ax34.set_ylabel("$\Delta V_{cmax}$ (%)")

ax33.set_xlabel("$[O_3]$ (ppb)", x=1)
"""
# Show it
plt.show(block=False)
