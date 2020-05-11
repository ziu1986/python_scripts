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
case = ('brazil_2000', 'brazil_2000_ozone', 'brazil_2000_ozone_luna_100','brazil_2000_wohydr', 'brazil_2000_wohydr_ozone', 'brazil_2000_wohydr_ozone_luna_100')#, 'brazil_2000_ozone_luna_100_pwu')
#case = ('brazil_2000_wohydr', 'brazil_2000_wohydr_ozone', 'brazil_2000_wohydr_ozone_luna_100')

# Load data
ref_data = pd.read_csv(ref_data_src)

brazil_test = {}

for icase in case:
    brazil_src = run_archive + icase + land_hist
    brazil_test.update({icase:load_data(brazil_src)})


# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("plot_brazil_test_gs_photosyn_absolute")
ax11 = plt.subplot(221)
ax12 = plt.subplot(222)
ax13 = plt.subplot(223)
ax14 = plt.subplot(224)

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("plot_brazil_test_gs_photosyn_relative")
ax21 = plt.subplot(221)
ax22 = plt.subplot(222)
ax23 = plt.subplot(223)
ax24 = plt.subplot(224)

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("plot_brazil_test_gs_photosyn_relative_clim")
ax31 = plt.subplot(221)
ax32 = plt.subplot(222)
ax33 = plt.subplot(223)
ax34 = plt.subplot(224)

fig4 = plt.figure(4, figsize=(16,9))
fig4.canvas.set_window_title("plot_brazil_test_jmax_vcmax_absolute")
ax41 = plt.subplot(221)
ax42 = plt.subplot(222)
ax43 = plt.subplot(223)
ax44 = plt.subplot(224)

fig5 = plt.figure(5, figsize=(16,9))
fig5.canvas.set_window_title("plot_brazil_test_jmax_vcmax_relative")
ax51 = plt.subplot(221)
ax52 = plt.subplot(222)
ax53 = plt.subplot(223)
ax54 = plt.subplot(224)

fig6 = plt.figure(6, figsize=(16,9))
fig6.canvas.set_window_title("plot_brazil_test_jmax_vcmax_relative_clim")
ax61 = plt.subplot(221)
ax62 = plt.subplot(222)
ax63 = plt.subplot(223)
ax64 = plt.subplot(224)

fig7 = plt.figure(7)
fig7.canvas.set_window_title("plot_brazil_test_jmax_vs_vcmax")
ax71 = plt.subplot()

for icase in case:
    # Fig1
    brazil_test[icase]['GSSHA'].plot(ax=ax11, label='%s' % (icase))
    brazil_test[icase]['GSSUN'].plot(ax=ax12, label='%s' % (icase))
    #(1/brazil_test[icase]['RSSHA']).plot(ax=ax11, label='%s' % (icase))
    #(1/brazil_test[icase]['RSSUN']).plot(ax=ax12, label='%s' % (icase))
    brazil_test[icase]['PSNSHA'].plot(ax=ax13, label='%s' % (icase))
    brazil_test[icase]['PSNSUN'].plot(ax=ax14, label='%s' % (icase))

    for ax in fig1.axes:
        ax.set_xlabel("Time (years)")

    # Fig2
    (brazil_test[icase]['GSSHA']-brazil_test[case[0]]['GSSHA']).plot(ax=ax21, label='%s' % (icase))
    (brazil_test[icase]['GSSUN']-brazil_test[case[0]]['GSSUN']).plot(ax=ax22, label='%s' % (icase))
    #(1/brazil_test[icase]['RSSHA']-1/brazil_test[case[0]]['RSSHA']).plot(ax=ax21, label='%s' % (icase))
    #(1/brazil_test[icase]['RSSUN']-1/brazil_test[case[0]]['RSSUN']).plot(ax=ax22, label='%s' % (icase))
    (brazil_test[icase]['PSNSHA']-brazil_test[case[0]]['PSNSHA']).plot(ax=ax23, label='%s' % (icase))
    (brazil_test[icase]['PSNSUN']-brazil_test[case[0]]['PSNSUN']).plot(ax=ax24, label='%s' % (icase))

    for ax in fig2.axes:
        ax.set_xlabel("Time (years)")

    #Fig3
    (brazil_test[icase]['GSSHA']-brazil_test[case[0]]['GSSHA']).groupby('time.dayofyear').mean().plot(ax=ax31, label='%s' % (icase))
    plot_error_bands(ax31, np.arange(1,366), (brazil_test[icase]['GSSHA']-brazil_test[case[0]]['GSSHA']).groupby('time.dayofyear').mean().values.flatten(), error=(brazil_test[icase]['GSSHA']-brazil_test[case[0]]['GSSHA']).groupby('time.dayofyear').std().values.flatten(), color=ax31.lines[-1].get_color(), label='%s' % (icase))

    (brazil_test[icase]['GSSUN']-brazil_test[case[0]]['GSSUN']).groupby('time.dayofyear').mean().plot(ax=ax32, label='%s' % (icase))
    plot_error_bands(ax32, np.arange(1,366), (brazil_test[icase]['GSSUN']-brazil_test[case[0]]['GSSUN']).groupby('time.dayofyear').mean().values.flatten(), error=(brazil_test[icase]['GSSUN']-brazil_test[case[0]]['GSSUN']).groupby('time.dayofyear').std().values.flatten(), color=ax32.lines[-1].get_color(), label='%s' % (icase))

    #(1/brazil_test[icase]['RSSHA']-1/brazil_test[case[0]]['RSSHA']).groupby('time.dayofyear').mean().plot(ax=ax31, label='%s' % (icase))
    #plot_error_bands(ax31, np.arange(1,366), (1/brazil_test[icase]['RSSHA']-1/brazil_test[case[0]]['RSSHA']).groupby('time.dayofyear').mean().values.flatten(), error=(1/brazil_test[icase]['RSSHA']-1/brazil_test[case[0]]['RSSHA']).groupby('time.dayofyear').std().values.flatten(), color=ax31.lines[-1].get_color(), label='%s' % (icase))

    #(1/brazil_test[icase]['RSSUN']-1/brazil_test[case[0]]['RSSUN']).groupby('time.dayofyear').mean().plot(ax=ax32, label='%s' % (icase))
    #plot_error_bands(ax32, np.arange(1,366), (1/brazil_test[icase]['RSSUN']-1/brazil_test[case[0]]['RSSUN']).groupby('time.dayofyear').mean().values.flatten(), error=(1/brazil_test[icase]['RSSUN']-1/brazil_test[case[0]]['RSSUN']).groupby('time.dayofyear').std().values.flatten(), color=ax32.lines[-1].get_color(), label='%s' % (icase))

    (brazil_test[icase]['PSNSHA']-brazil_test[case[0]]['PSNSHA']).groupby('time.dayofyear').mean().plot(ax=ax33, label='%s' % (icase))
    plot_error_bands(ax33, np.arange(1,366), (brazil_test[icase]['PSNSHA']-brazil_test[case[0]]['PSNSHA']).groupby('time.dayofyear').mean().values.flatten(), error=(brazil_test[icase]['PSNSHA']-brazil_test[case[0]]['PSNSHA']).groupby('time.dayofyear').std().values.flatten(), color=ax33.lines[-1].get_color(), label='%s' % (icase))

    (brazil_test[icase]['PSNSUN']-brazil_test[case[0]]['PSNSUN']).groupby('time.dayofyear').mean().plot(ax=ax34, label='%s' % (icase))
    plot_error_bands(ax34, np.arange(1,366), (brazil_test[icase]['PSNSUN']-brazil_test[case[0]]['PSNSUN']).groupby('time.dayofyear').mean().values.flatten(), error=(brazil_test[icase]['PSNSUN']-brazil_test[case[0]]['PSNSUN']).groupby('time.dayofyear').std().values.flatten(), color=ax34.lines[-1].get_color(), label='%s' % (icase))

    for ax in fig3.axes:
        ax.set_xlabel("Time (days of year)")

    ##
    # Fig4
    brazil_test[icase]['JMX25T'].plot(ax=ax41, label='%s' % (icase))
    brazil_test[icase]['Jmx25Z'].plot(ax=ax42, label='%s' % (icase))
    brazil_test[icase]['VCMX25T'].plot(ax=ax43, label='%s' % (icase))
    brazil_test[icase]['Vcmx25Z'].plot(ax=ax44, label='%s' % (icase))

    for ax in fig4.axes:
        ax.set_xlabel("Time (years)")

    # Fig5
    (brazil_test[icase]['JMX25T']-brazil_test[case[0]]['JMX25T']).plot(ax=ax51, label='%s' % (icase))
    (brazil_test[icase]['Jmx25Z']-brazil_test[case[0]]['Jmx25Z']).plot(ax=ax52, label='%s' % (icase))
    (brazil_test[icase]['VCMX25T']-brazil_test[case[0]]['VCMX25T']).plot(ax=ax53, label='%s' % (icase))
    (brazil_test[icase]['Vcmx25Z']-brazil_test[case[0]]['Vcmx25Z']).plot(ax=ax54, label='%s' % (icase))

    for ax in fig5.axes:
        ax.set_xlabel("Time (years)")

    # Fig6
    (brazil_test[icase]['JMX25T']-brazil_test[case[0]]['JMX25T']).groupby('time.dayofyear').mean().plot(ax=ax61, label='%s' % (icase))
    plot_error_bands(ax61, np.arange(1,366), (brazil_test[icase]['JMX25T']-brazil_test[case[0]]['JMX25T']).groupby('time.dayofyear').mean().values.flatten(), error=(brazil_test[icase]['JMX25T']-brazil_test[case[0]]['JMX25T']).groupby('time.dayofyear').std().values.flatten(), color=ax61.lines[-1].get_color(), label='%s' % (icase))

    (brazil_test[icase]['Jmx25Z']-brazil_test[case[0]]['Jmx25Z']).groupby('time.dayofyear').mean().plot(ax=ax62, label='%s' % (icase))
    plot_error_bands(ax62, np.arange(1,366), (brazil_test[icase]['Jmx25Z']-brazil_test[case[0]]['Jmx25Z']).groupby('time.dayofyear').mean().values.flatten(), error=(brazil_test[icase]['Jmx25Z']-brazil_test[case[0]]['Jmx25Z']).groupby('time.dayofyear').std().values.flatten(), color=ax62.lines[-1].get_color(), label='%s' % (icase))

    (brazil_test[icase]['VCMX25T']-brazil_test[case[0]]['VCMX25T']).groupby('time.dayofyear').mean().plot(ax=ax63, label='%s' % (icase))
    plot_error_bands(ax63, np.arange(1,366), (brazil_test[icase]['VCMX25T']-brazil_test[case[0]]['VCMX25T']).groupby('time.dayofyear').mean().values.flatten(), error=(brazil_test[icase]['VCMX25T']-brazil_test[case[0]]['VCMX25T']).groupby('time.dayofyear').std().values.flatten(), color=ax63.lines[-1].get_color(), label='%s' % (icase))

    (brazil_test[icase]['Vcmx25Z']-brazil_test[case[0]]['Vcmx25Z']).groupby('time.dayofyear').mean().plot(ax=ax64, label='%s' % (icase))
    plot_error_bands(ax64, np.arange(1,366), (brazil_test[icase]['Vcmx25Z']-brazil_test[case[0]]['Vcmx25Z']).groupby('time.dayofyear').mean().values.flatten(), error=(brazil_test[icase]['Vcmx25Z']-brazil_test[case[0]]['Vcmx25Z']).groupby('time.dayofyear').std().values.flatten(), color=ax64.lines[-1].get_color(), label='%s' % (icase))

    for ax in fig6.axes:
        ax.set_xlabel("Time (days of year)")
        
    # Fig 7
    ax71.scatter((brazil_test[icase]['Jmx25Z']/brazil_test[case[0]]['Jmx25Z']).values,
                 (brazil_test[icase]['Vcmx25Z']/brazil_test[case[0]]['Vcmx25Z']).values, 
                 marker='x', label='%s' % (icase))

# Referance data
ax71.errorbar(ref_data['Jmax_ratio'], ref_data['Vcmax_ratio'], xerr=ref_data['Jmax_ratio_std'], yerr=ref_data['Vcmax_ratio_std'], ls='None', color='black', label='Ref. data')
# Unweighted fit
from scipy.optimize import curve_fit
def poly_origin(x, m):
    '''
    Line through origin
    '''
    return(m*x)
for icase in case:
    popt, pcov = curve_fit(poly_origin, 
                           (brazil_test[icase]['Jmx25Z']/brazil_test[case[0]]['Jmx25Z']).values.flatten(),
                           (brazil_test[icase]['Vcmx25Z']/brazil_test[case[0]]['Vcmx25Z']).values.flatten(), 
                           [1,])
    print(icase, popt[0], pcov[0])
    #ax71.plot(np.arange(0,2.1,0.1), poly_origin(np.arange(0,2.1,0.1), *popt), label="fit %s" % icase)


ax11.set_ylim(0, 1.5e5)
ax12.set_ylim(0, 1e6)
ax13.set_ylim(0, 2)
ax14.set_ylim(0, 6)

ax11.legend(ncol=2)
ax21.legend(ncol=2)
ax31.legend(ncol=2)
ax41.legend(ncol=2)
ax51.legend(ncol=2)
ax61.legend(ncol=2)
ax71.legend(ncol=2)

for ax in fig2.axes:
    ax.set_ylabel("$\Delta$%s" % ax.get_ylabel())
for ax in fig3.axes:
    ax.set_ylabel("$\Delta$%s" % ax.get_ylabel())


ax41.set_ylabel(ax41.get_ylabel().replace("canopy profile of", "top patch"))
ax42.set_ylabel(ax42.get_ylabel().replace("vcmax", "jmax"))
ax43.set_ylabel(ax43.get_ylabel().replace("canopy profile of", "top patch"))

for ax in fig5.axes:
    ax.set_ylabel("$\Delta$%s" % ax.get_ylabel())
for ax in fig6.axes:
    ax.set_ylabel("$\Delta$%s" % ax.get_ylabel())

ax71.set_xlim(0.3, 1.2)
ax71.set_ylim(0.3, 1.2)
ax71.set_xlabel("$J_{max}^{O_3}/J_{max}^0$")
ax71.set_ylabel("$V_{cmax}^{O_3}/V_{cmax}^0$")
ax71.plot(np.arange(2.1), np.arange(2.1), ls='--', color='grey')


# Show it
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig5.tight_layout()
fig6.tight_layout()
fig7.tight_layout()

plt.show(block=False)

# Pickle the plot to merge with data from metastudy
#import pickle
#with open('brazil_test_vcmax_jmax_ratio.pkl', 'wb') as tgt: # should be 'wb' rather than 'w'
#    pickle.dump(fig7, tgt) 

    
