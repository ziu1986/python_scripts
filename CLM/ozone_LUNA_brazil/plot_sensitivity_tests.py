import os, sys, glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all, plot_error_bands
from mytools.clm_tools import *

def plot_data(fig, data, iter_data, **karg):
    # option: ozone, threshold
    mode = karg.pop('mode', 'threshold')
    # option: abs, rel
    scale = karg.pop('scale', 'abs')
    comp_idx = karg.pop('compare_index', 0)
    b_3d = karg.pop('3d', False)
    # Plot adjustments
    label = karg.pop('label', '')
    color = karg.pop('color', 'blue')
    marker = karg.pop('marker', 'o')

    # Define subplots
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)

    # Iterate through data
    for ithresh in iter_data:
        if scale == 'abs':
            probe = (data[ithresh]-data[iter_data[comp_idx]])
        else:
            probe = (data[ithresh]-data[iter_data[comp_idx]])/data[iter_data[comp_idx]]*100

        ax1.errorbar(ithresh, probe['GSSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSHA'], color=color, marker=marker, fillstyle='none', label="shade%s" % label)
        ax1.errorbar(ithresh, probe['GSSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['GSSUN'], color=color, marker=marker, label='sun%s' % label)

        ax2.errorbar(ithresh, probe['PSNSHA'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSHA'], color=color, marker=marker, fillstyle='none', label='shade%s' % label)
        ax2.errorbar(ithresh, probe['PSNSUN'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['PSNSUN'], color=color, marker=marker, label='sun%s' % label)

        ax3.errorbar(ithresh, probe['Jmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Jmx25Z'], color=color, marker=marker, fillstyle='none', label='avg%s' % label)
        ax3.errorbar(ithresh, probe['JMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['JMX25T'], color=color, marker=marker, label='top%s' % label)

        ax4.errorbar(ithresh, probe['Vcmx25Z'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['Vcmx25Z'], color=color, marker=marker, fillstyle='none', label='avg%s' % label)
        ax4.errorbar(ithresh, probe['VCMX25T'].mean(), yerr=probe.apply(lambda x: x.std()/np.sqrt(x.size))['VCMX25T'], color=color, marker=marker, label='top%s' % label)

    if scale == 'abs':
        for ax in fig.axes:
            ax.set_yscale('symlog', linthreshy=1e-16)
            ax.axhline(1e-15, color='red', ls='--')
            ax.axhline(-1e-15, color='red', ls='--')
            ax.axhline(0, color='grey', ls=':')
            #ax.axhspan(-1e-15, 1e-15, edgecolor='black', facecolor='None', hatch="//")

        ax1.set_ylabel("$\Delta_{0.8, 100} G_{sto}$ ($\mu mol H_2O m^{-2}s^{-1}$)")
        ax2.set_ylabel("$\Delta_{0.8, 100} A_{n}$ ($\mu molCO_2 m^{-2}s^{-1}$)")
        ax3.set_ylabel("$\Delta_{0.8, 100} J_{max}$ ($\mu mol m^{-2}s^{-1}$)")
        ax4.set_ylabel("$\Delta_{0.8, 100} V_{cmax}$ ($\mu mol m^{-2}s^{-1}$)")
    else:
        for ax in fig.axes:
            ax.set_yscale('symlog', linthreshy=1e-9)

        ax1.set_ylabel("$\Delta_{0.8, 100} G_{sto}$ (%)")
        ax2.set_ylabel("$\Delta_{0.8, 100} A_{n}$ (%)")
        ax3.set_ylabel("$\Delta_{0.8, 100} J_{max}$ (%)")
        ax4.set_ylabel("$\Delta_{0.8, 100} V_{cmax}$ (%)")

    if mode == 'threshold':
        ax3.set_xlabel("$O_3^{TH}$ (nmol $m^{-2} s^{-1}$)", x=1)
    else:
        ax3.set_xlabel("$[O_3]$ (ppb)", x=1)

# 3D plot
def plot_3d(**karg):
    from mpl_toolkits.mplot3d import Axes3D # Register 3d projection
    import matplotlib as mpl
    variable = karg.pop("variable", "GSSUN")
    b_rel = karg.pop("relative", False)
    rotation = karg.pop('rot', (30, 157))

    fig5 = plt.figure(5)
    if b_rel:
        fig5.canvas.set_window_title("ozone_sensitivity_rel_3d_%s" % (variable))
    else:
        fig5.canvas.set_window_title("ozone_sensitivity_3d_%s" % (variable))
    ax51 = plt.subplot(projection='3d')
    cmap = plt.cm.get_cmap('BuPu_r')#('viridis')
    
    z_data = []
    x_data = []
    y_data = []
   
    if variable == 'C:N':
        # Define colormap norm
        vmin = min((np.log(brazil_test[0].mean()['TOTVEGC']/
                           brazil_test[0].mean()['TOTVEGN']),
                    np.log(brazil_test_20[5].mean()['TOTVEGC']/
                           brazil_test_20[5].mean()['TOTVEGN'])))
        vmax = max((np.log(brazil_test[0].mean()['TOTVEGC']/
                           brazil_test[0].mean()['TOTVEGN']),
                    np.log(brazil_test_20[5].mean()['TOTVEGC']/
                           brazil_test_20[5].mean()['TOTVEGN'])))
        
        #norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        for ithresh in threshold:
            x_data.append(ithresh)
            y_data.append(100)
            z_data.append((brazil_test[ithresh]['TOTVEGC']/brazil_test[ithresh]['TOTVEGN']).mean().values)
            
        for ithresh in threshold_2:
            x_data.append(ithresh)
            y_data.append(0)
            z_data.append((brazil_ref_ozone[0]['TOTVEGC']/brazil_ref_ozone[0]['TOTVEGN']).mean().values)

            x_data.append(ithresh)
            y_data.append(20)
            z_data.append((brazil_test_20[ithresh]['TOTVEGC']/brazil_test_20[ithresh]['TOTVEGN']).mean().values)
        
            x_data.append(ithresh)
            y_data.append(40)
            z_data.append((brazil_test_40[ithresh]['TOTVEGC']/brazil_test_40[ithresh]['TOTVEGN']).mean().values)
        
            x_data.append(ithresh)
            y_data.append(60)
            z_data.append((brazil_test_60[ithresh]['TOTVEGC']/brazil_test_60[ithresh]['TOTVEGN']).mean().values)
   
            x_data.append(ithresh)
            y_data.append(80)
            z_data.append((brazil_test_80[ithresh]['TOTVEGC']/brazil_test_80[ithresh]['TOTVEGN']).mean().values)

    elif variable.find("g:A")>=0:
        if variable.find("sun")>=0:
            light = "SUN"
        else:
            light = "SHA"
        print(['GS%s' % light])
        # Define colormap norm
        vmin = min((np.log(brazil_test[0].mean()['GS%s' % light]/
                           brazil_test[0].mean()['PSN%s' % light]),
                    np.log(brazil_test_20[5].mean()['GS%s' % light]/
                           brazil_test_20[5].mean()['PSN%s' % light])))
        vmax = max((np.log(brazil_test[0].mean()['GS%s' % light]/
                           brazil_test[0].mean()['PSN%s' % light]),
                    np.log(brazil_test_20[5].mean()['GS%s' % light]/
                           brazil_test_20[5].mean()['PSN%s' % light])))
        
        #norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        for ithresh in threshold:
            x_data.append(ithresh)
            y_data.append(100)
            z_data.append((brazil_test[ithresh]['GS%s' % light]/brazil_test[ithresh]['PSN%s' % light]).mean().values)
            
        for ithresh in threshold_2:
            x_data.append(ithresh)
            y_data.append(0)
            z_data.append((brazil_ref_ozone[0]['GS%s' % light]/brazil_ref_ozone[0]['PSN%s' % light]).mean().values)

            x_data.append(ithresh)
            y_data.append(20)
            z_data.append((brazil_test_20[ithresh]['GS%s' % light]/brazil_test_20[ithresh]['PSN%s' % light]).mean().values)
        
            x_data.append(ithresh)
            y_data.append(40)
            z_data.append((brazil_test_40[ithresh]['GS%s' % light]/brazil_test_40[ithresh]['PSN%s' % light]).mean().values)
        
            x_data.append(ithresh)
            y_data.append(60)
            z_data.append((brazil_test_60[ithresh]['GS%s' % light]/brazil_test_60[ithresh]['PSN%s' % light]).mean().values)
   
            x_data.append(ithresh)
            y_data.append(80)
            z_data.append((brazil_test_80[ithresh]['GS%s' % light]/brazil_test_80[ithresh]['PSN%s' % light]).mean().values)
    
    else:
        vmin = min((np.log(brazil_test[0].mean()[variable]),np.log(brazil_test_20[5].mean()[variable])))
        vmax = max((np.log(brazil_test[0].mean()[variable]),np.log(brazil_test_20[5].mean()[variable])))
        #norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        for ithresh in threshold:
            x_data.append(ithresh)
            y_data.append(100)
            z_data.append((brazil_test[ithresh][variable].mean().values))
            
        for ithresh in threshold_2:

            x_data.append(ithresh)
            y_data.append(0)
            z_data.append((brazil_ref_ozone[0][variable].mean().values))

            x_data.append(ithresh)
            y_data.append(20)
            z_data.append((brazil_test_20[ithresh][variable].mean().values))
            
            x_data.append(ithresh)
            y_data.append(40)
            z_data.append((brazil_test_40[ithresh][variable].mean().values))
           
            x_data.append(ithresh)
            y_data.append(60)
            z_data.append((brazil_test_60[ithresh][variable].mean().values))   
            
            x_data.append(ithresh)
            y_data.append(80)
            z_data.append((brazil_test_80[ithresh][variable].mean().values))
               
    
    print("Min %s\n Max %s" %(vmin, vmax))
    # Plot it
    x_data_plot = np.array(x_data)[16:].reshape(len(threshold_2),5)
    y_data_plot = np.array(y_data)[16:].reshape(len(threshold_2),5)
    z_data_plot = np.array(z_data)[16:].reshape(len(threshold_2),5)
    ax51.plot_wireframe(x_data_plot, y_data_plot, z_data_plot, color='grey')
    ax51.plot_surface(x_data_plot, y_data_plot, z_data_plot, cmap='BuPu')
    #z = cmap(norm(np.log(z_data)))
    ax51.scatter(x_data, y_data, z_data, c='black', s=50)
    ax51.ticklabel_format(axis='z', style='sci')
   
    
    #ax51.set_zscale('log')
    # Rotation
    ax51.view_init(*rotation)
    ax51.set_xlabel("$O_3^{TH}$ (nmol $m^{-2} s^{-1}$)")
    ax51.set_ylabel("$[O_3]$ (ppb)")
    if variable.find('GS')>=0:
        ax51.set_zlabel("$G_{sto}$ ($\mu mol H_2O m^{-2}s^{-1}$)")
    elif variable.find('PSN')>=0:
        ax51.set_zlabel("$A_{n}$ ($\mu molCO_2 m^{-2}s^{-1}$)")
    elif variable.find("Jmx")>=0:
        ax51.set_zlabel("$J_{max}$ ($\mu mol m^{-2}s^{-1}$)")
    elif variable.find("Vcmx")>=0:
        ax51.set_zlabel("$V_{cmax}$ ($\mu mol m^{-2}s^{-1}$)")
    elif variable.find("TOTVEGC")>=0:
        ax51.set_zlabel("$C_{tot}^{veg}$ ($g m^{-2}$)")
    elif variable.find("TOTVEGN")>=0:
        ax51.set_zlabel("$N_{tot}^{veg}$ ($g m^{-2}$)")
    elif variable.find("NPP")>=0:
        ax51.set_zlabel("$NPP$ ($g m^{-2}s^{-1}$)")
    elif variable.find("GPP")>=0:
        ax51.set_zlabel("$GPP$ ($g m^{-2}s^{-1}$)")
    elif variable.find("g:A")>=0:
        ax51.set_zlabel("$G_{sto}/A_{n}$ ($\mu mol H_2O/CO_2 m^{-2}s^{-1}$)")
        print(z_data)
    else:
        ax51.set_zlabel('C:N')

def load_date():
    # Load reference simulation
    brazil_src = run_archive + case[0] + land_hist
    brazil_test.update({threshold[0]:load_data(brazil_src)})
    brazil_test_40.update({threshold[0]:load_data(brazil_src.replace('100', '40'))})
    brazil_test_20.update({threshold[0]:load_data(brazil_src.replace('100', '20'))})
    brazil_test_60.update({threshold[0]:load_data(brazil_src.replace('100', '60'))})
    brazil_test_80.update({threshold[0]:load_data(brazil_src.replace('100', '80'))})
    brazil_test_ozone.update({ozone[0]:load_data(brazil_src)})
    brazil_test_ozone.update({ozone[1]:load_data(brazil_src.replace('_ozone_luna_100', ''))})
    brazil_ref_ozone.update({ozone[0]:load_data(brazil_src.replace('_luna_100', ''))})
    brazil_ref_ozone.update({ozone[1]:load_data(brazil_src.replace('_ozone_luna_100', ''))})

    # Load sensitivity tests
    for ithresh in threshold[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test.update({ithresh:load_data(brazil_src)})
    for ithresh in threshold_2[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test_40.update({ithresh:load_data(brazil_src.replace('100', '40'))})
    for ithresh in threshold_2[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test_20.update({ithresh:load_data(brazil_src.replace('100', '20'))})
    for ithresh in threshold_2[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test_60.update({ithresh:load_data(brazil_src.replace('100', '60'))})
    for ithresh in threshold_2[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test_80.update({ithresh:load_data(brazil_src.replace('100', '80'))})


    for iozone in ozone[2:]:
        brazil_src = run_archive + case[2] + "%s" % iozone + land_hist
        brazil_test_ozone.update({iozone:load_data(brazil_src)})
        brazil_ref_ozone.update({iozone:load_data(brazil_src.replace('luna_',''))})
    for iozone in ozone_deep:
        brazil_src = run_archive + case[2] + "%s" % iozone + land_hist
        brazil_test_ozone_deep.update({iozone:load_data(brazil_src)})

    # Merge fine resolution scan ozone deep
    brazil_test_ozone.update(brazil_test_ozone_deep)

def save_data(**karg):
    directory = karg.pop('dir', os.environ['DATA'] + "/preprocessed_data/CLM50_ozone_luna_brazil")
    for ithr in threshold:
        brazil_test[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_100_%s.nc" % (directory, ithr))
        brazil_test[ithr][['NPP','GPP', 'TLAI']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_100_%s.nc" % (directory, ithr))
    for ithr in threshold_2:
        brazil_test_40[ithr][['NPP','GPP', 'TLAI']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_40_%s.nc" % (directory, ithr))
        brazil_test_20[ithr][['NPP','GPP', 'TLAI']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_20_%s.nc" % (directory, ithr))
        brazil_test_60[ithr][['NPP','GPP', 'TLAI']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_60_%s.nc" % (directory, ithr))
        brazil_test_80[ithr][['NPP','GPP', 'TLAI']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_80_%s.nc" % (directory, ithr))
        brazil_test_40[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_40_%s.nc" % (directory, ithr))
        brazil_test_20[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_20_%s.nc" % (directory, ithr))
        brazil_test_60[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_60_%s.nc" % (directory, ithr))
        brazil_test_80[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_80_%s.nc" % (directory, ithr))
    for ioz in ozone:
        brazil_ref_ozone[ioz].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_%s_ref.nc" % (directory,ioz))
        brazil_ref_ozone[ioz].groupby('time.month').mean().to_netcdf("%s/brazil_%s_ref.nc" % (directory,ioz))

# Clean up
plt.close("all")

# Source
ref_data_src = os.environ['PY_SCRIPTS'] + '/plant_model/test.cvs'
try:
    run_archive = os.environ['CESM_RUN'] + '/archive/'
except KeyError:
    run_archive = os.environ['DATA'] + '/astra_data/clm_results/'
    
land_hist = '/lnd/hist/*.clm2.h0.*.nc'
case = ('brazil_2000_ozone_luna_100', 'brazil_2000_ozone_luna_100_thresh_', 'brazil_2000_ozone_luna_')
threshold = (0.8, 0.0, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.7, 0.85, 0.9, 1, 2, 3, 4, 5)
threshold_2 = (0.8, 0, 0.2, 0.5, 1, 2, 3, 4, 5)
ozone = (100, 0, 40, 60, 80)
ozone_deep = np.arange(42, 60, 2)

# Load data

ref_data = pd.read_csv(ref_data_src)

brazil_test = {}
brazil_test_40 = {}
brazil_test_20 = {}
brazil_test_60 = {}
brazil_test_80 = {}
brazil_test_ozone = {}
brazil_test_ozone_deep = {}
brazil_ref_ozone = {}

load_date()


# Plot it
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("ozone_threshold_sensitivity")
plot_data(fig1, brazil_test, threshold, label='_100', compare_index=1)
plot_data(fig1, brazil_test_40, threshold_2, label='_40', color='black', marker='s', compare_index=1)
plot_data(fig1, brazil_test_20, threshold_2, label='_20', color='black', marker='d', compare_index=1)
plot_data(fig1, brazil_test_60, threshold_2, label='_60', color='black', marker='^', compare_index=1)
plot_data(fig1, brazil_test_80, threshold_2, label='_80', color='black', marker='v', compare_index=1)

for ax in fig1.axes:
    handles, labels = ax.get_legend_handles_labels()
    labels_new = np.unique(labels)
    handles_new = [handles[i] for i in np.unique(labels, return_index=True)[1]]
    ax.legend(handles_new, labels_new, ncol=2, loc='center right')


fig2 = plt.figure(2,figsize=(16,9))
fig2.canvas.set_window_title("ozone_sensitivity")
plot_data(fig2, brazil_test_ozone, np.append(ozone, ozone_deep), label=' - OzoneLunaMod', mode='ozone')
plot_data(fig2, brazil_ref_ozone, ozone, label=' - OzoneMod', mode='ozone', color='black', marker='s')
#plot_data(fig2, brazil_test_ozone_deep, ozone_deep, label=' - OzoneLunaMod', mode='ozone')

for ax in fig2.axes:
    handles, labels = ax.get_legend_handles_labels()
    labels_new = np.unique(labels)
    handles_new = [handles[i] for i in np.unique(labels, return_index=True)[1]]
    ax.legend(handles_new, labels_new, ncol=2, loc='center right')


fig3 = plt.figure(3,figsize=(16,9))
fig3.canvas.set_window_title("ozone_threshold_sensitivity_rel")
plot_data(fig3, brazil_test, threshold, label='_100', scale='rel', compare_index=1)
plot_data(fig3, brazil_test_40, threshold_2, label='_40', color='black', marker='s', scale='rel', compare_index=1)
plot_data(fig3, brazil_test_20, threshold_2, label='_20', color='black', marker='d', scale='rel', compare_index=1)
plot_data(fig3, brazil_test_60, threshold_2, label='_60', color='black', marker='^', scale='rel', compare_index=1)
plot_data(fig3, brazil_test_80, threshold_2, label='_80', color='black', marker='v', scale='rel', compare_index=1)

for ax in fig3.axes:
    handles, labels = ax.get_legend_handles_labels()
    labels_new = np.unique(labels)
    handles_new = [handles[i] for i in np.unique(labels, return_index=True)[1]]
    ax.legend(handles_new, labels_new, ncol=2, loc='center right')


fig4 = plt.figure(4,figsize=(16,9))
fig4.canvas.set_window_title("ozone_sensitivity_rel")
plot_data(fig4, brazil_test_ozone, np.append(ozone, ozone_deep), label=' - OzoneLunaMod', mode='ozone', scale='rel', compare_index=1)
plot_data(fig4, brazil_ref_ozone, ozone, label=' - OzoneMod', mode='ozone', color='black', marker='s', scale='rel', compare_index=1)
#plot_data(fig4, brazil_test_ozone_deep, ozone_deep, label=' - OzoneLunaMod', mode='ozone', scale='rel', compare_index=1)

for ax in fig4.axes:
    handles, labels = ax.get_legend_handles_labels()
    labels_new = np.unique(labels)
    handles_new = [handles[i] for i in np.unique(labels, return_index=True)[1]]
    ax.legend(handles_new, labels_new, ncol=2, loc='center right')

plt.close('all')

plot_3d(variable='C:N')

# Show it
plt.show(block=False)
