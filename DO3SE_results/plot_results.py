import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mytools.plot_tools import print_all, plot_error_bands
from do3se_tools import *

# Clean up
plt.close('all')

# Source
src = os.environ['DATA'] + '/DO3SE_results/'
species = ("Birch", "Norway spruce", "Perennial grass")



# Read the data only once
try:
    date_list
except NameError:
    data_list = {}

    for spec in species:
        data_list.update({"%s" % spec:read_data(glob.glob(src+spec+'*')[0])})
 
# Plot it
stop_day = 300
start_day = 90
fig1 = plt.figure(1)
fig1.canvas.set_window_title("DO3SE_results")

ax11 = plt.subplot(211)
ax11.set_title("2018")
ax12 = plt.subplot(212)
ax12.set_title('2019')

for spec in species:
    data = data_list[spec]
    for ax, sheet, color in zip((ax11, ax12), data.sheet_names[1::2][:3], ('violet', 'purple')):
        date = pd.read_excel(data, sheet, header=2)
        date.index = date.index+(date['Day'].iloc[0]-1)*24
        date = date.reindex(np.arange(1,365*24))
        date['PODY (mmol/m^2 PLA)'].plot(ax=ax, linewidth=3, label=spec, color=color, use_index=False)

ax12.set_xlabel("Time (doy)")
ax12.set_ylabel('$POD_y$ ($mmol\,m^{-2}$ PLA)', y=1)
for ax in fig1.axes:
    ax.set_xlim(start_day*24,stop_day*24)
    ax.set_xticks(np.arange(start_day*24,stop_day*24,30*24))
    ax.set_xticklabels(np.arange(start_day*24,stop_day*24,30*24)/24)
    ax.legend()
    ax.set_ylim(0,15)

fig2 = plt.figure(2, figsize=(12,12))
b_delta = False
if b_delta:
    fig2.canvas.set_window_title("DO3SE_results_rel_clim")
else:
    fig2.canvas.set_window_title("DO3SE_results_rel")
ax21 = plt.subplot(311)
ax22 = plt.subplot(312)
ax23 = plt.subplot(313)
ax21.set_title("(a)", x=0.025, y=0.85)
ax22.set_title("(b)", x=0.025, y=0.85)
ax23.set_title("(c)", x=0.025, y=0.85)


#key_list = ('Gsto (mmol/m^2/s)', 'Vd (m/s)', 'f_temp', 'VPD (kPa)', 'Ts_C (C)', 'PAR (umol/m^2/s)', 'f_light', 'LAI', 'f_phen', 'O3_zR (ppb)', 'precip (mm)', 'AOT40 (ppm)', 'f_VPD', 'LWP (MPa)', 'Rsur (s/m)')
key_list = ('Gsto (mmol/m^2/s)', 'f_temp', 'VPD (kPa)', 'Ts_C (C)', 'PAR (umol/m^2/s)', 'f_light', 'f_phen', 'O3_zR (ppb)', 'precip (mm)', 'AOT40 (ppm)', 'f_VPD', 'Vd (m/s)')

correlation_pody_list = []
correlation_gsto_list = []

for ax, spec in zip((ax21, ax22, ax23), species):
    # Check species
    print(spec)
    # Load excel file for the species
    data = data_list[spec]
    # Extract climatology (disregard fSWP simulations)
    date_clim = pd.read_excel(data, data.sheet_names[1::2][2], header=2)
    # Max PODY value
    pody_max = date_clim['PODY (mmol/m^2 PLA)'].max()
    # Max PODY value for fSWP simulation
    pody_fswp_max = pd.read_excel(data, data.sheet_names[2::2][2], header=2)['PODY (mmol/m^2 PLA)'].max()
    # Check on systematic uncertainty
    print("Relative syst. uncertainty from fSWP: %1.3f" % ((pody_max-pody_fswp_max)/pody_max))
    # De-accumulate PODY and generate a DateFrame
    uncum_date_clim = uncumsum(date_clim,'PODY (mmol/m^2 PLA)')
    uncum_date_clim = pd.DataFrame({'PODY':uncum_date_clim, 'Month':date_clim['Month'], 'Doy':date_clim['Day']})
    if ~b_delta:
        #uncum_date_clim.groupby(['Doy']).agg([np.mean, np.std])['PODY'].plot(ax=ax, y = "mean", linewidth=3, label='_', color='black',ls='-')
        plot_error_bands(ax, uncum_date_clim.groupby(['Doy']).mean().index, uncum_date_clim.groupby(['Doy']).mean()['PODY'], uncum_date_clim.groupby(['Doy']).std()['PODY'], alpha=0.25, color='black')
            
    # Compute correlation coefficient for all variables (for daily means) for climatology
    # Sort dictionary: sorted(correlations(date_clim, "PODY (mmol/m^2 PLA)", uncum=True, keys=key_list).items(), key=lambda x: x[1], reverse=True)
    correlation_pody_list.append((correlations(date_clim.where((date_clim['Day']>100)&(date_clim['Day']<270)), "PODY (mmol/m^2 PLA)", uncum=True, keys=key_list, daily=True)))
    correlation_gsto_list.append((correlations(date_clim.where((date_clim['Day']>100)&(date_clim['Day']<270)), "Gsto_l (mmol/m^2/s)", keys=key_list, daily=True)))
    
    for sheet, color, year, marker in zip(data.sheet_names[1::2][:2], ('violet', 'purple'), ("2018", "2019"), ('^','v')):
        # Check year
        print(year)
        # Load excel sheet for each year
        date = pd.read_excel(data, sheet, header=2)
        # Max PODY value
        pody_max = date['PODY (mmol/m^2 PLA)'].max()
        # Max PODY value for fSWP simulation
        pody_fswp_max = pd.read_excel(data, sheet+'_SWP', header=2)['PODY (mmol/m^2 PLA)'].max()
        # Check on systematic uncertainty
        print("Relative syst. uncertainty from fSWP: %1.3f" % ((pody_max-pody_fswp_max)/pody_max))
        
        # Reindex to match climatology
        date.index = date.index+(date['Day'].iloc[0]-1)*24
        date = date.reindex(np.arange(1,365*24))

        # Compute correlation coefficience for all variables (2018->2019)
        correlation_pody_list.append((correlations(date.where((date['Day']>100)&(date['Day']<270)), "PODY (mmol/m^2 PLA)", uncum=True, keys=key_list, daily=True)))
        correlation_gsto_list.append((correlations(date.where((date['Day']>100)&(date['Day']<270)), "Gsto_l (mmol/m^2/s)", keys=key_list, daily=True)))
        
        # De-accumulate and create DataFrame
        uncum_date = pd.DataFrame({'PODY':uncumsum(date, 'PODY (mmol/m^2 PLA)'), 'Month':date['Month'], 'Doy':date['Day']})

        # Compute difference with climatology
        delta_uncum_date = (uncum_date-uncum_date_clim)
        #
        # Compute deviation from monthly standart deviation for each hour 
        #delta_siguncum_date = delta_uncum_date[['PODY']].div(uncum_date_clim.groupby('Month').transform('std')).join(date_clim['Month']) 
        
        if b_delta:
            delta_uncum_date['PODY'].plot(ax=ax, linewidth=3, label=year, color=color, use_index=False)
        else:
            # Add mean and std to table using aggiate function and plot with errorbars
            uncum_date.groupby(['Doy']).agg([np.mean, np.std])['PODY'].plot(ax=ax, y = "mean", yerr = "std", label=year, color=color,ls='None', marker=marker, fillstyle='none')

ax23.set_xlabel("Time (doy)")
ax22.set_ylabel('$\Delta POD_y / \Delta t$ ($mmol\,m^{-2}\,s^{-1}$)')
for ax in fig2.axes[:-1]:
    ax.set_xlabel('')
    ax.set_xticklabels(())
    
for ax in fig2.axes:
    if b_delta:
        ax.set_xlim(start_day*24,stop_day*24)
        ax.set_xticks(np.arange(start_day*24,stop_day*24,30*24))
        ax.set_xticklabels(np.arange(start_day*24,stop_day*24,30*24)/24)
        ax.set_ylim(-0.02,0.02)
    else:
        ax.set_xlim(start_day, stop_day)
        ax.set_xticks(np.arange(start_day, stop_day, 30))
        ax.set_ylim(0,0.02)
    ax.legend()


for spec, figi  in zip(species, np.arange(3,6)):
    fig = plt.figure(figi, figsize=(16,6))
    fig.canvas.set_window_title("DO3SE_results_pody_gsto_o3_%s" % spec.replace(' ', '_'))

    data = data_list[spec]
    for iax, sheet, color, ititle in zip(np.arange(1,10), data.sheet_names[1::2][:3], ('violet', 'purple', 'blueviolet'), char_range('a', 'c')):
        ax = plt.subplot(1,3,iax)
        ax.set_title("(%s)" % ititle)

        date = pd.read_excel(data, sheet, header=2)
        date.index = date.index+(date['Day'].iloc[0]-1)*24
        date = date.reindex(np.arange(1,365*24))
        # Plot data
        plot_pody_gsto_o3(ax, date, o3color=color)

    for ax in fig.axes[1:-2]:
        ax.set_ylabel("")
        ax.set_yticklabels("")

# Plot correlation coefficients
corr_birch = pd.DataFrame(correlation_pody_list[:3], index=('clim', '2018', '2019'))
corr_spruce = pd.DataFrame(correlation_pody_list[3:6], index=('clim', '2018', '2019'))
corr_grass = pd.DataFrame(correlation_pody_list[6:], index=('clim', '2018', '2019'))

fig6 = plt.figure(6, figsize=(12,12))
fig6.canvas.set_window_title("DO3SE_results_pody_corr")
ax61 = plt.subplot(311)
ax62 = plt.subplot(312)
ax63 = plt.subplot(313)

corr_birch.transpose().plot.bar(ax=ax61, color=('blueviolet','violet','purple'))
corr_spruce.transpose().plot.bar(ax=ax62, color=('blueviolet','violet','purple'))
corr_grass.transpose().plot.bar(ax=ax63, color=('blueviolet','violet','purple'))

for ax, ititle in zip(fig6.axes, char_range('a','c')):
    ax.set_ylim(-0.7,1)
    ax.set_title("(%s)" % ititle, x=0.025, y=0.875)
    ax.legend(loc='lower right', ncol=3)
    
ax61.set_xticklabels(())
ax62.set_xticklabels(())

ax63.set_xticklabels([label.get_text()[:label.get_text().find(' ')] for label in ax63.get_xticklabels()])
    
ax62.set_ylabel("$\\rho$")

# Show it
plt.show(block=False)

