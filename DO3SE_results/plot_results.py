import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mytools.plot_tools import print_all

# Clean up
plt.close('all')

# Source
src = os.environ['DATA'] + '/DO3SE_results/'
species = ("Birch", "Norway spruce", "Perennial grass")

def read_data(src):
    # Read as Excel object
    data = pd.ExcelFile(src)
    return(data)

def plot_pody_gsto_o3(ax, date, **karg):
    o3_color = karg.pop('o3color', 'blueviolet')
    # Define additional axis
    ax2 = ax.twinx()
    ax3 = ax.twinx()

    date['Gsto_l (mmol/m^2/s)'].plot(ax=ax, ls='None', marker='o', label='$G_{sto}^{leaf}$')
    date['O3_zR (ppb)'].plot(ax=ax2, ls='None', marker='x', color=o3_color, alpha=0.5, label='$[O_3]$')
    date['PODY (mmol/m^2 PLA)'].plot(ax=ax3, ls='None', marker='s', color='darkgreen', label='$POD_y$')
    
    # Adjust axis
    ax.set_ylim(0,250)
    ax.set_ylabel('$G_{sto}^{leaf}$ ($mmol\,m^{-2}s^{-1})$')
    
    ax2.set_ylim(0,250)
    ax2.set_yticklabels("")
    ax2.set_yticks(())
    
    ax3.set_ylim(0,16)
    ax3.set_ylabel('$POD_y$ ($mmol\,m^{-2}$ PLA)')
    ax3.grid(b=False)

    ax.set_xlabel("Time (doy)")
    ax.set_xlim(start_day*24,stop_day*24)
    ax.set_xticks(np.arange(start_day*24,stop_day*24,30*24))
    ax.set_xticklabels(np.arange(start_day*24,stop_day*24,30*24)/24)

    # Legend
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    lines3, labels3 = ax3.get_legend_handles_labels()
    ax.legend(lines + lines2 + lines3, labels + labels2 + labels3, loc="upper left")

def uncumsum(date, var):
    test = date[var].values
    # Undo cum sum and return new pandas series
    return(pd.Series(np.append(test[0], (test[2:]-test[1:-1])), index=date.index[:-1]))

def biomass(uptake, **karg):
    species = karg.pop('species', 'Norway spruce')
    b_above_gr = karg.pop('above_ground', False)

    biomass_tot = {'Birch': (100.2, 0.93), 'Norway spruce': (99.8, 0.22), "Perennial grass": (94.7, 0.62)}
    biomass_above_gr = {'Perennial grass': (93.9, 0.99)}


    if (b_above_gr) & (species == 'Perennial grass'):
        return(np.ones_like(uptake)*biomass_above_gr[species][0]-biomass_above_gr[species][1]*uptake)
    else:
        return(np.ones_like(uptake)*biomass[species][0]-biomass[species][1]*uptake)

def correlations(date, var, **karg):
    '''
    Compute correlation for each veriable in DO3SE model.
    
    '''
    b_uncum = karg.pop('uncum', False)
    keys = karg.pop('keys', date.keys())
    output = {}
    if b_uncum:
        date_var = uncumsum(date, var)
    else:
        date_var = date[var]
    for each in keys:
        tmp = date_var.corr(date[each])
        if ~np.isnan(tmp):
            #print(each, date_var.corr(date[each]))
            output.update({each: tmp})

    return(output)


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
    data = data_list[spec]
    date_clim = pd.read_excel(data, data.sheet_names[1::2][2], header=2)

    # Compute correlation coefficience for all variables
    # ASort dictionary: sorted(correlations(date_clim, "PODY (mmol/m^2 PLA)", uncum=True, keys=key_list).items(), key=lambda x: x[1], reverse=True)
    correlation_pody_list.append((correlations(date_clim.where((date_clim['Day']>100)&(date_clim['Day']<270)), "PODY (mmol/m^2 PLA)", uncum=True, keys=key_list)))
    correlation_gsto_list.append((correlations(date_clim.where((date_clim['Day']>100)&(date_clim['Day']<270)), "Gsto_l (mmol/m^2/s)", keys=key_list)))
    print(spec)
    for sheet, color, year in zip(data.sheet_names[1::2][:2], ('violet', 'purple'), ("2018", "2019")):
        print(year)
        date = pd.read_excel(data, sheet, header=2)
        date.index = date.index+(date['Day'].iloc[0]-1)*24
        date = date.reindex(np.arange(1,365*24))

        # Compute correlation coefficience for all variables
        correlation_pody_list.append((correlations(date.where((date['Day']>100)&(date['Day']<270)), "PODY (mmol/m^2 PLA)", uncum=True, keys=key_list)))
        correlation_gsto_list.append((correlations(date.where((date['Day']>100)&(date['Day']<270)), "Gsto_l (mmol/m^2/s)", keys=key_list)))
        
        uncum_date = pd.DataFrame({'FstoY':uncumsum(date, 'PODY (mmol/m^2 PLA)'), 'Month':date['Month']})

        uncum_date_clim = uncumsum(date_clim,'PODY (mmol/m^2 PLA)')
                
        delta_uncum_date = (uncumsum(date, 'PODY (mmol/m^2 PLA)')-uncum_date_clim)
        delta_uncum_date = pd.DataFrame({'FstoY':delta_uncum_date, 'Month':date_clim['Month']})

        uncum_date_clim = pd.DataFrame({'FstoY':uncum_date_clim, 'Month':date_clim['Month']})
        
        delta_siguncum_date = delta_uncum_date[['FstoY']].div(uncum_date_clim.groupby('Month').transform('std')).join(date_clim['Month']) 
        print(">+1sigma:") 
        print(delta_siguncum_date.where(delta_siguncum_date['FstoY']>1).groupby('Month').count()/(delta_siguncum_date.groupby('Month').count())*100)
        print("<-1sigma:")
        print(delta_siguncum_date.where(delta_siguncum_date['FstoY']<-1).groupby('Month').count()/(delta_siguncum_date.groupby('Month').count())*100)

        if b_delta:
            delta_uncum_date['FstoY'].plot(ax=ax, linewidth=3, label=year, color=color, use_index=False)
        else:
            uncum_date['FstoY'].plot(ax=ax, linewidth=3, label=year, color=color,ls='None', marker='o', fillstyle='none', use_index=False)

ax23.set_xlabel("Time (doy)")
ax22.set_ylabel('$\Delta POD_y$ ($mmol\,m^{-2}\,s^{-1}$)')
for ax in fig2.axes:
    ax.set_xlim(start_day*24,stop_day*24)
    ax.set_xticks(np.arange(start_day*24,stop_day*24,30*24))
    ax.set_xticklabels(np.arange(start_day*24,stop_day*24,30*24)/24)
    ax.legend()
    if b_delta:
        ax.set_ylim(-0.02,0.02)
    

fig3 = plt.figure(3, figsize=(16,6))
fig4 = plt.figure(4, figsize=(16,6))
fig5 = plt.figure(5, figsize=(16,6))

for spec, fig in zip(species, (fig3, fig4, fig5)):
    fig.canvas.set_window_title("DO3SE_results_pody_gsto_o3_%s" % spec.replace(' ', '_'))
    ax1, ax2, ax3 = fig.subplots(1,3)
    ax1.set_title("(a)")
    ax2.set_title("(b)")
    ax3.set_title("(c)")

    data = data_list[spec]
    for ax, sheet, color in zip(fig.axes, data.sheet_names[1::2][:3], ('violet', 'purple', 'blueviolet')):
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

for ax, num in zip(fig6.axes, ('(a)', '(b)', '(c)')):
    ax.set_ylim(-0.7,1)
    ax.set_title(num, x=0.025, y=0.875)
    ax.legend(loc='lower right', ncol=3)
    #ax.axhline(0.8, color='grey', ls='--')
    
ax61.set_xticklabels(())
ax62.set_xticklabels(())

ax63.set_xticklabels([label.get_text()[:label.get_text().find(' ')] for label in ax63.get_xticklabels()])
    
ax62.set_ylabel("$\\rho$")



"""
ax32 = plt.subplot(312)
ax33 = plt.subplot(313)

for spec,ax in zip(species, (ax31,ax32,ax33)):
   data = data_list[spec]
   gmax = pd.read_excel(data, data.sheet_names[0],header=2).iloc[14,1]
   for sheet, color in zip(data.sheet_names[1::2][:3], ('violet', 'purple', 'grey')):
       date = pd.read_excel(data, sheet, header=2)
       date.index = date.index+(date['Day'].iloc[0]-1)*24
       date = date.reindex(np.arange(1,365*24))
       date = date.where((date['Day']>=110)&(date['Day']<=270)).dropna()
       (1-date[u'Gsto (mmol/m^2/s)']/gmax).plot.hist(ax=ax, bins=np.arange(0,1.1,0.1), linewidth=3, label=spec, color=color, histtype='step', use_index=False)
"""
# Show it
plt.show(block=False)

