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

fig2 = plt.figure(2)
fig2.canvas.set_window_title("DO3SE_results_rel")

ax21 = plt.subplot(211)
ax22 = plt.subplot(212)

for spec in species:
    data = data_list[spec]
    for ax, sheet, color in zip((ax21, ax22), data.sheet_names[1::2][:3], ('violet', 'purple')):
        date_clim = pd.read_excel(data, data.sheet_names[2::2][2], header=2)
        date = pd.read_excel(data, sheet, header=2)
        date.index = date.index+(date['Day'].iloc[0]-1)*24
        date = date.reindex(np.arange(1,365*24))
        (date['PODY (mmol/m^2 PLA)']-date_clim['PODY (mmol/m^2 PLA)']).plot(ax=ax, linewidth=3, label=spec, color=color, use_index=False)

ax22.set_xlabel("Time (doy)")
ax22.set_ylabel('$POD_y$ ($mmol\,m^{-2}$ PLA)', y=1)
for ax in fig2.axes:
    ax.set_xlim(start_day*24,stop_day*24)
    ax.set_xticks(np.arange(start_day*24,stop_day*24,30*24))
    ax.set_xticklabels(np.arange(start_day*24,stop_day*24,30*24)/24)
    ax.legend()
    ax.set_ylim(-6,6)


fig3 = plt.figure(3, figsize=(16,9))
ax31 = plt.subplot()
ax32 = ax31.twinx()
ax33 = ax31.twinx()

date['Gsto_l (mmol/m^2/s)'].plot(ax=ax31, ls='None', marker='o', label='$G_{sto}^{leaf}$')
date['PODY (mmol/m^2 PLA)'].plot(ax=ax32, ls='None', marker='s', color='darkgreen', label='$POD_y$')
date['O3_zR (ppb)'].plot(ax=ax33, ls='None', marker='x', color='blueviolet', alpha=0.5, label='$[O_3]$')

ax31.set_ylim(0,250)
ax31.set_ylabel('$G_{sto}^{leaf}$ ($mmol\,m^{-2}s^{-1})$')
ax32.set_ylim(0,16)
ax32.set_ylabel('$POD_y$ ($mmol\,m^{-2}$ PLA)')
ax33.set_ylim(0,250)
ax33.set_yticklabels("")
#ax33.set_yticks()

ax31.set_xlabel("Time (doy)")
ax31.set_xlim(start_day*24,stop_day*24)
ax31.set_xticks(np.arange(start_day*24,stop_day*24,30*24))
ax31.set_xticklabels(np.arange(start_day*24,stop_day*24,30*24)/24)

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

