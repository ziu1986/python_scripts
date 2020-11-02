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

# Plot it
stop_day = 300
start_day = 90
fig1 = plt.figure(1)
selection = data.sheet_names[1][:data.sheet_names[1].find('_')]
fig1.canvas.set_window_title("DO3SE_results")
ax11 = plt.subplot(211)
ax11.set_title("2018")
ax12 = plt.subplot(212)
ax12.set_title('2019')

for spec in species:
    data = read_data(glob.glob(src+spec+'*')[0])
    for ax, sheet, color in zip((ax11, ax12), data.sheet_names[1::2][:-2], ('violet', 'purple')):
        date = pd.read_excel(data, sheet, header=2)
        date.index = date.index+(date['Day'].iloc[0]-1)*24
        date = date.reindex(np.arange(1,365*24))
        date['PODY (mmol/m^2 PLA)'].plot(ax=ax, linewidth=3, label=spec, color=color, use_index=False)

ax12.set_xlabel("Time (doy)")
ax12.set_ylabel('$POD_y$ ($mmol m^{-2}$ PLA)', y=1)
for ax in fig1.axes:
    ax.set_xlim(start_day*24,stop_day*24)
    ax.set_xticks(np.arange(start_day*24,stop_day*24,30*24))
    ax.set_xticklabels(np.arange(start_day*24,stop_day*24,30*24)/24)
    ax.legend()
    ax.set_ylim(0,15)
# Show it
plt.show(block=False)

