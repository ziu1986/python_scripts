import os, glob, sys
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl

from mytools.plot_tools import print_all, get_month_name

# Clean up
plt.close('all')

# Source
src = os.environ['DATA'] + "/astra_data/observations/soildata_svanvik/*.txt"

# Read data
try:
    data_list
except NameError:
    data_list = []
    for file in sorted(glob.glob(src)):
        data = pd.read_table(file, header=7, delim_whitespace=True, na_values=-9999)
        datetime = ["%s %s" % (data.index[i], data.iloc[:,0][i]) for i in range(len(data))]
        data.index = pd.to_datetime(datetime)
        data = data.drop(columns=[u'Tid'])

        data_list.append(data)

    data = pd.concat(data_list, axis=1)

selection = data.iloc[:,[1,3,5,7,9,11]].astype(float)
selection.columns = data.iloc[0][np.arange(0,12,2)].values
vmax = selection.where((selection < 30) & (selection >= 0)).max().max()
vmin = 0.

# Plot it
fig1 = plt.figure(1, figsize=(9,12))
fig1.canvas.set_window_title("svanhovd_soilwater_content")
for i, j in zip(range(6), range(1,7)):
    ax = plt.subplot(6,1,j)
    ax.set_title("%s cm" % data.iloc[:,i*2][0], x=0.1, y=0.8)
    z_data = selection.where((selection < 30) & (selection >= 0)).iloc[:,i]
    z_data.plot(ax=ax, marker='.', ls='None')

    ax.set_ylim(0,30)
    ax.set_xlim([dt.date(2018, 1, 1), dt.date(2018, 12, 31)])

ax.set_ylabel("Volumetric soilwater content (%)", y=3.25)
ax.set_xlabel("Time (months)")
ax.set_xticklabels([get_month_name(i, length=3) for i in range(1,13 )])

for ax in fig1.axes[:-1]:
    ax.set_xticklabels("")


fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("svanhovd_soilwater_content_3d")
Y = selection.columns.values
X = selection.index.values
Z = selection.values

x,y = np.meshgrid(X, Y)
levels = np.arange(vmin,np.round(vmax,-1), 1)
cp1 = plt.contourf(x, y, Z.transpose(), levels=levels, cmap=plt.cm.viridis_r)
cbar = plt.colorbar(cp1)

ax21 = plt.gca()
ax21.set_xlabel("Time (months)")
ax21.set_ylabel("Height (cm)")
ax21.set_xticklabels([get_month_name(i, length=3) for i in range(1,13)])
cbar.set_label("Volumetric soilwater content (%)")

# Show it
plt.show(block=False)

    
