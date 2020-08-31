import os, glob, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from mytools.plot_tools import print_all

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
        data = pd.read_table(file, header=7, delim_whitespace=True)
        datetime = ["%s %s" % (data.index[i], data.iloc[:,0][i]) for i in range(len(data))]
        data.index = pd.to_datetime(datetime)
        data = data.drop(columns=[u'Tid'])

        data_list.append(data)

    data = pd.concat(data_list, axis=1)


# Plot it
fig1 = plt.figure(1, figsize=(9,12))
fig1.canvas.set_window_title("svanhovd_soilwater_content")
for i, j in zip((1,3,4,5), range(1,5)):
    ax = plt.subplot(4,1,j)
    ax.set_title("%s cm" % data.iloc[:,2*i][0], x=0.1, y=0.8)
    data.iloc[:,2*i+1].plot(ax=ax)

    ax.set_ylim(0,20)

ax.set_ylabel("Volumetric soilwater content", y=2.25)
ax.set_xlabel("Time (months)")

for ax in fig1.axes[:-1]:
    ax.set_xticklabels("")
# Show it
plt.show(block=False)

    
