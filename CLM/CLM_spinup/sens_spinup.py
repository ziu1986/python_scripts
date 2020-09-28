import os, sys, glob, calendar
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all
from mytools.clm_tools import *

# Clean up
plt.close('all')

case = ("spin-up_brazil_2000",
        "spin-up_brazil_2000_ozone",
        "spin-up_brazil_2000_ozone_luna_100",
        "spin-up_brazil_2000_ozone_luna_100_pwu",
        "spin-up_brazil_2000_wohydr",
        "spin-up_brazil_2000_wohydr_ozone",
        "spin-up_brazil_2000_wohydr_ozone_luna_100")

basedir = os.environ['CESM_RUN']
subdir1 = ('work', 'archive')
subdir2 = {'work':'run/', 'archive':'lnd/hist/'}
filename = "*.clm2.h0.*"

vars = ["GPP", "NPP", "TLAI", "TOTVEGC", "TOTVEGN"]
case_name = ["ctrl."] + ["OzoneMod", "OzoneLunaMod par1", "OzoneLunaMod par2", "ctrl. w/o hydr.", "OzoneMod w/o hydr.", "OzoneLunaMod w/o hydr."] #[icase.replace("spin-up_brazil_2000_", "") for icase in case[1:]]

colors = ['#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b']
data_list = []

# Loop through all cases
for icase in case:
    # Select right directory
    try:
        src = basedir + '/' + subdir1[1] + '/' + icase + '/' + subdir2[subdir1[1]] + filename
        files = sorted(glob.glob(src))[-1]
        print(subdir1[1], icase)
    except IndexError:
        src = basedir + '/' + subdir1[0] + '/' + icase + '/' + subdir2[subdir1[0]] + filename
        files = sorted(glob.glob(src))[-1]
        print(subdir1[1], icase)

    # Extract data
    data_list.append(load_data(files, var=vars))

data = xr.concat(data_list, dim='case')
data.coords['case'] = np.arange(len(case))

# Generate pandas series for bar plotting
pd_data = (data.mean(dim='time')/data.mean(dim='time').sel(case=0)-1).to_dataframe()
pd_data.index = case_name


# Plot it
fig1 = plt.figure(1, figsize=(10,8))
fig1.canvas.set_window_title("sens_spinup")
ax11 = plt.subplot()
(pd_data.transpose().iloc[:,1:]*100).plot.bar(ax=ax11, rot=0, width=0.95, color=colors[2:])
for p in ax11.patches:
    if p.get_height() < 0:
        color='white'
    else:
        color='black'
    
    ax11.annotate("%d" % (np.abs(np.round(p.get_height(),2))), (p.get_x(), p.get_height()), size='x-large', color=color)


ax11.legend(ncol=3)
ax11.set_ylim(-35, 35)
ax11.set_ylabel("$(X_{case}/X_{ctrl})_i -1 $ (%)")


# Show it
plt.show(block=False)
