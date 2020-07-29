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
        "spin-up_brazil_2000_ozone_luna_0",
        "spin-up_brazil_2000_ozone_luna_100",
        "spin-up_brazil_2000_ozone_luna_100_pwu",
        "spin-up_brazil_2000_wohydr",
        "spin-up_brazil_2000_wohydr_ozone",
        "spin-up_brazil_2000_wohydr_ozone_luna_100")

basedir = os.environ['CESM_RUN']
subdir1 = ('work', 'archive')
subdir2 = {'work':'run/', 'archive':'lnd/hist/'}
filename = "*.clm2.h0.*"

data_list = []

# Loop through all cases
for icase in case:
    # Select right directory
    try:
        src = basedir + '/' + subdir1[1] + '/' + icase + '/' + subdir2[subdir1[1]] + filename
        files = sorted(glob.glob(src))[-1]
    except IndexError:
        src = basedir + '/' + subdir1[0] + '/' + icase + '/' + subdir2[subdir1[0]] + filename
        files = sorted(glob.glob(src))[-1]

    # Extract data
    data_list.append(load_data(files, var=["NPP", "TOTVEGC", "TOTVEGN"]))

data = xr.concat(data_list, dim='case')
data.coords['case'] = np.arange(len(case))
# Plot it
fig1 = plt.figure(1)
ax11 = plt.subplot()
ax11.bar(data.mean(dim='time')['case'].data, data.mean(dim='time')['NPP'].data.flatten())

fig2 = plt.figure(2)
ax21 = plt.subplot()
ax21.bar(data.mean(dim='time')['case'].data, data.mean(dim='time')['TOTVEGC'].data.flatten())


fig3 = plt.figure(3)
ax31 = plt.subplot()
ax31.bar(data.mean(dim='time')['case'].data, data.mean(dim='time')['TOTVEGN'].data.flatten())


# Show it
plt.show(block=False)
