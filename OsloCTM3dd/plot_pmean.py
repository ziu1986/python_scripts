import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *
from collections import OrderedDict


# Data source
data_dir = os.environ['HOME']+"/bin/P_Mean/*"
p_mean_raw_data = []
try:
    p_mean
except NameError:
    # Read the files
    for file in sorted(glob.glob(data_dir)):
        p_mean_raw_data.append(read_data(file))

    p_mean = xr.concat(p_mean_raw_data,dim='time')

# Clean up
plt.close('all')
# Plot it
fig1 = plt.figure(1)
fig1.canvas.set_window_title("pmean-variation")
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)
p_mean.mean(dim='time')["Pmean"].plot(ax=ax11, cmap=plt.cm.Spectral_r)
(p_mean.std(dim='time')["Pmean"]/p_mean.mean(dim='time')["Pmean"]).plot(ax=ax12)

for ax in fig1.axes:
    ax.set_xlabel("")
    ax.set_xticks(np.arange(0,361,60))
    ax.set_ylabel("")
    ax.set_yticks(np.arange(-90,91,45))
ax12.set_xlabel("Longitude (deg)")
ax12.set_ylabel("Latitude (deg)", y=1)
fig1.axes[-2].set_ylabel("$<P_{sfc}> (hPa)$")
fig1.axes[-1].set_ylabel("$\sigma_{P_{sfc}} / <P_{sfc}> (hPa)$")


# Show it
plt.show(block=False)


    
