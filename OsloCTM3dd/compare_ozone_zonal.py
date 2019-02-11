import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Data source
nc_bo = os.environ['DATA']+'/BodekerScientific/BSVerticalOzone_ND_PRS_Tier0.5_V1.0.nc'
nc_src = os.environ['DATA']+'/astra_data/ctm_results/C3RUN_emep_ppgs_2005/moleccm3/mean_ozone2005.nc'
try:
    data
except NameError:
    # Open dataset
    data = xr.open_dataset(nc_src,decode_times=False)
    # Define new time coordinates and drop the old one
    print("Setting time...")
    date = [dt.datetime(int(nc_src[-7:-3]), 1, 1) + dt.timedelta(each) for each in np.arange(len(data['time']))]
    data['time'].reset_coords(drop=True)
    data.coords['time'] = (date)
    data_bo = xr.open_dataset(nc_bo)
       
# Concatinating the list
ozone_zonalmean = data.mean(dim='lon')
ozone_zonalmean.attrs['units'] = "molec/m3"

# Plotting
plt.close('all')

fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(122)
ax12 = plt.subplot(121)
#levels = np.arange(0,6.1e12,5e11)
levels = np.concatenate((np.arange(0,1.2e12,1e11),np.arange(2.5e12,6.1e12,5e11)))
(data['O3'].mean(dim='time').mean(dim='lon')*N_A*1e-6).plot(ax=ax11, levels=levels) #.contourf
ax11.set_title("EMEP_ppgs_2005")
#ax11.axhspan(1000,data_bo.level.max().data,alpha=0.75,color='white')
data_bo['o3_mean'].where(data_bo.time.dt.year==2005).mean(dim='time').transpose().plot(ax=ax12, levels=levels) #.contourf
ax12.set_title("BSVerticalOzone")

for ax in fig1.axes[:2]:
    ax.invert_yaxis()
    ax.set_ylim(1000,10)
    ax.set_ylabel("")
    ax.set_xlabel("Latitude (deg)")
for ax in fig1.axes[2:]:
    ax.set_ylabel("")

ax12.set_ylabel("Pressure (hPa)")
fig1.axes[-2].set_ylabel("[$O_3$] (molec cm$^{-3}$)")

for each in data_bo.level:
    ax11.axhline(each,ls=':', color='black', alpha=0.5)
    ax12.axhline(each,ls=':', color='black', alpha=0.5)

# Show it
plt.show(block=False)
