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

# Clean up
plt.close('all')

# Data source
data_dir = os.environ['DATA']+'/astra_data/ctm_results/'
scav_dir = "scavenging_daily/"

experiment = (#'C3RUN_default/',
              'C3RUN_mOSaic/',
              #'C3RUN_mOSaic_offLight/',
              #'C3RUN_mOSaic_offPhen/',
              #'C3RUN_mOSaic_SWVL1/',
              #'C3RUN_mOSaic_ice/',
              #'C3RUN_mOSaic_desert/',
              #'C3RUN_mOSaic_emis2014/',
              #'C3RUN_mOSaic_hough/'
)


labels = (#'Wesely_type',
          'mOSaic',
          'mOSaic_offLight','mOSaic_offPhen','mOSaic_SWVL1',
          'mOSaic_ice',
          'mOSaic_desert',
          'mOSaic_emis2014',
          'mOSaic_hough')
colors = np.concatenate((#('black',),
                         ('blue',),
                         np.repeat('cornflowerblue',3),
                         np.repeat('purple',2),
                         ('grey',),
                         ('orange',)) )

try:
    data
except NameError:
    vdd_raw_data = []
    for iexp in experiment:
        subdir = data_dir+iexp+scav_dir+'scavenging_daily_2d*.nc'
        print("Reading from path %s" % (os.path.abspath(subdir)))
        data_list = []
        # Open dataset
        for file in sorted(glob.glob(subdir)):
            print("Reading %s" % (os.path.basename(file)))
            data = xr.open_dataset(file)
            # Define new time coordinates and drop the old one
            #data['time'].reset_coords(drop=True)
            year = data['YEAR']
            month = data['MONTH']
            day = data['DAY']
            data.coords['time'] = ([dt.datetime(year, month, day),])
            if (data.lat.ndim > 1):
                print("Changing coordinates...")
                data.coords['x'] = data.lat.lon[0].data
                data.coords['y'] = data.lat[:,0].data
                data = data.drop(('lon','lat'))
                data = data.rename({'x':'lon','y':'lat'})
            data_list.append(data.get(['VRaO3_avg','VRbO3_avg','VRcO3_avg']))
        # Concatenating the list
        vdd_raw_data.append(xr.concat(data_list, dim='time'))


# Plot it
for icat in vdd_raw_data[0]['NLCAT']:
    fig1 = plt.figure(figsize=(10,9))
    fig1.canvas.set_window_title("drydep_velo_nlcat_%d" % (icat))
    ax11 = plt.subplot(311)
    ax12 = plt.subplot(312)
    ax13 = plt.subplot(313)
    levels = np.arange(0.1,6.1,0.5)
    (1e2*vdd_raw_data[0])['VRaO3_avg'].sel(NLCAT=icat).mean(dim='time').plot(ax=ax11, levels=levels , cmap=plt.cm.hot_r)
    (1e2*vdd_raw_data[0])['VRbO3_avg'].sel(NLCAT=icat).mean(dim='time').plot(ax=ax12, levels=levels , cmap=plt.cm.hot_r)
    (1e2*vdd_raw_data[0])['VRcO3_avg'].sel(NLCAT=icat).mean(dim='time').plot(ax=ax13, levels=levels , cmap=plt.cm.hot_r)

    ax11.set_title('$v_{Ra}$')
    ax12.set_title('$v_{Rb}$')
    ax13.set_title('$v_{Rc}$')

    for ax in fig1.axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
    ax13.set_xlabel("Latitudes (deg)")
    ax12.set_ylabel("Longitudes (deg)")

ax11.legend()
# Show it
plt.show(block=False)
