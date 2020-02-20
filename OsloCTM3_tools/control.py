import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import get_vmr, set_pressure_axis, print_all

### Functions ###
def load_data(path, tracer):
    '''
    Load monthly averages from OsloCTM3 run.
    '''
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src)):
        print(file)
        data = xr.open_dataset(file)
        data = get_vmr(data,tracer=tracer, unit='ppm')
        data_list.append(data)
    # Concatinating the list
    data = xr.concat(data_list, dim='time')
    return(data)

def plot_global_avg(tracer):
    '''
    Plot some zonal average.
    '''
    
    weights = np.cos(tracer.lat*np.pi/180)
    zonal_mean = ((tracer.mean(dim='lon')*weights).mean(dim='lat').transpose())
    levels = 40

    
    fig = plt.figure(figsize=(16,9))
    fig.canvas.set_window_title("global_averages_%s" % (tracer['name'].values))

    # Global zonal weighted average
    ax1 = plt.subplot(221)
    ax1.set_title("Global zonal weighted average")
    zonal_mean.plot.contourf(ax=ax1, levels=levels)
    
    # Global zonal weighted average difference with start of simulation
    ax2 = plt.subplot(222)
    ax2.set_title("Delta global zonal weighted average")
    (zonal_mean-zonal_mean.isel(time=0)).plot.contourf(ax=ax2, levels=levels)

    for ax in fig.axes[::2]:
        ax.set_xlabel("Time (months)")

    # Zonal longitudinal average
    ax3 = plt.subplot(223)
    ax3.set_title("Zonal longitudinal average")
    (tracer.mean(dim='lon').mean(dim='time')).plot.contourf(ax=ax3, levels=levels)

    # Zonal longitudinal average difference
    ax4 = plt.subplot(224)
    ax4.set_title("Delta zonal longitudinal average end-begin")
    (tracer.mean(dim='lon').isel(time=-1)-tracer.mean(dim='lon').isel(time=0)).plot.contourf(ax=ax4, levels=levels)

    for ax in fig.axes[::2]:
        set_pressure_axis(ax)
        ax.set_title('')
    for ax in fig.axes[1::2]:  
        ax.set_ylabel("%s (%s)" % (tracer['name'].values, tracer.units))

    ax4.set_title("month_end - month_begin")


### MAIN ###
# Path to monthly means from OsloCTM3 run
#nc_src = os.environ['DATA']+'/nird_data/results/OsloCTM3/drydepdevel/version2/C3RUN_mOSaic/monthly_means/*.nc'
nc_src = os.environ['DATA']+'/nird_data/models/results/OsloCTM3/ozone25/C3RUN_ozone25_*/monthly_means/*.nc'

# Load the data only once (if run in ipython interpreter)
try:
    data
except NameError:
    data = load_data(nc_src,'CH4')

# Close plots from previous executions of the script (if run in ipython interpreter)
plt.close('all')
# Plot it
plot_global_avg(data)

# Show it (comment this out if you don't run interactivly)
plt.show(block=False)

# Print it
print_all()
