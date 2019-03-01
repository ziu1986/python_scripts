import os, glob, sys
import numpy as np
import pandas as pd               # Analysis of timeseries data
import xarray as xr               # Reading netcdf files and dealing with global fields
import matplotlib.pyplot as plt   # Plotting data
from matplotlib.ticker import NullFormatter
from scipy.constants import *     # Get physics constants
import datetime as dt             # Time-objects and maipulation
from windrose import WindroseAxes # Drawing wind roses
from mytools.met_tools import print_all
from station_info import station_location

def hist_hist_diagram(wd, ws, title, xlabel, ylabel):
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    fig1 = plt.figure(1)
    fig1.canvas.set_window_title(title)

    nullfmt = NullFormatter()         # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    axScatter.hist2d(wd, ws, bins=(range(-180,180),np.arange(0,10.1,0.1)), cmap=plt.cm.hot_r)
    axHistx.hist(wd, bins=range(-180,180), normed=True)
    axHisty.hist(ws, bins=np.arange(0,25,0.1), normed=True, orientation='horizontal')

    # Adjust the axes ranges
    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    # Decorate the plot
    axScatter.set_xlabel(xlabel)
    axScatter.set_ylabel(ylabel)
    axHistx.set_ylabel("Probability density")
    axHisty.set_xlabel("Probability density")

def windrose(wd, ws, title):
    fig2 = plt.figure(2)
    fig2.canvas.set_window_title(title)
    # Shift wind directions [-180:180] -> [0:360]
    wd[np.where(wd<0)] = 360+wd[np.where(wd<0)]
    ax21 = WindroseAxes.from_ax(fig=fig2)
    ax21.bar(wd, ws, normed=True, opening=0.8, edgecolor='white', cmap=plt.cm.viridis)
    ax21.set_legend()
    
# Closing plots from previous runs
plt.close('all')

# Path to the file
#src_dir = os.environ['DATA'] + '/CTM3_input_data' + '/metdataOpenIFS/cy38r1nc4_2012/T159N80L60/*/*.nc'
src_dir = os.environ['HOME'] +"/bin/wind*1990.nc"

# Read the data only once in interactive mode
try:
    data
except NameError:
    #u10_fenno = []
    #v10_fenno = []
    for each in sorted(glob.glob(src_dir)):
        print("Reading file %s" % (each))
        data = xr.open_dataset(each)
        year = int(each[-7:-3])
        time = pd.date_range(start='%d-01-01' % (year), end='%d-01-01' % (year+1), freq='3H')[:-1]
        lons = np.linspace(data.lon[0], data.lon[-1], data.lon.size)
        # Selecting 10 m winds      
        u10 = data['U10M']
        v10 = data['V10M']
        # Correcting the longitudinal coordinates
        u10.coords['lon'] = lons
        u10.coords['time'] = time
        v10.coords['lon'] = lons
        v10.coords['time'] = time
        # Select an area
        u10 = u10.sel(lat=slice(67.6,71.4),lon=slice(19.4,31.4))
        v10 = v10.sel(lat=slice(67.6,71.4),lon=slice(19.4,31.4))
selection = "Pallas"
u10_selection = u10.sel(lat=station_location[selection].lat, lon=station_location[selection].lon, method='nearest')
v10_selection = v10.sel(lat=station_location[selection].lat, lon=station_location[selection].lon, method='nearest')
# Plot it
X = (np.arctan2(v10, u10)/np.pi*180).data.ravel()
Y = (np.sqrt(u10**2+v10**2)).data.ravel()
wd_selection = (np.arctan2(v10_selection, u10_selection)/np.pi*180).data.ravel()
ws_selection = (np.sqrt(u10_selection**2+v10_selection**2)).data.ravel()

xlabel = "Wind direction (deg)"
ylabel = "Wind strength ($ms^{-1}$)"

#hist_hist_diagram(X, Y, "fennoscandic_wind", xlabel, ylabel)
hist_hist_diagram(wd_selection, ws_selection, "%s_wind" % selection, xlabel, ylabel)
#windrose(X, Y, "fennoscandic_windrose")
windrose(wd_selection, ws_selection, "%s_windrose" % selection)

# Show it
plt.show(block=False)

