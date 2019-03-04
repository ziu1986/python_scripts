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

def hist_hist_diagram(wd, ws, title, xlabel, ylabel, fig1, **kwargs):
    model = kwargs.pop("model", True)
    if model:
        wd[np.where(wd>360)] = wd[np.where(wd>360)]-360
    bin_range_wd = range(0,361)
    bin_range_ws = np.arange(-0.05, 13.1, 0.1)
    
    # Definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    fig1.canvas.set_window_title(title)

    nullfmt = NullFormatter()         # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    axScatter.hist2d(wd, ws, bins=(bin_range_wd, bin_range_ws), cmap=plt.cm.hot_r)
    axHistx.hist(wd, bins=bin_range_wd, normed=True)
    axHisty.hist(ws, bins=bin_range_ws, normed=True, orientation='horizontal')

    # Adjust the axes ranges
    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    # Decorate the plot
    axScatter.set_xlabel(xlabel)
    axScatter.set_ylabel(ylabel)
    axScatter.set_title("%s" % (title[:title.find('_')]), y=1.15, x=1.1, size='xx-large')
    axScatter.text(87,12.7,"E", size='x-large')
    axScatter.text(177,12.7,"S", size='x-large')
    axScatter.text(267,12.7,"W", size='x-large')
    axScatter.text(0,12.7,"N", size='x-large')
    axScatter.text(353,12.7,"N", size='x-large')
    axScatter.set_xticks(np.arange(0,361,45))
    axHistx.set_xticks(np.arange(0,361,45))
    axHistx.set_ylabel("Probability density")
    axHisty.set_xlabel("Probability density")

def windrose(wd, ws, title, fig2, **kwargs):
    model = kwargs.pop('model', True)
    fig2.canvas.set_window_title(title)

    if model:
        wd[np.where(wd>360)] = wd[np.where(wd>360)]-360
    # Shift wind directions [-180:180] -> [0:360]
    ax21 = WindroseAxes.from_ax(fig=fig2)
    ax21.bar(wd, ws, normed=True, opening=0.8, edgecolor='white', cmap=plt.cm.viridis)
    ax21.set_legend()
    ax21.set_title("%s" % (title[:title.find('_')]), y=0.98, x=0., size='xx-large')


# Closing plots from previous runs
plt.close('all')

# Path to the file
#src_dir = os.environ['DATA'] + '/CTM3_input_data' + '/metdataOpenIFS/cy38r1nc4_2012/T159N80L60/*/*.nc'
src_dir = os.environ['HOME'] +"/bin/wind*1990.nc"
src_dir_obs = os.environ['DATA'] + '/astra_data/observations/wind/Svanvik_2018.txt'

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

# Observation data
data_obs = pd.read_csv(src_dir_obs,';', header=13)

# Select station
selection = "Svanvik"
u10_selection = u10.sel(lat=station_location[selection].lat, lon=station_location[selection].lon, method='nearest')
v10_selection = v10.sel(lat=station_location[selection].lat, lon=station_location[selection].lon, method='nearest')
# Plot it
wd = 270-(np.arctan2(v10, u10)/np.pi*180).data.ravel()
ws = (np.sqrt(u10**2+v10**2)).data.ravel()

wd_selection = 270-(np.arctan2(v10_selection, u10_selection)/np.pi*180).data.ravel()
ws_selection = (np.sqrt(u10_selection**2+v10_selection**2)).data.ravel()


 
xlabel = "Wind direction (deg)"
ylabel = "Wind strength ($ms^{-1}$)"

fig1 = plt.figure(1)
#hist_hist_diagram(wd, ws, "fennoscandic_wind", xlabel, ylabel, fig1)
hist_hist_diagram(wd_selection, ws_selection, "%s_wind" % selection, xlabel, ylabel, fig1)

fig2 = plt.figure(2)
#windrose(wd, ws, "fennoscandic_windrose")
windrose(wd_selection, ws_selection, "%s_windrose" % selection, fig2)


fig3 = plt.figure(3)
hist_hist_diagram(data_obs['DD'].where(data_obs['FF']>0.8), data_obs['FF'].where(data_obs['FF']>0.8), "%s_obs_wind" % selection, xlabel, ylabel, fig3, model=False)


fig4 = plt.figure(4)
windrose(data_obs['DD'].where(data_obs['FF']>0.8), data_obs['FF'].where(data_obs['FF']>0.8), "%s_obs_windrose" % selection, fig4, model=False)


# Show it
plt.show(block=False)

