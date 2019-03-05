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
    # Shift wind directions
    if model:
        wd[np.where(wd>360)] = wd[np.where(wd>360)]-360
   
    ax21 = WindroseAxes.from_ax(fig=fig2)
    ax21.bar(wd, ws, normed=True, opening=0.8, edgecolor='white', cmap=plt.cm.viridis)
    ax21.set_legend()
    ax21.set_title("%s" % (title[:title.find('_')]), y=0.98, x=0., size='xx-large')


# Closing plots from previous runs
plt.close('all')

# Path to the file
#src_dir = os.environ['DATA'] + '/CTM3_input_data' + '/metdataOpenIFS/cy38r1nc4_2012/T159N80L60/*/*.nc'
src_dir = os.environ['HOME'] +"/bin/wind*200[0-8].nc"
src_dir_obs = os.environ['DATA'] + '/astra_data/observations/wind/Svanvik/Svanvik_*.txt'
obs_header = {'Svanvik':{2009:16,2010:13,2011:13,2012:13,2013:13,2014:13,2015:13,2016:13,2017:13,2018:13}}
obs_footer = {'Svanvik':{2009:8,2010:2,2011:2,2012:2,2013:2,2014:2,2015:2,2016:2,2017:2,2018:2}}
# Read the data only once in interactive mode
try:
    data
except NameError:
    u10_list = []
    v10_list = []
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
        u10_list.append(u10.sel(lat=slice(67.6,71.4),lon=slice(19.4,31.4)))
        v10_list.append(v10.sel(lat=slice(67.6,71.4),lon=slice(19.4,31.4)))
    u10 = xr.concat(u10_list, dim='time')
    v10 = xr.concat(v10_list, dim='time')
# Observation data
try:
    data_obs
except NameError:
    data_obs = []
    for each in sorted(glob.glob(src_dir_obs)):
        year = int(each[each.rfind('_')+1:-4])
        data = pd.read_csv(each,';', header=obs_header['Svanvik'][year], skipfooter=obs_footer['Svanvik'][year], na_values=(-9999,-9999.0,), parse_dates=[[1,2,3,4]], date_parser=lambda y,m,d,h : pd.datetime(int(y), int(m), int(d), int(h)-1))
        data = data.set_index(['Year_Mnth_Date_Time(NMT)'])
        data_obs.append(data) #
        
    data_obs = pd.concat(data_obs)

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
windrose(data_obs.dropna().where(data_obs.dropna()['FF']>0.8).dropna()['DD'], data_obs.dropna().where(data_obs.dropna()['FF']>0.8).dropna()['FF'], "%s_obs_windrose" % selection, fig4, model=False)


# Show it
plt.show(block=False)

