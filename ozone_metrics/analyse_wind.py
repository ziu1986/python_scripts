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
    ff_threashold = kwargs.pop('ff_threashold',0)
    
    if model:
        wd[np.where(wd>360)] = wd[np.where(wd>360)]-360
    bin_range_wd = np.arange(-0.5, 361)
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
    axScatter.set_title("%s" % (title[:title.find('_')]), y=1.3, x=1.1, size='xx-large')
    if model:
        axHistx.text(1.05, 0.8,"%s" % ('model'), size='x-large', transform=axHistx.transAxes)
    else:
        axHistx.text(1.05, 0.8,"%s" % ('observation'), size='x-large', transform=axHistx.transAxes)

        
    axScatter.set_xlabel(xlabel)
    axScatter.set_ylabel(ylabel)
    axScatter.text(87,12.7,"E", size='x-large')
    axScatter.text(177,12.7,"S", size='x-large')
    axScatter.text(267,12.7,"W", size='x-large')
    axScatter.text(0,12.7,"N", size='x-large')
    axScatter.text(353,12.7,"N", size='x-large')
    axScatter.set_xticks(np.arange(0,361,45))
    axHistx.set_xticks(np.arange(0,361,45))
    axHistx.set_ylabel("Probability density")
    axHisty.set_xlabel("Probability density")
    if ff_threashold>0:
        axHisty.axhspan(0,ff_threashold, facecolor='grey', hatch='///')
        axScatter.axhspan(0,ff_threashold, facecolor='grey', hatch='///')
    # Add some stats
    axHistx.text(1.05, 0.65, "$N = %d$" % (wd.size), size='x-large', transform=axHistx.transAxes)
    if ff_threashold:
        axScatter.text(-25, ff_threashold, "%1.2f" % (ff_threashold), size='x-large')
    
def windrose(wd, ws, title, fig2, **kwargs):
    model = kwargs.pop('model', True)
    fig2.canvas.set_window_title(title)
    # Shift wind directions
    if model:
        wd[np.where(wd>360)] = wd[np.where(wd>360)]-360
           
    ax21 = WindroseAxes.from_ax(fig=fig2)
    ax21.bar(wd, ws, bins=np.arange(-0.05,14,2), normed=True, opening=0.8, edgecolor='white', cmap=plt.cm.viridis)
    ax21.set_legend()
    ax21.set_title("%s" % (title[:title.find('_')]), y=0.98, x=0., size='xx-large')

    if model:
        ax21.text(0.05, 0.95,"%s" % ('model'), size='x-large', transform=ax21.transAxes)
    else:
        ax21.text(0.05, 0.95,"%s" % ('observation'), size='x-large', transform=ax21.transAxes)

def test_wind_treashold(data, **kwargs):
    start = kwargs.pop('start',0)
    end = kwargs.pop('end',2)
    steps = kwargs.pop('stepwidth', 0.5)
    result = []
    for min_wind_speed in np.arange(start, end, steps):
        hist_wd = np.histogram(data['DD'].where(data['FF']>=min_wind_speed), bins=np.arange(-0.5,361), normed=True)
        bin_0 = hist_wd[0][0]
        bin_360 = hist_wd[0][-1]
        nn_0 = hist_wd[0][1:11]
        nn_360 = hist_wd[0][-11:-1]
        result.append((bin_0-nn_0.mean(), bin_360-nn_360.mean()))
    
    return(np.arange(start, end, steps), (np.array(result).T.reshape(2,len(result))))

# Closing plots from previous runs
plt.close('all')
#*************************************************************Begin of script*********************************************
# Path to the file
selection = "Karasjok"
#src_dir = os.environ['DATA'] + '/CTM3_input_data' + '/metdataOpenIFS/cy38r1nc4_2012/T159N80L60/*/*.nc'
src_dir = os.environ['HOME'] +"/bin/wind*.nc"
src_dir_obs = os.environ['DATA'] + '/astra_data/observations/wind/%s/%s_*.txt' % (selection, selection)
obs_header = {'Svanvik':{2009:16,'default':13},
              'Karasjok':{2004:16, 'default':13}}#2010:13,2011:13,2012:13,2013:13,2014:13,2015:13,2016:13,2017:13,2018:13}}
obs_footer = {'Svanvik':{2009:8, 'default':2},
              'Karasjok':{2004:8, 'default':2}}#2010:2,2011:2,2012:2,2013:2,2014:2,2015:2,2016:2,2017:2,2018:2}}
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
        if year == 2016:
            time = time.where(~((time>=time[28*24/3+31*24/3])&(time<time[29*24/3+31*24/3]))).dropna()
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
    data_ob
except NameError:
    data_obs = []
    for each in sorted(glob.glob(src_dir_obs)):
        year = int(each[each.rfind('_')+1:-4])
        if each==sorted(glob.glob(src_dir_obs))[0]:
            data = pd.read_csv(each,';', header=obs_header[selection][year], skipfooter=obs_footer[selection][year], na_values=(-9999,-9999.0,), parse_dates=[[1,2,3,4]], date_parser=lambda y,m,d,h : pd.datetime(int(y), int(m), int(d), int(h)-1))
        else:
            data = pd.read_csv(each,';', header=obs_header[selection]['default'], skipfooter=obs_footer[selection]['default'], na_values=(-9999,-9999.0,), parse_dates=[[1,2,3,4]], date_parser=lambda y,m,d,h : pd.datetime(int(y), int(m), int(d), int(h)-1))
        data = data.set_index(['Year_Mnth_Date_Time(NMT)'])
        data_obs.append(data) #
        
    data_obs = pd.concat(data_obs)
    if selection=='Svanvik':
        data_obs = data_obs['2012-02-22':]

# Select station

u10_selection = u10.sel(lat=station_location[selection].lat, lon=station_location[selection].lon, method='nearest')
v10_selection = v10.sel(lat=station_location[selection].lat, lon=station_location[selection].lon, method='nearest')
# Plot it
wd = 270-(np.arctan2(v10, u10)/np.pi*180).data.ravel()
ws = (np.sqrt(u10**2+v10**2)).data.ravel()

wd_selection = 270-(np.arctan2(v10_selection, u10_selection)/np.pi*180).data.ravel()
ws_selection = (np.sqrt(u10_selection**2+v10_selection**2)).data.ravel()

# There are anomalous values in wind direction in Svanvik in 2012 and 2011. Need to filter them out.
# Furthermore if the wind strength is below a certain threashold the direction is set to 0 or 360 deg.
# By checking the probability density function of wind directions
# (probability of 0 or 360 should be about the same as for the directions close-by),
# I set the threshold to 0.8 ms-1.
# ff_threashold = 0.8
# Used brute-force minimazation to get a more precise answer test_wind_threashold
test = test_wind_treashold(data_obs,end=1.1,stepwidth=0.01)
ff_threashold = test[0][np.where(((test[1][0]*test[1][1])**2)==((test[1][0]*test[1][1])**2).min())][0]
data_obs_filtered = data_obs.where(data_obs['FF']>=ff_threashold).dropna()
if selection == 'Svanvik':
    mask = (data_obs_filtered.index>='2012-01') & (data_obs_filtered.index<='2012-02-21')
elif selection == 'Karasjok':
    mask = ((data_obs_filtered.index>='2010-07-13') & (data_obs_filtered.index<='2010-09-01')) | ((data_obs_filtered.index>='2012-10-13') & (data_obs_filtered.index<='2013-02-07'))
    data_obs_filtered = data_obs_filtered.loc[~mask]
 
xlabel = "Met. wind direction (deg)"
ylabel = "Wind strength ($ms^{-1}$)"

fig1 = plt.figure(1)
#hist_hist_diagram(wd, ws, "fennoscandic_wind", xlabel, ylabel, fig1)
hist_hist_diagram(wd_selection, ws_selection, "%s_wind" % selection, xlabel, ylabel, fig1)

fig2 = plt.figure(2)
#windrose(wd, ws, "fennoscandic_windrose")
windrose(wd_selection, ws_selection, "%s_windrose" % selection, fig2)


fig3 = plt.figure(3)
hist_hist_diagram(data_obs_filtered['DD'], data_obs_filtered.dropna()['FF'], "%s_obs_wind" % selection, xlabel, ylabel, fig3, model=False, ff_threashold=ff_threashold) #


fig4 = plt.figure(4)
windrose(data_obs_filtered['DD'], data_obs_filtered['FF'], "%s_obs_windrose" % selection, fig4, model=False)

#fig5 = plt.figure(5)
#fig5.canvas.set_window_title("%s_manually_identify_bad_periodes" % (selection))
#ax51 = plt.subplot()
#data_obs_filtered.dropna()['DD']['2010':'2013'].plot(ax=ax51, ls='None', marker='.')
#data_obs_filtered.where(data_obs_filtered['DD']==236).dropna()['DD']['2010':'2013'].plot(ax=ax51, ls='None', marker='x')
#data_obs_filtered.loc[~mask].dropna()['DD']['2010':'2013'].plot(ax=ax51, ls='None', marker='+')

#for i, iyear in zip(range(5,5+len(glob.glob(src_dir_obs))), data_obs_filtered.index.year.unique()):
#  fig = plt.figure(i)  
#  hist_hist_diagram(data_obs_filtered['DD']['%d' % (iyear)], data_obs_filtered['FF']['%d' % (iyear)], "%s_obs_wind_%d" % (selection, iyear), xlabel, ylabel, fig, model=False)

#fig5 = plt.figure(5, figsize=(16,9))
#ax51 = plt.subplot(211)
#ax52 = plt.subplot(212)

#ax51.plot(test[0], test[1][0], label="bin 0deg")
#ax51.plot(test[0], test[1][1], label="bin 360deg")
#ax52.plot(test[0], , label="bin 360deg")
#ax51.set_xlabel("Windspeed ($ms^{-1}$)")
#ax51.set_ylabel("Residual prob. denisity")
#ax51.legend()
# Show it
plt.show(block=False)

