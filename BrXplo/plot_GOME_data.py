import os, glob
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import * # Get physics constants
import matplotlib.pyplot as plt
import matplotlib.dates as mdates # handling of dates in plotting
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime as dt  # Python standard library datetime module
import sys  
sys.path.append("../mytools/") 
from med_tools import *
#from mytools.netcdf_tools import * # ncdump implementation for python

execfile("read_GOME_data.py")

# Clean up
plt.close('all')
b_hist = False
# Access environment variable for directory
data_dir = os.environ['DATA']
subd = '/BrXplo/GOME_total_BrO/'
src = 'gome_*.nc'
#month = (3,4,5,8,9,10)
month = np.arange(1,13)
month_name = ("January", "February", "March", "April", "May", "June", 
              "July", "August", "September", "October", "November", "December")
try: 
    BrO_data
except NameError:
    print("Reading data...")
    BrO_data = []
    for each in sorted(glob.glob(src)):
        BrO_data.append(xr.open_dataset(each))

    BrO_data = xr.concat(BrO_data, dim='time')
BrO_data_mean = BrO_data.groupby('time.month').mean(dim='time')

# Set map
meridians = np.arange(-180.,181.,20.)
parallels = np.arange(-80.,81.,20.)
m_nh = Basemap(projection='npstere', boundinglat=45., lon_0=0, resolution='l') #, round=True
m_sh = Basemap(projection='spstere', boundinglat=-45., lon_0=0, resolution='l')#, round=True

plt.subplots_adjust(wspace=0.3, hspace=0.1)
fig1 = plt.figure(1, figsize=(11.69,8.27))
fig1.canvas.set_window_title("gome_BrO_2000")

for imonth in range(len(month)):
    BrO_sel = BrO_data_mean.sel(month=month[imonth])#BrO_data.sel(time='2000-04-10')
    # Plot it
    #data_cyclic, lon_cyclic = addcyclic(BrO_data_mean['BrO'][imonth-1], BrO_data_mean['lon']) # 
    data_cyclic, lon_cyclic = addcyclic(BrO_sel['BrO']-3, BrO_sel['lon'])
    # Shift grid so lon runs from -180 to +180
    #data_cyclic, lon_cyclic = shiftgrid(180, data_cyclic, lon_cyclic, start=False)
    # Create "D lat/lon arrays for Basemap
    lon2d, lat2d = np.meshgrid(lon_cyclic, BrO_sel['lat'])
    # Transform lat/lon into plotting coordinates for projection
    x, y = m_nh(lon2d, lat2d)
    ax11 = plt.subplot(4,6,imonth*2+1)
    ax11.set_title(month_name[month[imonth]-1], ha='left', y=0.85, x=0.02, color='white')
    cp11 = ax11.contourf(x, y, data_cyclic, np.arange(0,5.1,0.5), extend='both')
    m_nh.drawcoastlines()
    m_nh.drawmapboundary()
    # draw parallels and meridians.
    m_nh.drawparallels(parallels)
    m_nh.drawmeridians(meridians, labels=[True,False,False,True])
    
    # Transform lat/lon into plotting coordinates for projection
    x, y = m_sh(lon2d, lat2d)
    ax12 = plt.subplot(4,6,imonth*2+2)
    cp12 = ax12.contourf(x, y, data_cyclic, np.arange(0,5.1,0.5), extend='both')
    m_sh.drawcoastlines()
    m_sh.drawmapboundary()
    # draw parallels and meridians.
    m_sh.drawparallels(parallels)
    m_sh.drawmeridians(meridians, labels=[True,False,False,True])
# Set color bar
cp_bar = fig1.colorbar(cp12, ax=fig1.axes, orientation='horizontal', fraction=0.046, pad=0.04, aspect=30) #
cp_bar.set_ticks(np.arange(0,10.1))
cp_bar.set_label("%s (%s %s)" % (BrO_data['BrO'].attrs['long_name'], BrO_data['BrO'].attrs['scaling'], BrO_data['BrO'].attrs['units']))

if b_hist:
    fig2 = plt.figure(2)
    if month < 10:
        fig2.canvas.set_window_title("gome_BrO_2000%s_hist" % (month))
        BrO_hist = (BrO_data_mean['BrO'][month-1]).where(BrO_data_mean.lat>45).plot.hist(bins=100, normed=True, range=(0,8))
    else:
        fig2.canvas.set_window_title("gome_BrO_200%s_hist" % (month))
        BrO_hist = (BrO_data_mean['BrO'][month-1]).where(BrO_data_mean.lat<-45).plot.hist(bins=100, normed=True, range=(0,8))
    ax21 = plt.gca()
    ax21.set_xlabel("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'))
    ax21.set_ylabel("Probability Density")
    ax21.set_ylim(0,1)


plt.show(block=False)
