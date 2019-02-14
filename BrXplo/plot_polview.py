import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import cartopy as cp
import cartopy.util as ccrs_util  # Add cyclic
from scipy.constants import * # Get physics constants
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates # handling of dates in plotting
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime as dt  # Python standard library datetime module
sys.path.append("../mytools/") 
from met_tools import *
#from mytools.netcdf_tools import * # ncdump implementation for python
from local_solar_time import local_solar_time


# Clean up
plt.close('all')
# Control
b_anom = True
b_spring = True
b_emac = False

nc_src = os.environ['DATA']
subd = '/BrXplo_data/'
if b_emac:
    src = 'EMAC_total_BrO/mm_BrO_col_2000*_BrXplo_mysic.nc'
else:
    src = 'GOME_total_BrO/gome_bro_*.nc'
    #src = './gome_bro_*.nc'
month = np.arange(1,13)
month_name = ("January", "February", "March", "April", "May", "June", 
              "July", "August", "September", "October", "November", "December")
try:
    data
except NameError:
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src+subd+src)): #
        print(file)
        data = xr.open_dataset(file)
        data_list.append(data)

if b_spring:
    month = (4, 9)

# Plot it
#
if b_spring:
    fig1 = plt.figure(1, figsize=(10,10))
else:
    fig1 = plt.figure(1, figsize=(16,10))
if b_emac:
    if b_anom:
        fig1.canvas.set_window_title("emac_BrO_2000_anom")
    else:
        fig1.canvas.set_window_title("emac_BrO_2000")
else:
    if b_anom:
        fig1.canvas.set_window_title("gome_BrO_2000_anom")
    else:
        fig1.canvas.set_window_title("gome_BrO_2000")
plt.subplots_adjust(wspace=0.3, hspace=0.1)
if b_spring:
    fig1.canvas.set_window_title(fig1.canvas.get_window_title()+"_spring")

for imonth in np.array(month)-1:
    # Set map
    meridians = np.arange(-180.,181.,20.)
    parallels = np.arange(-80.,81.,20.)
    m_nh = Basemap(projection='npstere', boundinglat=45., lon_0=0, resolution='l') #, round=True
    m_sh = Basemap(projection='spstere', boundinglat=-45., lon_0=0, resolution='l')#, round=True
    if b_anom:
        plot_data = data_list[imonth]-data_list[imonth].mean(dim='lon')
    else:
        plot_data = data_list[imonth]
    if b_emac:
        data_cyclic, lon_cyclic = addcyclic(1e-13*plot_data['BrO_column'].data.squeeze(), 
                                            plot_data['lon'].data)
    else:
        data_cyclic, lon_cyclic = addcyclic(plot_data['BrO'].mean(dim='time').data, 
                                            plot_data['lon'].data)
    #data_cyclic, lon_cyclic = ccrs_util.add_cyclic_point(data_cyclic, lon_cyclic)
    # Shift grid so lon runs from -180 to +180
    #data_cyclic, lon_cyclic = shiftgrid(180, data_cyclic, lon_cyclic, start=False)
    # Create "D lat/lon arrays for Basemap
    lon2d, lat2d = np.meshgrid(lon_cyclic, plot_data['lat'])
    # Transform lat/lon into plotting coordinates for projection
    x, y = m_nh(lon2d, lat2d)
    if b_spring:
        ax11 = plt.subplot(2,2,imonth/4+1)
        if b_emac:
            levels = np.arange(0,5.1,0.1)
        else:
            levels = np.arange(0,9.1,0.1)
        
    else:
        ax11 = plt.subplot(4,6,imonth*2+1)
        if b_emac:
            levels = np.arange(0,7.1,0.1)
        else:
            levels = np.arange(0,9.1,0.1)
    if b_anom:
        ax11.set_title(month_name[imonth], ha='left', y=0.85, x=0.02, color='black')
        cp11 = ax11.contourf(x, y, data_cyclic, np.arange(-2,2.1,0.1), extend='both', cmap=plt.cm.RdYlBu_r)
    else:
        ax11.set_title(month_name[imonth], ha='left', y=0.85, x=0.02, color='white')
        cp11 = ax11.contourf(x, y, data_cyclic, levels, extend='max')
    m_nh.drawcoastlines()
    m_nh.drawmapboundary()
    # draw parallels and meridians.
    m_nh.drawparallels(parallels)
    m_nh.drawmeridians(meridians, labels=[True,False,False,True])
    # Transform lat/lon into plotting coordinates for projection
    x, y = m_sh(lon2d, lat2d)
    if b_spring:
        ax12 = plt.subplot(2,2,imonth/4+2)
    else:
        ax12 = plt.subplot(4,6,imonth*2+2)
    if b_anom:
        ax12.set_title(month_name[imonth], ha='left', y=0.85, x=0.02, color='black')
        cp12 = ax12.contourf(x, y, data_cyclic, np.arange(-2,2.1,0.1), extend='both', cmap=plt.cm.RdYlBu_r)
    else:
        ax12.set_title(month_name[imonth], ha='left', y=0.85, x=0.02, color='white')
        cp12 = ax12.contourf(x, y, data_cyclic, levels, extend='max')
    m_sh.drawcoastlines()
    m_sh.drawmapboundary()
    # draw parallels and meridians.
    m_sh.drawparallels(parallels)
    m_sh.drawmeridians(meridians, labels=[True,False,False,True])
cp_bar = fig1.colorbar(cp12, ax=fig1.axes, orientation='horizontal', fraction=0.046, pad=0.04, aspect=30)
if b_anom:
    cp_bar.set_ticks(np.arange(-2,2.1))
else:
    cp_bar.set_ticks(np.arange(0,10.1))
cp_bar.set_label("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'))

# Show it
plt.show(block=False)
#plt.savfig(fig1.canvas.get_window_title()+".png", format='png')
