import os, glob
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import * # Get physics constants
import matplotlib.pyplot as plt
import matplotlib.dates as mdates # handling of dates in plotting
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime as dt  # Python standard library datetime module
from mytools.med_tools import *
from mytools.netcdf_tools import * # ncdump implementation for python
from local_solar_time import local_solar_time

def resample_lst(data):
    '''
    Compute local solar time and resample data accordingly.
    '''
    time = np.arange(-1,23) # shifted time, since 23 is not valid
    lon_shifted = np.concatenate((data.lon.data[np.where(data.lon.data<180)],data.lon.data[np.where(data.lon.data>=180)]-360))
    lst = local_solar_time(time, lon_shifted)
    x1, y1 = np.where((lst<11) & (lst>10))
    lst_lon = data.lon[x1]
    lst_time = time[y1]
    lst_time[np.where(lst_time==-1)] = 23 # "shift" time back
    data_sample = []
    for iday in range(int(len(data.time)/24.)):
        data_test = []
        for i in range(len(lst_time[::5])):
            data_test.append(data.sel(time=data.time[iday*24+lst_time[::5][i]]).where(data.lon==lst_lon[5*i:5*(i+1)]))
        data_test = xr.concat((data_test), dim='lon')
        data_test.coords['time'] = dt.datetime(data_test.coords['time.year'][0].data,data_test.coords['time.month'][0].data,iday+1,10)
        data_sample.append(data_test)
    data_sample = xr.concat((data_sample), dim='time')
    return data_sample
# Clean up
plt.close('all')
nc_src = os.environ['DATA']
subd = '/BrXplo/EMAC_total_BrO/'
#src = 'BrO_col_2000*_BrXplo_mysic_corr.nc'
src = 'BrO_col_2000*_BrXplo_mysic.nc'
time_sel = '2000-05-10'
#src_ref = 'BrO_col_200004_BrXplo_ref.nc'
month = np.arange(1,13)
month_name = ("January", "February", "March", "April", "May", "June", 
              "July", "August", "September", "October", "November", "December")
data_list = []
# Open dataset
for file in sorted(glob.glob(nc_src+subd+src)):
    data = xr.open_dataset(file)
    data_list.append(data)
#data_ref = xr.open_dataset(nc_src+subd+src_ref)
b_check = False
b_hist = False

# Set map
meridians = np.arange(-180.,181.,20.)
parallels = np.arange(-80.,81.,20.)
m_nh = Basemap(projection='npstere', boundinglat=45., lon_0=0, resolution='l') #, round=True
m_sh = Basemap(projection='spstere', boundinglat=-45., lon_0=0, resolution='l')#, round=True

# Plot it
#plot_data = data_sample.sel(time=time_sel)
plt.subplots_adjust(wspace=0.3, hspace=0.1)
fig1 = plt.figure(1, figsize=(12,8))
fig1.canvas.set_window_title("emac_BrO_2000")

for imonth in range(len(month)):
    # Resample data to local solar time
    data_sample = resample_lst(data_list[imonth])
    # Claculate daily mean
    plot_data = data_sample.groupby("time.month").mean(dim='time')
    data_cyclic, lon_cyclic = addcyclic(plot_data['BrO_column'].squeeze(), 
                                        plot_data['lon'].squeeze())
    data_cyclic, lon_cyclic = addcyclic(data_cyclic, lon_cyclic)
    # Shift grid so lon runs from -180 to +180
    #data_cyclic, lon_cyclic = shiftgrid(180, data_cyclic, lon_cyclic, start=False)
    # Create "D lat/lon arrays for Basemap
    lon2d, lat2d = np.meshgrid(lon_cyclic, plot_data['lat'].squeeze())
    # Transform lat/lon into plotting coordinates for projection
    x, y = m_nh(lon2d, lat2d)
    ax11 = plt.subplot(4,6,imonth*2+1)
    ax11.set_title(month_name[month[imonth]-1], ha='left', y=0.85, x=0.02, color='white')
    cp11 = ax11.contourf(x, y, data_cyclic*1e-13, np.arange(0,5.1,0.5), extend='both')
    m_nh.drawcoastlines()
    m_nh.drawmapboundary()
    # draw parallels and meridians.
    m_nh.drawparallels(parallels)
    m_nh.drawmeridians(meridians, labels=[True,False,False,True])
    # Transform lat/lon into plotting coordinates for projection
    x, y = m_sh(lon2d, lat2d)
    ax12 = plt.subplot(4,6,imonth*2+2)
    cp12 = ax12.contourf(x, y, data_cyclic*1e-13, np.arange(0,5.1,0.5), extend='both')
    m_sh.drawcoastlines()
    m_sh.drawmapboundary()
    # draw parallels and meridians.
    m_sh.drawparallels(parallels)
    m_sh.drawmeridians(meridians, labels=[True,False,False,True])
cp_bar = fig1.colorbar(cp12, ax=fig1.axes, orientation='horizontal', fraction=0.046, pad=0.04, aspect=30)
cp_bar.set_ticks(np.arange(0,10.1))
cp_bar.set_label("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'))

if b_hist:
    fig2 = plt.figure(2)
    fig2.canvas.set_window_title('%s_hist' % src[:src.find('.')])
    BrO_hist = (data['BrO_column']*1e-13).where(data.lat>45).plot.hist(bins=100, normed=True, range=(0,8))
    ax21 = plt.gca()
    ax21.set_xlabel("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'))
    ax21.set_ylabel("Probability Density")
    ax21.set_ylim(0,1)
    fig3 = plt.figure(3)
    fig3.canvas.set_window_title('%s_hist_resampled' % src[:src.find('.')])
    (data_sample['BrO_column']*1e-13).where(data.lat>45).plot.hist(bins=100, normed=True, range=(0,8)) #.sel(time='2000-04-01')
    ax31 = plt.gca()
    ax31.set_xlabel("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'))
    ax31.set_ylabel("Probability Density")
    ax31.set_ylim(0,1)

if b_check:
    fig4 = plt.figure(4)
    from matplotlib.colors import BoundaryNorm
    cmap = plt.get_cmap('viridis')
    levels = np.arange(round(lst.min()),round(lst.max()),2)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cm41 = plt.pcolormesh(time,data.lon,lst,norm=norm)
    plt.scatter(lst_time, data.lon[x1], marker='x', color='black')
    ax41 = plt.gca()
    ax41.set_xlabel('Time (hours)')
    ax41.set_xlim(0,24)
    ax41.set_xticks(np.arange(24))
    ax41.set_ylabel('Longitude (deg)')
    ax41.set_ylim(0,360)
    ax41.set_yticks(np.arange(0,370,10))
    cbar = plt.colorbar(cm41)
    cbar.set_label("Local solar time")
# Show it
plt.show(block=False)
