import os, glob
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import * # Get physics constants
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from mytools.med_tools import print_all

nc_src = os.environ['DATA']
src = 'sic_multi-year_1980_2013.nc'
nc_subd = '/ESCIMO'
nc_src_ref = '/RC1SD-base-08/ECHAM5/RC1SD-base-08__200004_ECHAM5.nc'
nc_src_ref_sh = '/RC1SD-base-08/ECHAM5/RC1SD-base-08__200009_ECHAM5.nc'
# Open dataset
data = xr.open_dataset(nc_src+'/'+src)
data_nh = xr.open_dataset(nc_src+nc_subd+nc_src_ref) 
data_sh = xr.open_dataset(nc_src+nc_subd+nc_src_ref_sh)

sic_nh = data_nh.icecov.mean(dim='time')#sel(time='2000-04-01')[0]
sic_sh = data_sh.icecov.mean(dim='time')#sel(time='2000-09-01')[0]

# Set map
meridians = np.arange(-180.,181.,20.)
parallels = np.arange(-80.,81.,20.)
m_nh = Basemap(projection='npstere', boundinglat=45., lon_0=0, resolution='l') 
m_sh = Basemap(projection='spstere', boundinglat=-45., lon_0=0, resolution='l')

# Clean up
plt.close('all')
# Plot it
fig1 = plt.figure(1, figsize=(12,8))
fig1.canvas.set_window_title(src[:src.find('.')])
# Claculate daily mean
data_cyclic, lon_cyclic = addcyclic(data['my_sic'].sel(time=2000), data['lon'])
# Shift grid so lon runs from -180 to +180
data_cyclic, lon_cyclic = shiftgrid(180, data_cyclic, lon_cyclic, start=False)
# Create "D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lon_cyclic, data['lat'])
# Transform lat/lon into plotting coordinates for projection
x, y = m_nh(lon2d, lat2d)
ax11 = plt.subplot(121)
cp11 = ax11.contourf(x, y, data_cyclic, np.arange(0,1.1,0.1), cmap=plt.cm.CMRmap_r, extend='both')
m_nh.drawcoastlines()
m_nh.drawmapboundary()
# draw parallels and meridians.
m_nh.drawparallels(parallels)
m_nh.drawmeridians(meridians, labels=[True,False,False,True])
# Transform lat/lon into plotting coordinates for projection
x, y = m_sh(lon2d, lat2d)
ax12 = plt.subplot(122)
cp12 = ax12.contourf(x, y, data_cyclic, np.arange(0,1.1,0.1), cmap=plt.cm.CMRmap_r, extend='both')
m_sh.drawcoastlines()
m_sh.drawmapboundary()
# draw parallels and meridians.
m_sh.drawparallels(parallels)
m_sh.drawmeridians(meridians, labels=[True,False,False,True])
# Add color bar
cp_bar = fig1.colorbar(cp12, ax=(ax11,ax12), orientation='horizontal', aspect=30)
cp_bar.set_ticks(np.arange(0,1.1,0.1))
cp_bar.set_label("%s" % ('Multi-year Sea Ice Cover'))


fig2 = plt.figure(2, figsize=(12,8))
# Set color of hatches
mpl.rcParams['hatch.color'] = 'c'
fig2.canvas.set_window_title("sic_mysic_200004_09")
# Claculate daily mean
data_cyclic_nh, lon_cyclic = addcyclic(sic_nh, data_nh['lon']) 
# Shift grid so lon runs from -180 to +180
data_cyclic_nh, lon_cyclic = shiftgrid(180, data_cyclic_nh, lon_cyclic, start=False)
# Create "D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lon_cyclic, data['lat'])
# Transform lat/lon into plotting coordinates for projection
x, y = m_nh(lon2d, lat2d)
ax21 = plt.subplot(121)
cp21 = ax21.contourf(x, y, data_cyclic_nh, np.arange(0,1.1,0.1), cmap=plt.cm.CMRmap_r, extend='both')
cp23 = ax21.contourf(x, y, data_cyclic, 
                     levels=np.arange(0,1.1,0.2), hatches=[ None, '....', '////', '\\\\','----'], 
                     colors='none')
m_nh.drawcoastlines()
m_nh.drawmapboundary()
# draw parallels and meridians.
m_nh.drawparallels(parallels)
m_nh.drawmeridians(meridians, labels=[True,False,False,True])
# Claculate daily mean
data_cyclic_sh, lon_cyclic = addcyclic(sic_sh, data_sh['lon']) 
# Shift grid so lon runs from -180 to +180
data_cyclic_sh, lon_cyclic = shiftgrid(180, data_cyclic_sh, lon_cyclic, start=False)
# Create "D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lon_cyclic, data['lat'])
# Transform lat/lon into plotting coordinates for projection
x, y = m_sh(lon2d, lat2d)
ax22 = plt.subplot(122)
cp22 = ax22.contourf(x, y, data_cyclic_sh, np.arange(0,1.1,0.1), cmap=plt.cm.CMRmap_r, extend='both')
cp24 = ax22.contourf(x, y, data_cyclic, 
                     levels=np.arange(0,1.1,0.2), hatches=[ None, '....', '////', '\\\\','----'], 
                     colors='none')
m_sh.drawcoastlines()
m_sh.drawmapboundary()
# draw parallels and meridians.
m_sh.drawparallels(parallels)
m_sh.drawmeridians(meridians, labels=[True,False,False,True])
# Add hatched color bar
cp_bar2 = fig2.colorbar(cp24, ax=(ax21,ax22), orientation='vertical', aspect=30, shrink=0.6)
cp_bar2.set_ticks(np.arange(0,1.1,0.2))
cp_bar2.set_label("%s" % ('Multi-year Sea Ice Cover'))
# Add color bar
cp_bar = fig2.colorbar(cp22, ax=(ax21,ax22), orientation='horizontal', aspect=30)
cp_bar.set_ticks(np.arange(0,1.1,0.1))
cp_bar.set_label("%s" % ('Sea Ice Cover'))


# Show it
plt.show(block=False)
