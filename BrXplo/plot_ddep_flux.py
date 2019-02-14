import os, glob
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import * # Get physics constants
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime as dt  # Python standard library datetime module
from mytools.med_tools import print_all, seconds_in_month, find_nearest
from mytools.netcdf_tools import average


plt.close('all')
# Read data
nc_src = os.environ['DATA']
nc_subd = "/BrXplo/BrXplo_ref/"
nc_ref = "/BrXplo_ref_____2000*ddep_gp.nc"

try: 
    ddepflux_HBr
except NameError:
    print("Reading BrXplo data...")
    ddepflux_HBr = []
    ddepflux_HOBr = []
    ddepflux_BrNO3 = []
    ddepflux_BrCl = []
    ddepflux_CH3Br = []
    ddepflux_CHBr3 = []
    ddepflux_CH2Br2 = []
    ddepflux_O3 = []
    for file in sorted(glob.glob(nc_src+nc_subd+nc_ref)):
        data = xr.open_dataset(file)
        ddepflux_HBr.append(data['ddepflux_HBr'])
        ddepflux_HOBr.append(data['ddepflux_HOBr'])
        ddepflux_BrNO3.append(data['ddepflux_BrNO3'])
        ddepflux_BrCl.append(data['ddepflux_BrCl'])
        ddepflux_CH3Br.append(data['ddepflux_CH3Br'])
        ddepflux_CHBr3.append(data['ddepflux_CHBr3'])
        ddepflux_CH2Br2.append(data['ddepflux_CH2Br2'])
        ddepflux_O3.append(data['ddepflux_O3'])
    ddepflux_HBr = xr.concat(ddepflux_HBr, dim='time')
    ddepflux_HOBr = xr.concat(ddepflux_HOBr, dim='time')
    ddepflux_BrNO3 = xr.concat(ddepflux_BrNO3, dim='time')
    ddepflux_BrCl = xr.concat(ddepflux_BrCl, dim='time')
    ddepflux_CH3Br = xr.concat(ddepflux_CH3Br, dim='time')
    ddepflux_CHBr3 = xr.concat(ddepflux_CHBr3, dim='time')
    ddepflux_CH2Br2 = xr.concat(ddepflux_CH2Br2, dim='time')
    ddepflux_O3 = xr.concat(ddepflux_O3, dim='time')

# Read ESCiMo data
nc_subd = "/ESCIMO/RC1-base-08/"
nc_ref = "/ddep_gp/RC1-base-08____2000*ddep_gp.nc"
try: 
    rc1_ddepflux_HBr
except NameError:
    print("Reading RC1-base-08 data...")
    rc1_ddepflux_HBr = []
    rc1_ddepflux_HOBr = []
    rc1_ddepflux_BrNO3 = []
    rc1_ddepflux_BrCl = []
    rc1_ddepflux_CH3Br = []
    rc1_ddepflux_CHBr3 = []
    rc1_ddepflux_CH2Br2 = []
    rc1_ddepflux_O3 = []
    for file in sorted(glob.glob(nc_src+nc_subd+nc_ref)):
        data = xr.open_dataset(file)
        rc1_ddepflux_HBr.append(data['ddepflux_HBr'])
        rc1_ddepflux_HOBr.append(data['ddepflux_HOBr'])
        rc1_ddepflux_BrNO3.append(data['ddepflux_BrNO3'])
        rc1_ddepflux_BrCl.append(data['ddepflux_BrCl'])
        rc1_ddepflux_CH3Br.append(data['ddepflux_CH3Br'])
        rc1_ddepflux_CHBr3.append(data['ddepflux_CHBr3'])
        rc1_ddepflux_CH2Br2.append(data['ddepflux_CH2Br2'])
        rc1_ddepflux_O3.append(data['ddepflux_O3'])
    rc1_ddepflux_HBr = xr.concat(rc1_ddepflux_HBr, dim='time')
    rc1_ddepflux_HOBr = xr.concat(rc1_ddepflux_HOBr, dim='time')
    rc1_ddepflux_BrNO3 = xr.concat(rc1_ddepflux_BrNO3, dim='time')
    rc1_ddepflux_BrCl = xr.concat(rc1_ddepflux_BrCl, dim='time')
    rc1_ddepflux_CH3Br = xr.concat(rc1_ddepflux_CH3Br, dim='time')
    rc1_ddepflux_CHBr3 = xr.concat(rc1_ddepflux_CHBr3, dim='time')
    rc1_ddepflux_CH2Br2 = xr.concat(rc1_ddepflux_CH2Br2, dim='time')
    rc1_ddepflux_O3 = xr.concat(rc1_ddepflux_O3, dim='time')

nc_subd = "/ESCIMO/RC2-base-05/"
nc_ref = "/ddep_gp/RC2-base-05____2000*ddep_gp.nc"
try: 
    rc2_ddepflux_HBr
except NameError:
    print("Reading RC2-base-05 data...")
    rc2_ddepflux_HBr = []
    rc2_ddepflux_HOBr = []
    rc2_ddepflux_BrNO3 = []
    rc2_ddepflux_BrCl = []
    rc2_ddepflux_CH3Br = []
    rc2_ddepflux_CHBr3 = []
    rc2_ddepflux_CH2Br2 = []
    rc2_ddepflux_O3 = []
    for file in sorted(glob.glob(nc_src+nc_subd+nc_ref)):
        data = xr.open_dataset(file)
        rc2_ddepflux_HBr.append(data['ddepflux_HBr'])
        rc2_ddepflux_HOBr.append(data['ddepflux_HOBr'])
        rc2_ddepflux_BrNO3.append(data['ddepflux_BrNO3'])
        rc2_ddepflux_BrCl.append(data['ddepflux_BrCl'])
        rc2_ddepflux_CH3Br.append(data['ddepflux_CH3Br'])
        rc2_ddepflux_CHBr3.append(data['ddepflux_CHBr3'])
        rc2_ddepflux_CH2Br2.append(data['ddepflux_CH2Br2'])
        rc2_ddepflux_O3.append(data['ddepflux_O3'])
    rc2_ddepflux_HBr = xr.concat(rc2_ddepflux_HBr, dim='time')
    rc2_ddepflux_HOBr = xr.concat(rc2_ddepflux_HOBr, dim='time')
    rc2_ddepflux_BrNO3 = xr.concat(rc2_ddepflux_BrNO3, dim='time')
    rc2_ddepflux_BrCl = xr.concat(rc2_ddepflux_BrCl, dim='time')
    rc2_ddepflux_CH3Br = xr.concat(rc2_ddepflux_CH3Br, dim='time')
    rc2_ddepflux_CHBr3 = xr.concat(rc2_ddepflux_CHBr3, dim='time')
    rc2_ddepflux_CH2Br2 = xr.concat(rc2_ddepflux_CH2Br2, dim='time')
    rc2_ddepflux_O3 = xr.concat(rc2_ddepflux_O3, dim='time')

# Computations
sec_in_year = np.array([seconds_in_month(each,2000) for each in np.arange(1,13)]).sum()
ddep_range = np.arange(0, 2, 0.1)
# Set map
meridians = np.arange(-180.,181.,20.)
parallels = np.arange(-80.,81.,20.)
m_nh = Basemap(projection='npstere', boundinglat=45., lon_0=0, resolution='l') #, round=True
m_sh = Basemap(projection='spstere', boundinglat=-45., lon_0=0, resolution='l')#, round=True

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("integrated_ddep_Bry")
ddepflux_Bry = (((ddepflux_HBr+ddepflux_HOBr+ddepflux_BrNO3+ddepflux_BrCl)*60**2).groupby("time.year").sum(dim='time'))
data_cyclic, lon_cyclic = addcyclic(ddepflux_Bry.data.squeeze(), ddepflux_Bry['lon'].squeeze())
lon2d, lat2d = np.meshgrid(lon_cyclic, ddepflux_Bry['lat'].squeeze())
# Transform lat/lon into plotting coordinates for projection
x, y = m_nh(lon2d, lat2d)
ax11 = plt.subplot(121)
ax11.set_title("Br$_y$")
cp11 = ax11.contourf(x, y, np.log10(data_cyclic*1e-17), ddep_range, extend='both')
m_nh.drawcoastlines()
m_nh.drawmapboundary()
# draw parallels and meridians.
m_nh.drawparallels(parallels)
m_nh.drawmeridians(meridians, labels=[True,False,False,True])
# Transform lat/lon into plotting coordinates for projection
x, y = m_sh(lon2d, lat2d)
ax12 = plt.subplot(122)
cp12 = ax12.contourf(x, y, np.log10(data_cyclic*1e-17), ddep_range, extend='both')
m_sh.drawcoastlines()
m_sh.drawmapboundary()
# draw parallels and meridians.
m_sh.drawparallels(parallels)
m_sh.drawmeridians(meridians, labels=[True,False,False,True])
cp_bar = fig1.colorbar(cp11, ax=fig1.axes, orientation='horizontal', fraction=0.046, pad=0.04, aspect=30)
cp_bar.set_label("%s (%s)" % ("log$_{10}$ dry deposition flux", "10$^{17}$ molec m$^{-2}$"))

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("integrated_ddep_BrOrg")
ddepflux_BrOrg = (((ddepflux_CH3Br+ddepflux_CHBr3+ddepflux_CH2Br2)*60**2).groupby("time.year").sum(dim='time'))
data_cyclic, lon_cyclic = addcyclic(ddepflux_BrOrg.data.squeeze(), ddepflux_BrOrg['lon'].squeeze())
lon2d, lat2d = np.meshgrid(lon_cyclic, ddepflux_BrOrg['lat'].squeeze())
# Transform lat/lon into plotting coordinates for projection
x, y = m_nh(lon2d, lat2d)
ax21 = plt.subplot(121)
ax21.set_title("Br$_{org}$")
cp21 = ax21.contourf(x, y, np.log10(data_cyclic*1e-16), ddep_range, extend='both')
m_nh.drawcoastlines()
m_nh.drawmapboundary()
# draw parallels and meridians.
m_nh.drawparallels(parallels)
m_nh.drawmeridians(meridians, labels=[True,False,False,True])
# Transform lat/lon into plotting coordinates for projection
x, y = m_sh(lon2d, lat2d)
ax22 = plt.subplot(122)
cp22 = ax22.contourf(x, y, np.log10(data_cyclic*1e-16), ddep_range, extend='both')
m_sh.drawcoastlines()
m_sh.drawmapboundary()
# draw parallels and meridians.
m_sh.drawparallels(parallels)
m_sh.drawmeridians(meridians, labels=[True,False,False,True])
cp_bar = fig2.colorbar(cp21, ax=fig2.axes, orientation='horizontal', fraction=0.046, pad=0.04, aspect=30)
cp_bar.set_label("%s (%s)" % ("log$_{10}$ dry deposition flux", "10$^{16}$ molec m$^{-2}$"))

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("integrated_ddep_ozone")
ddepflux_ozone = ((ddepflux_O3*60**2).groupby("time.year").sum(dim='time'))
data_cyclic, lon_cyclic = addcyclic(ddepflux_ozone.data.squeeze(), ddepflux_ozone['lon'].squeeze())
lon2d, lat2d = np.meshgrid(lon_cyclic, ddepflux_ozone['lat'].squeeze())
# Transform lat/lon into plotting coordinates for projection
x, y = m_nh(lon2d, lat2d)
ax31 = plt.subplot(121)
ax31.set_title("Ozone")
cp31 = ax31.contourf(x, y, np.log10(data_cyclic*1e-21), ddep_range, extend='both')
m_nh.drawcoastlines()
m_nh.drawmapboundary()
# draw parallels and meridians.
m_nh.drawparallels(parallels)
m_nh.drawmeridians(meridians, labels=[True,False,False,True])
# Transform lat/lon into plotting coordinates for projection
x, y = m_sh(lon2d, lat2d)
ax32 = plt.subplot(122)
cp32 = ax32.contourf(x, y, np.log10(data_cyclic*1e-21), ddep_range, extend='both')
m_sh.drawcoastlines()
m_sh.drawmapboundary()
# draw parallels and meridians.
m_sh.drawparallels(parallels)
m_sh.drawmeridians(meridians, labels=[True,False,False,True])
cp_bar = fig3.colorbar(cp31, ax=fig3.axes, orientation='horizontal', fraction=0.046, pad=0.04, aspect=30)
cp_bar.set_label("%s (%s)" % ("log$_{10}$ dry deposition flux", "10$^{21}$ molec m$^{-2}$"))

fig4 = plt.figure(4)
fig4.canvas.set_window_title("ddep_br_timeseries_nh")
# Shift ESCiMo times to middle of month
times = pd.date_range('1999-12-01', '2000-11-30', name='time', freq='1M').shift(15, freq=pd.datetools.day)
ddepflux_BrOrg_mean = average((ddepflux_CH3Br+ddepflux_CHBr3+ddepflux_CH2Br2).mean(dim='lon').where(ddepflux_CH3Br.lat>50), weights=np.cos(ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
ddepflux_Bry_mean = average((ddepflux_HBr+ddepflux_HOBr+ddepflux_BrNO3+ddepflux_BrCl).mean(dim='lon').where(ddepflux_CH3Br.lat>50), weights=np.cos(ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
#ddepflux_O3_mean = average((ddepflux_O3).mean(dim='lon').where(ddepflux_CH3Br.lat>50), weights=np.cos(ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-15
#ddepflux_Bry_alert = (ddepflux_HBr+ddepflux_HOBr+ddepflux_BrNO3+ddepflux_BrCl)[:,find_nearest(82.5, ddepflux_HBr['lat'], index=True).data,find_nearest(360-62.3, ddepflux_HBr['lon'].data, index=True)]*1e-10
ddepflux_Bry_mean.plot(label='Br$_y$')
ddepflux_BrOrg_mean.plot(label='Br$_{org}$')
#ddepflux_Bry_alert.plot()
#ddepflux_O3_mean.plot(label='O$_3$', color='dimgrey')
rc1_ddepflux_BrOrg_mean = average((rc1_ddepflux_CH3Br+rc1_ddepflux_CHBr3+rc1_ddepflux_CH2Br2).mean(dim='lon').where(rc1_ddepflux_CH3Br.lat>50), weights=np.cos(rc1_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
rc1_ddepflux_Bry_mean = average((rc1_ddepflux_HBr+rc1_ddepflux_HOBr+rc1_ddepflux_BrNO3+rc1_ddepflux_BrCl).mean(dim='lon').where(rc1_ddepflux_CH3Br.lat>50), weights=np.cos(rc1_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
rc1_ddepflux_O3_mean = average((rc1_ddepflux_O3).mean(dim='lon').where(rc1_ddepflux_CH3Br.lat>50), weights=np.cos(rc1_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-15
rc1_ddepflux_BrOrg_mean.coords['time'] = times
rc1_ddepflux_Bry_mean.coords['time'] = times
#rc1_ddepflux_O3_mean.coords['time'] = times
rc1_ddepflux_Bry_mean.plot(label='Br$_y$ (RC1-base-08)', color='cornflowerblue', ls='--')
rc1_ddepflux_BrOrg_mean.plot(label='Br$_{org}$ (RC1-base-08)', color='tomato', ls='--')
#rc1_ddepflux_O3_mean.plot(label='O$_3$ (RC1-base-08)', color='black', ls='--')
rc2_ddepflux_BrOrg_mean = average((rc2_ddepflux_CH3Br+rc2_ddepflux_CHBr3+rc2_ddepflux_CH2Br2).mean(dim='lon').where(rc2_ddepflux_CH3Br.lat>50), weights=np.cos(rc2_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
rc2_ddepflux_Bry_mean = average((rc2_ddepflux_HBr+rc2_ddepflux_HOBr+rc2_ddepflux_BrNO3+rc2_ddepflux_BrCl).mean(dim='lon').where(rc2_ddepflux_CH3Br.lat>50), weights=np.cos(rc2_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
rc2_ddepflux_O3_mean = average((rc2_ddepflux_O3).mean(dim='lon').where(rc2_ddepflux_CH3Br.lat>50), weights=np.cos(rc2_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-15
rc2_ddepflux_BrOrg_mean.coords['time'] = times
rc2_ddepflux_Bry_mean.coords['time'] = times
#rc2_ddepflux_O3_mean.coords['time'] = times
rc2_ddepflux_Bry_mean.plot(label='Br$_y$ (RC2-base-05)', color='darkblue', ls='-.')
rc2_ddepflux_BrOrg_mean.plot(label='Br$_{org}$ (RC2-base-05)', color='darkred', ls='-.')
#rc2_ddepflux_O3_mean.plot(label='O$_3$ (RC2-base-05)', color='black', ls='-.')
ax41 = plt.gca()
ax41.set_title("$50^\circ - 90^\circ N$")
ax41.set_ylabel("%s (%s)" % ("Dry deposition flux", "10$^{10}$ molec m$^{-2}$ s$^{-1}$"))
ax41.set_ylim(0,10)
ax41.legend(frameon=False, loc='best', ncol=3)
#plt.text(dt.date(2000, 12, 25), 0.2,"x10$^{5}$")

fig5 = plt.figure(5)
fig5.canvas.set_window_title("ddep_br_timeseries_sh")
# Shift ESCiMo times to middle of month
times = pd.date_range('1999-12-01', '2000-11-30', name='time', freq='1M').shift(15, freq=pd.datetools.day)
ddepflux_BrOrg_mean = average((ddepflux_CH3Br+ddepflux_CHBr3+ddepflux_CH2Br2).mean(dim='lon').where(ddepflux_CH3Br.lat<-50), weights=np.cos(ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
ddepflux_Bry_mean = average((ddepflux_HBr+ddepflux_HOBr+ddepflux_BrNO3+ddepflux_BrCl).mean(dim='lon').where(ddepflux_CH3Br.lat<-50), weights=np.cos(ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
ddepflux_Bry_mean.plot(label='Br$_y$')
ddepflux_BrOrg_mean.plot(label='Br$_{org}$')
rc1_ddepflux_BrOrg_mean = average((rc1_ddepflux_CH3Br+rc1_ddepflux_CHBr3+rc1_ddepflux_CH2Br2).mean(dim='lon').where(rc1_ddepflux_CH3Br.lat<-50), weights=np.cos(rc1_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
rc1_ddepflux_Bry_mean = average((rc1_ddepflux_HBr+rc1_ddepflux_HOBr+rc1_ddepflux_BrNO3+rc1_ddepflux_BrCl).mean(dim='lon').where(rc1_ddepflux_CH3Br.lat<-50), weights=np.cos(rc1_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
rc1_ddepflux_BrOrg_mean.coords['time'] = times
rc1_ddepflux_Bry_mean.coords['time'] = times
rc1_ddepflux_Bry_mean.plot(label='Br$_y$ (RC1-base-08)', color='cornflowerblue', ls='--')
rc1_ddepflux_BrOrg_mean.plot(label='Br$_{org}$ (RC1-base-08)', color='tomato', ls='--')
rc2_ddepflux_BrOrg_mean = average((rc2_ddepflux_CH3Br+rc2_ddepflux_CHBr3+rc2_ddepflux_CH2Br2).mean(dim='lon').where(rc2_ddepflux_CH3Br.lat<-50), weights=np.cos(rc2_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
rc2_ddepflux_Bry_mean = average((rc2_ddepflux_HBr+rc2_ddepflux_HOBr+rc2_ddepflux_BrNO3+rc2_ddepflux_BrCl).mean(dim='lon').where(rc2_ddepflux_CH3Br.lat<-50), weights=np.cos(rc2_ddepflux_CH3Br.lat*np.pi/180), dim='lat')*1e-10
rc2_ddepflux_BrOrg_mean.coords['time'] = times
rc2_ddepflux_Bry_mean.coords['time'] = times
rc2_ddepflux_Bry_mean.plot(label='Br$_y$ (RC2-base-05)', color='darkblue', ls='-.')
rc2_ddepflux_BrOrg_mean.plot(label='Br$_{org}$ (RC2-base-05)', color='darkred', ls='-.')
ax51 = plt.gca()
ax51.set_title("$50^\circ - 90^\circ S$")
ax51.set_ylabel("%s (%s)" % ("Dry deposition flux", "10$^{10}$ molec m$^{-2}$ s$^{-1}$"))
ax51.set_ylim(0,2)
ax51.legend(frameon=False, loc='best', ncol=3)
# Show it
plt.show(block=False)
