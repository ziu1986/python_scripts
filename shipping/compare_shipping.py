import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
#import nc_time_axis
#import cftime
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Clean-up
plt.close('all')

species = 'SO2'
label = 'SO_2'

# Data source
nc_src_historic = os.environ['DATA']+'/CTM3_input_data/EMIS/CEDS_CICERO/'+species+'*200001-201412*.nc'
nc_src_ssp585 = os.environ['DATA']+'/astra_data/input_data/ctm_input/EMIS/CEDS_SCENARIOS/SSP585/'+species+'*201501-210012.nc'

try:
    data
except NameError:
    for file in sorted(glob.glob(nc_src_historic)):
        print(file)
        data_historic = xr.open_dataset(file)
    for file in sorted(glob.glob(nc_src_ssp585)):
        print(file)
        data_ssp585= xr.open_dataset(file) 

    ssp585_global_sum = []
    ssp585 = []
    time_ssp585 = []
    for iyear in np.unique(data_ssp585['time'].dt.year):
        ssp585_global_sum.append((data_ssp585['%s_em_anthro_SHI'%(species)].sel(time='%d' % (iyear))*
                                  data_historic['Gridcell area']).sum(dim=('lat','lon'))*
                                 [seconds_in_month(imonth,iyear) for imonth in range(1,13)])
        ssp585.append(xr.concat([data_ssp585['%s_em_anthro_SHI'%(species)].sel(time='%d-%s' % (iyear,str(imonth).zfill(2)))*seconds_in_month(imonth,iyear) for imonth in range(1,13)],dim='time'))
        time_ssp585.append([dt.datetime(iyear,imonth,15) for imonth in range(1,13)])
    ssp585_global_sum = xr.concat(ssp585_global_sum,dim='time').drop('sector')
    ssp585 = xr.concat(ssp585,dim='time').drop('sector')              

    historic_global_sum = []
    historic = []                
    time_historic = []
    for iyear in np.unique(data_historic['time'].dt.year):
        historic_global_sum.append((data_historic['%s_em_anthro_SHI'%(species)].sel(time='%d' % (iyear))*
                                    data_historic['Gridcell area']).sum(dim=('lat','lon'))*
                                   [seconds_in_month(imonth,iyear) for imonth in range(1,13)])
        historic.append(xr.concat([data_historic['%s_em_anthro_SHI'%(species)].sel(time='%d-%s' % (iyear,str(imonth).zfill(2)))*seconds_in_month(imonth,iyear) for imonth in range(1,13)],dim='time'))
        time_historic.append([dt.datetime(iyear,imonth,15) for imonth in range(1,13)])
    historic_global_sum = xr.concat(historic_global_sum,dim='time')
    historic = xr.concat(historic,dim='time')

ssp585_nordic_sea = (ssp585*data_historic['Gridcell area']).sel(lat=slice(60,90)).sum(dim=('lat','lon'))
historic_nordic_sea = (historic*data_historic['Gridcell area']).sel(lat=slice(60,90)).sum(dim=('lat','lon'))

#Plotting
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("CEDS_shipping_emissions_time_%s" % (species))
ax11 = plt.subplot()

ax11.plot(np.array(time_ssp585).flatten(), ssp585_global_sum*1e-6, color='blue', alpha=0.5, label='global: SSP585')
ax11.plot(np.array(time_historic).flatten(), historic_global_sum*1e-6, color='blue', label='global: historical')

ax11.plot(np.array(time_ssp585).flatten(), ssp585_nordic_sea*1e-6, ls='--', color='blue', alpha=0.5, label='North Atlantic: SSP585')
ax11.plot(np.array(time_historic).flatten(), historic_nordic_sea*1e-6, ls='--', color='blue', label='North Atlantic: historical')

ax11.set_xlabel("Time (year)")
ax11.set_ylabel("$%s$ (Gg)" %(label))

ax11.legend()

# Show it
plt.show(block=False)

fig2 = plt.figure(2, figsize=(16,6))
fig2.canvas.set_window_title("CEDS_shipping_emissions_polar_%s" % (species))
ax21 = plt.subplot(131, projection=cp.crs.NorthPolarStereo())
ax22 = plt.subplot(132, projection=cp.crs.NorthPolarStereo())
ax23 = plt.subplot(133, projection=cp.crs.NorthPolarStereo())


# Set the map and axis attributes
for ax in fig2.axes:
    ax.set_global()   # Expands the map to fit the tick labels!
    ax.coastlines()
    #ax.gridlines()

    # Limit the map to 60 degrees latitude and above.
    ax.set_extent([-180, 180, 90, 60], cp.crs.PlateCarree())
    # Adding lines
    draw_parallels(ax, np.arange(60,91,10))
    draw_meridians(ax, np.arange(-180,181,45))

levels = np.arange(0,(ssp585.sel(time='2050').sum(dim='time')).max().round(4),(ssp585.sel(time='2050').sum(dim='time')).max().round(4)/10.)
(historic.sel(time='2005').sum(dim='time')).plot(ax=ax21, transform=cp.crs.PlateCarree(), levels=levels, cmap=plt.cm.hot_r, cbar_kwargs={'label':'$%s$ (kg m$^{-2}$)'%(label),'fraction':0.046, 'pad':0.04,'aspect':30, 'orientation':'horizontal'})
(ssp585.sel(time='2050').sum(dim='time')).plot(ax=ax22, transform=cp.crs.PlateCarree(), levels=levels, cmap=plt.cm.hot_r, cbar_kwargs={'label':'$%s$ (kg m$^{-2}$)'%(label),'fraction':0.046, 'pad':0.04,'aspect':30, 'orientation':'horizontal'})
((1-historic.sel(time='2005').sum(dim='time')/ssp585.sel(time='2050').sum(dim='time'))*100).plot(ax=ax23, transform=cp.crs.PlateCarree(), levels=np.arange(-100,100.1,1), cmap=plt.cm.seismic, cbar_kwargs={'label':'$\Delta %s (%s)$'%(label, '\%'),'fraction':0.046, 'pad':0.04,'aspect':30, 'orientation':'horizontal'})

ax21.set_title('2005')
ax22.set_title('2050')
ax23.set_title('2050-2005')



