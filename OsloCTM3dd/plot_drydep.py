import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *

R_earth = 6371*kilo #m

# Data source
nc_src = os.environ['DATA']+'/abel/C3RUN_nitrate_oDD_2005/scavenging_monthly/sum_scavenging*.nc'
nc_src_ndd = os.environ['DATA']+'/abel/C3RUN_nitrate_nDD_2005/mm_scavenging*.nc'
nc_src_mstc = os.environ['DATA']+'/abel/C3RUN_nitrate_mstc_2005/scavenging_monthly/sum_scavenging*.nc'

try:
    data
except NameError:
    data_list = []
    raw_data = []
    for subdir in (nc_src,nc_src_mstc): #, nc_src_ndd, nc_src_mstc
        print("Reading from path %s" % (os.path.abspath(subdir)))
        # Open dataset
        for file in sorted(glob.glob(subdir)):
            print("Reading %s" % (os.path.basename(file)))
            data = xr.open_dataset(file)
            # Defining new time coordinates
            data['time'].reset_coords(drop=True)
            data.coords['time'] = ([dt.datetime(data['YEAR'], data['MONTH'], 15),])
            #data = data.swap_dims({'time':'date'})
            data_list.append(data)
        # Concatinating the list
        raw_data.append(xr.concat(data_list, dim='time'))
# Extract ozone drydeposition and rescale units
ozone_data = [data['dry_O3']/mega for data in raw_data]
for data in ozone_data:
    data.attrs['unit'] = 'Gg'
# Monthly zonal total and rescale units
ozone_zonal = [data.sum(dim='lon')/kilo for data in ozone_data] #sum(dim='time').
for data in ozone_zonal:
    data.attrs['unit'] = 'Tg'
# Extracting some general information
gridarea = raw_data[0]['gridarea'].isel(time=0)
tot_surface_area = gridarea.sum()
gridarea_faction = gridarea/tot_surface_area
molarweight = get_molarweight(raw_data[0].isel(time=0))

#Plotting
# Clean-up
plt.close('all')

ozone_data_cyclic_da = [addcyclicpoint(data, data['lon']) for data in ozone_data]
for data in ozone_data_cyclic_da:
    data.attrs['units'] = ozone_data[0].attrs['unit']
ozone_dd_cyclic = [addcyclicpoint(data/gridarea, data['lon']) for data in ozone_data]
for data in ozone_dd_cyclic:
    data.attrs['units'] = ozone_data[0].attrs['unit']+'/m2'

# Set the map and axis attributes
def plot_global_ozone_drydep(ax, data, **kargs):
    mode = kargs.pop('stats', 'sum')
    title = kargs.pop('title',"")
    ax.set_global()   # Expands the map to fit the tick labels!
    ax.coastlines()
    ax.set_xticks(np.arange(-180, 181, 45), crs=cp.crs.PlateCarree())
    ax.set_yticks(np.arange(-80, 81, 20), crs=cp.crs.PlateCarree())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    if mode=='sum':
        data.sum(dim='time').plot(ax=ax,
                                  levels=np.arange(0,151, 25),
                                  transform=cp.crs.PlateCarree(),
                                  cbar_kwargs={'label':'%s Dry Deposition (%s/month)' %
                                               ('O$_3$', data.attrs['units']),'fraction':0.046, 'pad':0.04,'aspect':30})
    ax.set_title(title)
    
fig1 = plt.figure(1, figsize=(8,12))
fig1.canvas.set_window_title("global_ozone_drydep")
ax11 = plt.subplot(311, projection=cp.crs.PlateCarree())
ax12 = plt.subplot(312, projection=cp.crs.PlateCarree())
ax13 = plt.subplot(313, projection=cp.crs.PlateCarree())

plot_global_ozone_drydep(ax11, ozone_data_cyclic_da[0], title='Old Scheme')
plot_global_ozone_drydep(ax12, ozone_data_cyclic_da[1], title='Modified EMEP') 
#(100*((ozone_data_cyclic_da[each].sum(dim='time')-ozone_data_cyclic_da[0].sum(dim='time'))/ozone_data_cyclic_da[0].sum(dim='time'))).plot(ax=ax13, levels=np.arange(-100,100.1,10), extend='both', transform=cp.crs.PlateCarree(), cbar_kwargs={'label':'%s Dry Deposition (%s)' % ('O$_3$', '%') ,'fraction':0.046, 'pad':0.04,'aspect':30})

# Subplots on a map
fig2 = plt.figure(2, figsize=(16,10))
fig2.canvas.set_window_title("ozone_drydep_scandinavia")
ax21 = plt.subplot(131, projection=cp.crs.NorthPolarStereo())
ax22 = plt.subplot(132, projection=cp.crs.NorthPolarStereo())
ax23 = plt.subplot(133, projection=cp.crs.NorthPolarStereo())
for ax in fig2.axes:
    # Limit the map to 60 degrees latitude and above.
    ax.set_extent([0, 22, 56, 80], cp.crs.PlateCarree())
    # Adding lines
    draw_parallels(ax, np.arange(60,91,10))
    draw_meridians(ax, np.arange(-180,181,45))    
    ax.coastlines(resolution='50m')

for each in ozone_data_cyclic_da:
   each.sum(dim='time').plot(ax=ax21,
                             transform=cp.crs.PlateCarree(),
                             levels=np.arange(0,20.1), extend='max',
                             cbar_kwargs={'label':'%s Dry Deposition (%s)' %
                                          ('O$_3$', each.attrs['units']),'fraction':0.046, 'pad':0.04,'aspect':30})
    
#(100*((ozone_data_cyclic_da[each].sum(dim='time')-ozone_data_cyclic_da[0].sum(dim='time'))/ozone_data_cyclic_da[0].sum(dim='time'))).plot(ax=ax23, transform=cp.crs.PlateCarree(), levels=np.arange(-100,100.1,10), extend='both', cbar_kwargs={'label':'%s Dry Deposition (%s)' % ('O$_3$', '%'),'fraction':0.046, 'pad':0.04,'aspect':30})

ax21.set_title("Old Scheme")
ax22.set_title("EMAP Scheme")
ax23.set_title("$\Delta O_3^{dd}$")


fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("zonal_ozone_drydep")
ax31 = plt.subplot()

for each in ozone_zonal:
        each.plot()
ax31.set_xlabel("Latitude (deg)")
ax31.set_ylabel("O$_3$ Dry Deposition (Tg Lat$^{-1}$ month$^{-1}$)")
ax31.legend()

#fig4 = plt.figure(4, figsize=(16,9))
#ax41 = plt.subplot()
#(ozone_data[0]*giga/(24*60**2*gridarea*molarweight.sel(name='O3'))/nano).sel(lat=slice(70,71), lon=slice(339,341)).plot()
#ax41.set_ylabel("$\Phi_{O_3}^{DD} (nmol m^{-2} s^{-1})$")
#ax41.set_xlabel("Time (days)")

# Show it
plt.show(block=False)



