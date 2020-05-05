import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

def compute_column_density(press, atm_var, **kargs):
    '''
    Compute the column density.
    Parameters
    ----------
    press : xarray
        The atmospheric pressure field in hPa!
    atm_var : xarray
        The atmospheric tracer field in mol/mol.
    Keyword arguments
    -----------------
    unit : sting
        Define the unit for the output.
        Standard: molecules/cm3
        Option: DU
    accumulate : boolean
        Accumulate all levels.
        Standard: True
    
    '''
    # Keyword arguments
    unit = kargs.pop('unit', 'molecules/cm2')
    b_accu = kargs.pop('accumulate', True)
    # Conversion factors
    Mdry = 0.0289645    # molec. wt. of dry air, kg/mol
    constant = N_A/(g*Mdry) # molecules*s2/(m*kg)
    cm2 = 1e-4
    DU = 1e-2/2.69e16
    # Computation
    # If press.data is not used, labels may interfer in regrouping.
    # Shift the pressure field and multiply it with the unshifted atmospheric variable,
    # but start with second and end one before last.
    p1 = press.data[2:]
    p2 = press.data[:-2]
    if np.all(np.equal(data['O3'].shape, data['lev'].shape)):
        atm_new = 0.5*(p1-p2)*atm_var[1:-1]
    else:
        atm_new = xr.DataArray.copy(atm_var)
        for i in range(len(p1)):
            atm_new[:,i,:,:] = np.fabs(0.5*(p1[i]-p2[i]))*atm_var[:,i+1,:,:]
    # Sum it up == integration
    if b_accu:
        atm_new = atm_new.sum(dim='lev')
    if unit == 'DU':
        atm_new = atm_new*constant*DU
    else:
        atm_new = atm_new*constant*cm2
    # Set the new units
    atm_new.attrs['units'] = unit
    return atm_new

def draw_parallels(ax, parallels, **kargs):
    pc = 66.57 # deg
    b_pc = kargs.pop('polarcircle', False)
    resolution = kargs.pop('resolution', 100)
    lon = np.linspace(-180, 180, resolution)
    for ilat in parallels:
        lat = np.linspace(ilat, ilat, resolution)
        ax.plot(lon, lat, color='lightgrey', linewidth=.8, transform=cp.crs.PlateCarree())
    if b_pc:
        ax.plot(lon, np.linspace(pc,pc,resolution), color='lightgrey', linewidth=.8, ls='--', transform=cp.crs.PlateCarree())
        ax.plot(lon, np.linspace(-pc,-pc,resolution), color='lightgrey', linewidth=.8, ls='--', transform=cp.crs.PlateCarree())
def draw_meridians(ax, meridians, **kargs):
    polarcap = kargs.pop('polarcap', 80)
    resolution = kargs.pop('resolution', 100)
    lat = np.linspace(-polarcap, polarcap, resolution)
    for ilon in meridians:
        lon = np.linspace(ilon,ilon, resolution)
        ax.plot(lon, lat, color='lightgrey', linewidth=.8, transform=cp.crs.PlateCarree())
    ax.plot(np.linspace(-180, 180, resolution),
            np.linspace(polarcap, polarcap, resolution),
            color='lightgrey', linewidth=.8, transform=cp.crs.PlateCarree())
    ax.plot(np.linspace(-180, 180, resolution),
            np.linspace(-polarcap, -polarcap, resolution),
            color='lightgrey', linewidth=.8, transform=cp.crs.PlateCarree())
    # Adding text
    for ilon in np.unique(np.fabs(meridians)):
        ax.text(-ilon, 50, '%d$^\circ$W' % (ilon), horizontalalignment='center', transform=cp.crs.Geodetic())

    
    
# Data source
nc_src = os.environ['DATA']+'/nird_data/models/results/OsloCTM3/drydepdevel/version2/C3RUN_default/monthly_means/avgsav_20050101_20050201.nc'

try:
    data
except NameError:
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src)):
        print(file)
        data = xr.open_dataset(file)
        data_list.append(data)
# Concatinating the listq
data = xr.concat(data_list, dim='time')
molweights = get_molarweight(data)
ozone = get_vmr(data,tracer='O3', unit='ppm')
ozone_zonalmean = ozone.mean(dim='lon')
ozone_zonalmean.attrs['units'] = ozone.attrs['units']

ozone_column = compute_column_density(data['lev'], get_vmr(data,tracer='O3'), unit='DU')
# Unaccumulated ozone field in DU.
ozone_du = compute_column_density(data['lev'], get_vmr(data,tracer='O3'), unit='DU', accumulate=False)

#Plotting
# Clean-up
plt.close('all')
levels = np.arange(0,10.1,0.25)
fig1 = plt.figure(1, figsize=(16,9))
ozone_zonalmean.isel(time=0).plot.contourf(cmap=plt.cm.afmhot_r, levels=levels)
ax11 = plt.gca()
set_pressure_axis(ax11, limits=(1000,0.1), ticks=(1000,500,200,100,50,20,10,5,2,1,.5,.2,.1))
ax11.set_xlabel("Lat (deg)")

# Subplots on a map
levelsDU = np.arange(200,501,10)
fig2 = plt.figure(2, figsize=(12,12))
ax21 = plt.subplot(211, projection=cp.crs.PlateCarree())
ax22 = plt.subplot(212, projection=cp.crs.NorthPolarStereo())
#fig2, (ax21,ax22) = plt.subplots(nrows=2, subplot_kw={'projection': cp.crs.PlateCarree()}) #
#fig2.set_size_inches(16,9)

ozone_column_cyclic, cyclic_lon = ccrs_util.add_cyclic_point(ozone_column.data, ozone_column['lon'])
ozone_column_cyclic_da = xr.DataArray(ozone_column_cyclic, coords=[ozone_column['time'], ozone_column['lat'], cyclic_lon], dims=['time','lat', 'lon'])
ozone_column_cyclic_da.attrs['units'] = ozone_column.attrs['units']

# Set the map and axis attributes
for ax in fig2.axes:
    ax.set_global()   # Expands the map to fit the tick labels!
    ax.coastlines()
    #ax.gridlines()
ax21.set_xticks(np.arange(-180, 181, 45), crs=cp.crs.PlateCarree())
ax21.set_yticks(np.arange(-80, 81, 20), crs=cp.crs.PlateCarree())
lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()
ax21.xaxis.set_major_formatter(lon_formatter)
ax21.yaxis.set_major_formatter(lat_formatter)

# Limit the map to 60 degrees latitude and above.
ax22.set_extent([-180, 180, 90, 60], cp.crs.PlateCarree())

# Adding lines
draw_parallels(ax22, np.arange(60,91,10))
draw_meridians(ax22, np.arange(-180,181,45))    
    
ozone_column_cyclic_da.isel(time=0).plot.contourf(ax=ax21, levels=levelsDU, transform=cp.crs.PlateCarree(), cbar_kwargs={'label':'O$_3$ (DU)','fraction':0.046, 'pad':0.04,'aspect':30})
ozone_column_cyclic_da.isel(time=0).plot.contourf(ax=ax22, levels=levelsDU, transform=cp.crs.PlateCarree(), cbar_kwargs={'label':'O$_3$ (DU)','fraction':0.046, 'pad':0.04,'aspect':30})

# Show it
plt.show(block=False)



