import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
#import cartopy.geodesic as cgeo
from matplotlib.transforms import offset_copy
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *


# Data source
nc_src = os.environ['DATA']+'/CTM3_input_data/EMIS/ECLIPSE_V5a/ship_CLE_emis_*.nc'
nc_src_historic = os.environ['DATA']+'/CTM3_input_data/EMIS/CEDS0517/CEDS_CICERO/'
nc_src_ssp585 = os.environ['DATA']+'/astra_data/input_data/ctm_input/EMIS/CEDS_SCENARIOS/SSP585/'
var = ('SO2_em_anthro',
       'NOx_em_anthro', 
       'CO_em_anthro', 
       #'CO2_em_anthro', 
       'NH3_em_anthro', 
       'VOC02_ethane_em_speciated_VOC_anthro', 
       'VOC03_propane_em_speciated_VOC_anthro', 
       'VOC04_butanes_em_speciated_VOC_anthro', 
       'VOC05_pentanes_em_speciated_VOC_anthro', 
       'VOC06_hexanes_pl_em_speciated_VOC_anthro', 
       'VOC07_ethene_em_speciated_VOC_anthro', 
       'VOC08_propene_em_speciated_VOC_anthro', 
       'VOC13_benzene_em_speciated_VOC_anthro', 
       'VOC14_toluene_em_speciated_VOC_anthro', 
       'VOC15_xylene_em_speciated_VOC_anthro', 
       'VOC16_trimethylb_em_speciated_VOC_anthro', 
       'VOC17_other_arom_em_speciated_VOC_anthro', 
       'VOC21_methanal_em_speciated_VOC_anthro', 
       'VOC22_other_alka_em_speciated_VOC_anthro', 
       'VOC23_ketones_em_speciated_VOC_anthro', 
       'BC_em_FOSSIL_FUEL_anthro', 
       'BC_em_SOLID_BIOFUEL_anthro', 
       'OC_em_FOSSIL_FUEL_anthro', 
       'OC_em_SOLID_BIOFUEL_anthro')

try:
    data
except NameError:
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src)):
        print(file)
        tmp_data = xr.open_dataset(file)
        data_list.append(tmp_data)
    # Concatinating the list
    data = xr.concat(data_list, dim='time')
    # Reding CEDs datasets
    data_scandic_slice_ceds = []
    data_finnmark_slice_ceds = []
    for ivar in var[:4]:
        tmp_list = []
        # Reading the historical data
        for file in sorted(glob.glob(nc_src_historic+ivar.replace('_','-')+'*201412*.nc')):
            print(os.path.basename(file))
            tmp_data = xr.open_dataset(file)
            # Defining new time coordinates
            tmp_data['time'].reset_coords(drop=True)
            tmp_data.coords['time'] = [dt.datetime(iyear, imonth, iday) for iyear, imonth, iday in zip(tmp_data.time.dt.year, tmp_data.time.dt.month, tmp_data.time.dt.day)]            
            # Get gridarea (there is supposed to be a file on the CEDS servers)
            try:
                gridarea
            except NameError:
                gridarea = tmp_data["Gridcell area"]
            tmp_list.append(tmp_data[ivar+'_SHI'])
        # Reading the scenario data
        for file in sorted(glob.glob(nc_src_ssp585+ivar.replace('_','-')+'*12.nc')):
            print(os.path.basename(file))
            tmp_data = xr.open_dataset(file)
            # Defining new time coordinates
            tmp_data['time'].reset_coords(drop=True)
            tmp_data.coords['time'] = [dt.datetime(iyear, imonth, iday) for iyear, imonth, iday in zip(tmp_data.time.dt.year, tmp_data.time.dt.month, tmp_data.time.dt.day)] 
            tmp_list.append(tmp_data[ivar+'_SHI'].drop('sector'))
        # Concat the species and ad it to list
        tmp_data = xr.concat(tmp_list, dim='time')
        data_scandic_slice_ceds.append((tmp_data*gridarea).sel(lat=slice(0,40),lon=slice(56,81),drop=True))
        data_finnmark_slice_ceds.append((tmp_data*gridarea).sel(lat=slice(16,33),lon=slice(66,71),drop=True))
        #data_list.append(tmp_data)

    # Select Scandinavia
    data_scandic_slice = data.sel(lat=slice(0,40),lon=slice(56,81),drop=True)
    data_finnmark_slice = data.sel(lat=slice(16,33),lon=slice(66,71),drop=True)
    x_scandic, y_scandic = [0,40,40,0,0], [56,56,81,81,56]
    x, y = [16,33,33,16,16], [66,66,71,71,66]

# Get seconds in month for all years
try:
    siy
except NameError:
    siy = []
    for iyear in np.unique(data_scandic_slice_ceds[0].time.dt.year):
        siy.append([seconds_in_month(imonth, iyear) for imonth in np.arange(1,13)]) 

# Clean-up
plt.close('all')

#Plotting

fig1 = plt.figure(1, figsize=(7,14))
fig1.canvas.set_window_title("map_scandic_domain")
stamen_terrain = cimgt.Stamen('terrain-background')
#ax11 = fig1.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax11 = fig1.add_subplot(1,1,1, projection=stamen_terrain.crs) #ccrs.PlateCarree()
#ax41.set_extent([18,31,66,70], crs=ccrs.PlateCarree())
ax11.set_extent([0,40,56,81], crs=ccrs.PlateCarree())
# Add the Stamen data at zoom level 8.
ax11.add_image(stamen_terrain, 4)
#ax11.stock_img()
ax11.coastlines(resolution='10m')

# Mark the inner domain
ax11.fill(x, y, color='coral', alpha=0.4, transform=ccrs.Geodetic())

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("accumulated_shipping_emissions_ECLIPSE")
#ax21 = plt.subplot()

for i, key in enumerate(data.data_vars):
    ax = plt.subplot(3,3,i+1)
    (100*(data_scandic_slice[key].sum(dim=('lat','lon'))-data_scandic_slice[key].sum(dim=('lat','lon')).sel(time=2005))/data_scandic_slice[key].sum(dim=('lat','lon')).sel(time=2005)).plot(ax=ax, label="Fennoscandia")
    (100*(data_finnmark_slice[key].sum(dim=('lat','lon'))-data_finnmark_slice[key].sum(dim=('lat','lon')).sel(time=2005))/data_finnmark_slice[key].sum(dim=('lat','lon')).sel(time=2005)).plot(ax=ax, label="Cap of the North")
    ax.set_xlabel("")
    ax.set_ylabel("")
    #ax.set_ylabel("%s (%s)" % (key, data_scandic_slice[key].units))
    ax.set_title("%s" % (key), y=0.85, x=0.15)

ax28 = fig2.axes[-2]
ax24 = fig2.axes[3]
ax28.legend(bbox_to_anchor=(0.35, -0.39), loc=8, borderaxespad=0.,ncol=4)
ax28.set_xlabel("Time (years)")
ax24.set_ylabel("Emissions relative to 2005 (%)") #(%s)" % (data_scandic_slice[key].units))

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("accumulated_shipping_emissions_CEDS")
for i in np.arange(len(data_scandic_slice_ceds)):
    ax = plt.subplot(2,2,i+1)
    #((data_scandic_slice_ceds[i].sum(dim=('lat','lon'))*1e-3)*np.array(siy).ravel()).plot(ax=ax)
    #((data_finnmark_slice_ceds[i].sum(dim=('lat','lon'))*1e-3)*np.array(siy).ravel()).plot(ax=ax)
    ref_scandic = data_scandic_slice_ceds[i].sum(dim=('lat','lon')).groupby("time.year").mean().sel(year=2005)
    ref_finnmark = data_finnmark_slice_ceds[i].sum(dim=('lat','lon')).groupby("time.year").mean().sel(year=2005)
    (100*(data_scandic_slice_ceds[i].sum(dim=('lat','lon'))-ref_scandic)/ref_scandic).plot(ax=ax, marker='x', ls='None', label="Fennoscandia")
    (100*(data_finnmark_slice_ceds[i].sum(dim=('lat','lon'))-ref_finnmark)/ref_finnmark).plot(ax=ax, marker='+', ls='None', label="Cap of the North")
    
    ax.set_title(var[i].replace('_em_anthro',''), y=0.85, x=0.1)
    ax.set_xlabel('')
    

ax33 = fig3.axes[-2]
ax33.set_xlabel("Time (years)", x=1.15)
ax33.set_ylabel("Emissions relative to 2005 (%)", y=1.15)
fig3.axes[0].legend(loc="upper right", borderaxespad=0.,ncol=4)


# Show it
plt.show(block=False)
