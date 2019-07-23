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
from mytools.line_styles import linestyles
from collections import OrderedDict

# Clean up
plt.close('all')

# Read the data points from Hardacre et al. 2015
execfile("hardacre_data.py")

b_norm = True

# Data source
data_dir = os.environ['DATA']+'/astra_data/ctm_results/' 

pft_data_dir = os.environ['DATA']+'/astra_data/processed_data/pft_landusedyn.3x3.nc'

scav_dir = "scavenging_monthly/hardacre_grid/"
mm_dir = "monthly_means/hardacre_grid/"

experiment = (#'C3RUN_default/',
              'C3RUN_mOSaic/',
              'C3RUN_mOSaic_offLight/',
              'C3RUN_mOSaic_offPhen/',
              'C3RUN_mOSaic_SWVL1/',
              'C3RUN_mOSaic_ice/',
              'C3RUN_mOSaic_desert/',
              'C3RUN_mOSaic_emis2014/',
              'C3RUN_mOSaic_hough/'
)


labels = (#'Wesely_type',
          'mOSaic',
          'mOSaic_offLight','mOSaic_offPhen','mOSaic_SWVL1',
          'mOSaic_ice',
          'mOSaic_desert',
          'mOSaic_emis2014',
          'mOSaic_hough')
colors = np.concatenate((#('black',),
                         ('blue',),
                         np.repeat('cornflowerblue',3),
                         np.repeat('purple',2),
                         ('grey',),
                         ('orange',)) )

# Read the data
try:
    data
except NameError:
    vdd_raw_data = []
    for iexp in experiment:
        subdir = data_dir+iexp+scav_dir+'mm_vo3*.nc'
        print("Reading from path %s" % (os.path.abspath(subdir)))
        data_list = []
        # Open dataset
        for file in sorted(glob.glob(subdir)):
            print("Reading %s" % (os.path.basename(file)))
            data = xr.open_dataset(file)
            # Define new time coordinates and drop the old one
            data['time'].reset_coords(drop=True)
            year = int(file[-12:-8])
            month = int(file[-8:-6])
            data.coords['time'] = ([dt.datetime(year, month, 1),])
            if (data.lat.ndim > 1):
                print("Changing coordinates...")
                data.coords['x'] = data.lat.lon[0].data
                data.coords['y'] = data.lat[:,0].data
                data = data.drop(('lon','lat'))
                data = data.rename({'x':'lon','y':'lat'})
            data_list.append(data['VO3'])
        # Concatenating the list
        vdd_raw_data.append(xr.concat(data_list, dim='time'))

    # Read plant function types
    pft_data = read_data(pft_data_dir)
    # Get PFT selection for use as masks
try:
    pft_sel
except NameError:
    pft_sel = []
    for each in pft_data['PFT_PCT'].sel(time=2005.00):
        print("Changing coordinates...")
        each.coords['x'] = each.lat.lon[0].data
        each.coords['y'] = each.lat[:,0].data
        each = each.drop(('lon','lat'))
        each = each.rename({'x':'lon','y':'lat'})
        pft_sel.append(each)

#"0: Barren land";
#"1: Needleaf evergreen temperate tree";
#"2: Needleaf evergreen boreal tree";
#"3: Needleaf deciduous boreal tree";
#"4: Broadleaf evergreen tropical tree";
#"5: Broadleaf evergreen temperate tree";
#"6: Broadleaf deciduous tropical tree";
#"7: Broadleaf deciduous temperate tree";
#"8: Broadleaf deciduous boreal tree";
#"9: Broadleaf evergreen temperate shrub";
#"10: Broadleaf deciduous temperate shrub";
#"11: Broadleaf deciduous boreal shrub";
#"12: Arctic C3 grass (cold)";
#"13: C3 grass (cool)";
#"14: C4 grass (warm)";
#"15: Crop1";
#"16: Crop2";

# mOSaic categories (reduced from Simpson et al., 2012)
    #  1. Needleleaftree (temperated/boreal)
    #  2. Deciduoustree (temperated/boral)
    #  3. Needleleaftree (mediterranean)
    #  4. Broadleaftree (mediterranean)
    #  5. Crops <a,b,c>
    #  6. Moorland (savanna++)
    #  7. Grassland
    #  8. Scrubs (med.)
    #  9. Wetlands
    #  10. Tundra
    #  11. Desert
    #  12. Water
    #  13. Urban
    #  14. Ice/Snow - is treated seperately

pft_de = pft_sel[0].where((pft_sel[0].lat >=-60) & (pft_sel[0].lat <= 70))
pft_is = pft_sel[0].where(~((pft_sel[0].lat >=-60) & (pft_sel[0].lat <= 70)))
pft_oc = (xr.concat(pft_sel,dim='pft').sum(dim='pft')-100)*-1
pft_cf = xr.concat(pft_sel[1:4], dim='pft').sum(dim='pft') #xr.concat(pft_sel[10], dim='pft')
pft_df = xr.concat(pft_sel[7:11], dim='pft').sum(dim='pft') #xr.concat(pft_sel[11:13], dim='pft')
pft_tf = xr.concat(pft_sel[4:7], dim='pft').sum(dim='pft')
pft_ac = pft_sel[15]
pft_gr = xr.concat(pft_sel[13:15], dim='pft').sum(dim='pft')
pft_tu = xr.concat(pft_sel[11:13], dim='pft').sum(dim='pft')

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("final-total_ozone_drydepvelo_monthly_pft_v2")

ax11 = plt.subplot(331)
ax12 = plt.subplot(332)
ax13 = plt.subplot(333)
ax14 = plt.subplot(334)
ax15 = plt.subplot(335)
ax16 = plt.subplot(336)
ax17 = plt.subplot(337)
ax18 = plt.subplot(338)
ax19 = plt.subplot(339)

# Coniferous forest
for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if j < len(vdd_raw_data):
        ((vdd_raw_data[j].sel(NLCAT=1)*1e2)).where((pft_cf>=96),drop=True).mean(dim=('lat','lon')).plot(ax=ax11, label=labels[j], color=colors[j], ls=linestyle)

# Decidous forest
for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if j < len(vdd_raw_data):
        ((vdd_raw_data[j].sel(NLCAT=2)*1e2)).where((pft_df>=72),drop=True).mean(dim=('lat','lon')).plot(ax=ax12, label=labels[j], color=colors[j], ls=linestyle)

# Tropical forest
for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if j < len(vdd_raw_data):
        ((vdd_raw_data[j].sel(NLCAT=4)*1e2)).where((pft_tf>=100),drop=True).mean(dim=('lat','lon')).plot(ax=ax13, label=labels[j], color=colors[j], ls=linestyle)
        print(((vdd_raw_data[j].sel(NLCAT=4)*1e2)).where((pft_tf>=100),drop=True).mean(dim=('lat','lon')).mean(dim='time'))
        
# Cropland
for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if j < len(vdd_raw_data):
        ((vdd_raw_data[j].sel(NLCAT=5)*1e2)).where((pft_ac>=70)&(vdd_raw_data[j].lat>=40),drop=True).mean(dim=('lat','lon')).plot(ax=ax14, label=labels[j], color=colors[j], ls=linestyle)
        
# Grassland
for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if j < len(vdd_raw_data):
        #((vdd_raw_data[j].sel(NLCAT=7)*1e2)).where((pft_gr>=77)&(vdd_raw_data[j].lat>=0),drop=True).mean(dim=('lat','lon')).plot(ax=ax15, label=labels[j], color=colors[j], ls=linestyle)
        ((vdd_raw_data[j].sel(NLCAT=7)*1e2)).where((vdd_raw_data[j].lat>=0),drop=True).max(dim=('lat','lon')).plot(ax=ax15, label=labels[j], color=colors[j], ls=linestyle)
        print(((vdd_raw_data[j].sel(NLCAT=7)*1e2)).where((vdd_raw_data[j].lat>=0),drop=True).max(dim=('lat','lon')).mean(dim='time'))
                
# Tundra
for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if j < len(vdd_raw_data):
        ((vdd_raw_data[j].sel(NLCAT=10)*1e2)).where((pft_tu>=80)&(vdd_raw_data[j].lat>=0),drop=True).mean(dim=('lat','lon')).plot(ax=ax16, label=labels[j], color=colors[j], ls=linestyle)
        
# Snow and ice (some other barren ground?)
for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if j < len(vdd_raw_data):
        ((vdd_raw_data[j].sel(NLCAT=14)*1e2)).where((pft_is>=100)&(vdd_raw_data[j].lat>=0),drop=True).mean(dim=('lat','lon')).plot(ax=ax17, label=labels[j], color=colors[j], ls=linestyle)
            
# Ocean
for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if j < len(vdd_raw_data):
        ((vdd_raw_data[j].sel(NLCAT=12)*1e2)).where((pft_oc>=100)&(vdd_raw_data[j].lat>=0),drop=True).mean(dim=('lat','lon')).plot(ax=ax18, label=labels[j], color=colors[j], ls=linestyle)
        
# Desert
for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if j < len(vdd_raw_data):
        ((vdd_raw_data[j].sel(NLCAT=11)*1e2)).where((pft_de>=100),drop=True).mean(dim=('lat','lon')).plot(ax=ax19, label=labels[j], color=colors[j], ls=linestyle)
        print(((vdd_raw_data[j].sel(NLCAT=11)*1e2)).where((pft_de>=100),drop=True).mean(dim=('lat','lon')).mean(dim='time'))
        


ax11.scatter(vdd_raw_data[0].time.data, hardacre_cf, label='Hardacre (2015)', color='red')
ax12.scatter(vdd_raw_data[0].time.data, hardacre_df, label='Hardacre (2015)', color='red')
ax13.scatter(vdd_raw_data[0].time.data, hardacre_tf, label='Hardacre (2015)', color='red')
ax14.scatter(vdd_raw_data[0].time.data, hardacre_ac, label='Hardacre (2015)', color='red')
ax15.scatter(vdd_raw_data[0].time.data, hardacre_gr, label='Hardacre (2015)', color='red')
ax16.scatter(vdd_raw_data[0].time.data, hardacre_tu, label='Hardacre (2015)', color='red')
ax17.scatter(vdd_raw_data[0].time.data, hardacre_is, label='Hardacre (2015)', color='red')
ax18.scatter(vdd_raw_data[0].time.data, hardacre_oc, label='Hardacre (2015)', color='red')
ax19.scatter(vdd_raw_data[0].time.data, hardacre_de, label='Hardacre (2015)', color='red')

ax11.set_title("Coniferous forest", loc='right', x=0.99, y=0.85)
ax12.set_title("Deciduous forest", loc='right', x=0.99, y=0.85)
ax13.set_title("Tropical forest", loc='right', x=0.99, y=0.85)
ax14.set_title("Cropland", loc='right', x=0.99, y=0.85)
ax15.set_title("Grassland", loc='right', x=0.99, y=0.85)
ax16.set_title("Tundra", loc='right', x=0.99, y=0.85)
ax17.set_title("Ice/Snow", loc='right', x=0.99, y=0.85)
ax18.set_title("Ocean", loc='right', x=0.99, y=0.85)
ax19.set_title("Desert", loc='right', x=0.99, y=0.85)

for ax in fig1.axes:
    ticks = ax.get_xticks()
    start = ticks[0]-(ticks[1]-ticks[0])/2.
    end = ticks[-1]
    ax.set_title("")
    ax.set_xticks(np.linspace(start, end, 12))
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_ylim(0,1.5)
    ax.set_xticklabels("")
for ax in fig1.axes[-3:]:
    ax.set_ylim(0,0.3)
    
ax14.set_ylabel("$v^{O_3}_{DD}$ (%s)" % ("$cm\,s^{-1}$"))
#ax18.set_xlabel("Time (months)")
ax17.set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
ax18.set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
ax19.set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
ax18.legend(bbox_to_anchor=(0.45, -0.4), loc=8, borderaxespad=0.,ncol=5)

# Zonal bands
fig2 = plt.figure(2, figsize=(16,9))

ax2 = plt.subplot()

dd_velo_zonal = [each.sum(dim='NLCAT').mean(dim='lon')*1e2 for each in vdd_raw_data]

if b_norm:
    fig2.canvas.set_window_title("final-norm_total_ozone_drydepvelo_v2")
else:
    fig2.canvas.set_window_title("final-total_ozone_drydepvelo_v2")

for i, (name,linestyle) in enumerate(linestyles.items()[1:]):
    if i < len(dd_velo_zonal):
        if b_norm:
            (dd_velo_zonal[i].mean(dim='time')/dd_velo_zonal[i].mean(dim='time').sum()).plot(label=labels[i], color=colors[i], ls=linestyle)
        else:
            (dd_velo_zonal[i].mean(dim='time')).plot(label=labels[i], color=colors[i], ls=linestyle)
            
ax2.set_title('')
ax2.set_xlabel("Latitude (deg)")

if b_norm:
    ax2.set_ylabel("$v^{O_3}_{DD} / \Sigma v^{O_3}_{DD}$")
    ax2.scatter(hardacre_lat, hardacre_ddvel/hardacre_ddvel.sum(), label='mmm Hardacre (2015)', color='red')
    ax2.set_ylim(0,0.04)
    ax2.legend(ncol=1)
else:
    ax2.set_ylabel("$v^{O_3}_{DD}$ (%s)" % ("$cm\,s^{-1}$"))
    ax2.scatter(hardacre_lat, hardacre_ddvel, label='mmm Hardacre (2015)', color='red')
    ax2.set_ylim(0,0.4)
    ax2.legend(ncol=2)

# Total ozone dry deposition veloscity over latitudes split by month
fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("final-total_ozone_drydepvelo_monthly_v2")

rows, col = 4, 3
for i in np.arange(1,rows*col+1):
    plt.subplot(rows, col, i)
    ax =  plt.gca()
    for j, (name,linestyle) in enumerate(linestyles.items()[1:]):
        if j < len(dd_velo_zonal):
            ((dd_velo_zonal[j]).isel(time=i-1)).plot(label=labels[j], color=colors[j], ls=linestyle)
                    
    if i==2:
        ax.scatter(hardacre_lat, hardacre_ddvel_feb, label='Hardacre (2015)', color='red')
    elif i==8:
        ax.scatter(hardacre_lat, hardacre_ddvel_aug, label='Hardacre (2015)', color='red')
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_ylim(0,0.6)
    ax.set_title("")
    ax.set_title(get_month_name(i), x=0.05, y=0.8, loc='left')
fig3.axes[-2].set_xlabel("Latitude (deg)")
fig3.axes[6].set_ylabel("$v^{O_3}_{DD}$ (%s)" % ("$cm\,s^{-1}$"), y=1)

# Show it
plt.show(block=False)
