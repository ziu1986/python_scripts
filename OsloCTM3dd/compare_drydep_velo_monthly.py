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
from collections import OrderedDict

# Read the data points from Hardacre et al. 2015
execfile("hardacre_data.py")

b_norm = False

# Data source
scav_dir = "scavenging_monthly/regrid_hardacre/"
mm_dir = "monthly_means/regrid_hardacre/"

experiment = ('C3RUN_oDD/',
              'C3RUN_emep_full/',
              'C3RUN_emep_offLight/',
              'C3RUN_emep_offPhen/',
              'C3RUN_emep_SWVL4/',
              'C3RUN_emep_ppgs/',
              'C3RUN_emep_ppgssh/',
              'C3RUN_emep_ppgssh_ice/',
              'C3RUN_emep_ppgs_2005/')
data_dir = os.environ['DATA']+'/astra_data/ctm_results/' 

pft_data_dir = "./pft_landusedyn.3x3.nc"

labels = ('OsloCTM3: Wesely type',
          #'OsloCTM3: EMEP/MEGAN_corr','OsloCTM3: EMEP','OsloCTM3: EMEP_swgd',
          'OsloCTM3: EMEP_full',
          'OsloCTM3: EMEP_offLight','OsloCTM3: EMEP_offPhen','OsloCTM3: EMEP_SWVL4',
          'OsloCTM3: EMEP_ppgs','OsloCTM3: EMEP_ppgssh',
          'OsloCTM3: EMEP_ppgssh_ice',
          'OsloCTM3: EMEP_ppgs_2005')
colors = np.concatenate((('black',),
                         ('blue',),
                         np.repeat('cornflowerblue',3),
                         np.repeat('purple',2),
                         ('grey',),
                         ('orange',)) )

# Read the data
try:
    data
except NameError:
    dry_dep_raw_data = []
    for iexp in experiment:
        subdir = data_dir+iexp+scav_dir+'*.nc'
        print("Reading from path %s" % (os.path.abspath(subdir)))
        data_list = []
        # Open dataset
        for file in sorted(glob.glob(subdir)):
            print("Reading %s" % (os.path.basename(file)))
            data = xr.open_dataset(file)
            # Define new time coordinates and drop the old one
            data['time'].reset_coords(drop=True)
            data.coords['time'] = ([dt.datetime(data['YEAR'], data['MONTH'], 1),])
            if (data.lat.ndim > 1):
                print("Changing coordinates...")
                data.coords['x'] = data.lat.lon[0].data
                data.coords['y'] = data.lat[:,0].data
                data = data.drop(('lon','lat'))
                data = data.rename({'x':'lon','y':'lat'})
            data_list.append(data['dry_O3'])
        # Concatenating the list
        dry_dep_raw_data.append(xr.concat(data_list, dim='time'))
    # Extracting some general information
    # WARNING: cdo has summed these also up -> devide by number of days
    gridarea = data['gridarea'].isel(time=0)/31. 
    #molarweight = get_molarweight(data.isel(time=0))/31.
    
    ozone_raw_data = []
    height_raw_data = []
    for iexp in experiment:
        subdir = data_dir+iexp+mm_dir+'*.nc'
        print("Reading from path %s" % (os.path.abspath(subdir)))
        data_list = []
        data_list_height = []
        # Open dataset
        for file in sorted(glob.glob(subdir)):
            print("Reading %s" % (os.path.basename(file)))
            data = xr.open_dataset(file)
            # Define new time coordinates and drop the old one
            data.coords['time'] = ([dt.datetime.strptime(os.path.basename(file)[23:31], '%Y%M%d'),])
            if (data.lat.ndim > 1):
                print("Changing coordinates...")
                data.coords['x'] = data.lat.lon[0].data
                data.coords['y'] = data.lat[:,0].data
                data = data.drop(('lon','lat'))
                data = data.rename({'x':'lon','y':'lat'})
            data_list.append(data['O3'].isel(lev=0))
            data_list_height.append(data['height'].isel(ilev=[0,1]))
        # Concatenating the list
        ozone_raw_data.append(xr.concat(data_list, dim='time'))
        height_raw_data.append(xr.concat(data_list_height, dim='time'))
        
    # Read plant function types
    pft_data = read_data(pft_data_dir)

# Compute dry deposition velocity
sim = np.array([seconds_in_month(i, np.unique(dry_dep_raw_data[0].time.dt.year)) for i in range(1,13)])
dim = sim/(24*60**2)
layer_thickness = (height_raw_data[0].isel(ilev=1)-height_raw_data[0].isel(ilev=0))*100 #[1m -> 100cm] Lowermost level center 8 m!
dd_velo = []

#for i in np.arange(len(dry_dep_raw_data)):
    #dd_velo.append((dry_dep_raw_data[i]/ozone_raw_data[i]))
for i in np.arange(len(dry_dep_raw_data)):
    tmp_velo = dry_dep_raw_data[i]/ozone_raw_data[i]
    tmp_velo_list = []
    for imonth in range(0,12):
        tmp_velo_list.append(tmp_velo.where(tmp_velo.time.dt.month==imonth+1, drop=True)/sim[imonth]*layer_thickness.isel(time=imonth))
    dd_velo.append(xr.concat(tmp_velo_list, dim='time'))

    
dd_velo_zonal = [each.mean(dim='lon') for each in dd_velo]
for each in dd_velo_zonal:
    each.attrs['unit'] = 'cms-1'
    
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

pft_de = pft_sel[0].where((pft_sel[0].lat >=-60) & (pft_sel[0].lat <= 70))
pft_is = pft_sel[0].where(~((pft_sel[0].lat >=-60) & (pft_sel[0].lat <= 70)))
pft_oc = (xr.concat(pft_sel,dim='pft').sum(dim='pft')-100)*-1
pft_cf = xr.concat(pft_sel[1:4], dim='pft').sum(dim='pft') #xr.concat(pft_sel[10], dim='pft')
pft_df = xr.concat(pft_sel[7:11], dim='pft').sum(dim='pft') #xr.concat(pft_sel[11:13], dim='pft')
pft_tf = xr.concat(pft_sel[4:7], dim='pft').sum(dim='pft')
pft_ac = pft_sel[15]
pft_gr = xr.concat(pft_sel[13:15], dim='pft').sum(dim='pft')
pft_tu = xr.concat(pft_sel[11:13], dim='pft').sum(dim='pft')
       
# Clean up
plt.close('all')

# Line styles
linestyles = OrderedDict(
    [('solid',               (0, ())),
     #('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 3))),
     ('densely dotted',      (0, (1, 1))),

     #('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     #('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     #('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
if b_norm:
    fig1.canvas.set_window_title("final-norm_total_ozone_drydepvelo")
else:
    fig1.canvas.set_window_title("final-total_ozone_drydepvelo")
ax1 = plt.subplot()

for i, (name,linestyle) in enumerate(linestyles.items()):
    if i < len(dd_velo_zonal):
        if b_norm:
            (dd_velo_zonal[i].mean(dim='time')/dd_velo_zonal[i].mean(dim='time').sum()).plot(label=labels[i], color=colors[i], ls=linestyle)
        else:
            (dd_velo_zonal[i].mean(dim='time')).plot(label=labels[i], color=colors[i], ls=linestyle)
            
ax1.set_title('')
ax1.set_xlabel("Latitude (deg)")

if b_norm:
    ax1.set_ylabel("$v^{O_3}_{DD} / \Sigma v^{O_3}_{DD}$")
    ax1.scatter(hardacre_lat, hardacre_ddvel/hardacre_ddvel.sum(), label='mmm Hardacre (2015)', color='red')
    ax1.set_ylim(0,0.04)
    ax1.legend(ncol=1)
else:
    ax1.set_ylabel("$v^{O_3}_{DD}$ (%s)" % ("$cm\,s^{-1}$"))
    ax1.scatter(hardacre_lat, hardacre_ddvel, label='mmm Hardacre (2015)', color='red')
    ax1.set_ylim(0,0.4)
    ax1.legend(ncol=3)

# Total ozone dry deposition veloscity over latitudes split by month
fig2 = plt.figure(2, figsize=(16,9))
if b_norm:
    fig2.canvas.set_window_title("final-norm_total_ozone_drydepvelo_monthly")
else:
    fig2.canvas.set_window_title("final-total_ozone_drydepvelo_monthly")

rows, col = 4, 3
for i in np.arange(1,rows*col+1):
    plt.subplot(rows, col, i)
    ax =  plt.gca()
    for j, (name,linestyle) in enumerate(linestyles.items()):
        if j < len(dd_velo_zonal):
            if b_norm:
                ((dd_velo_zonal[j]).isel(time=i-1)/(dd_velo_zonal[j].isel(time=i-1).sum())).plot(label=labels[j], color=colors[j], ls=linestyle)
            else:
                ((dd_velo_zonal[j]).isel(time=i-1)).plot(label=labels[j], color=colors[j], ls=linestyle)
                    
    if i==2:
        if b_norm:
            ax.scatter(hardacre_lat, hardacre_ddvel_feb/hardacre_ddvel_feb.sum(), label='Hardacre (2015)', color='red')
        else:
            ax.scatter(hardacre_lat, hardacre_ddvel_feb, label='Hardacre (2015)', color='red')
    elif i==8:
        if b_norm:
            ax.scatter(hardacre_lat, hardacre_ddvel_aug/hardacre_ddvel_aug.sum(), label='Hardacre (2015)', color='red')
        else:
            ax.scatter(hardacre_lat, hardacre_ddvel_aug, label='Hardacre (2015)', color='red')
    ax.set_xlabel("")
    ax.set_ylabel("")
    if b_norm:
        ax.set_ylim(0,0.08)
    else:
        ax.set_ylim(0,0.6)
    ax.set_title("")
    ax.set_title(get_month_name(i), x=0.05, y=0.8, loc='left')
fig2.axes[-2].set_xlabel("Latitude (deg)")
if b_norm:
    fig2.axes[6].set_ylabel("$v^{O_3}_{DD}/\Sigma v^{O_3}_{DD}$", y=1)
else:
    fig2.axes[6].set_ylabel("$v^{O_3}_{DD}$ (%s)" % ("$cm\,s^{-1}$"), y=1)

# Seperation of different PFT types
fig4 = plt.figure(4, figsize=(16,9))
fig4.canvas.set_window_title("final-total_ozone_drydepvelo_monthly_pft")

ax41 = plt.subplot(331)
ax42 = plt.subplot(332)
ax43 = plt.subplot(333)
ax44 = plt.subplot(334)
ax45 = plt.subplot(335)
ax46 = plt.subplot(336)
ax47 = plt.subplot(337)
ax48 = plt.subplot(338)
ax49 = plt.subplot(339)


# Coniferous forest
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        (dd_velo[j].where((pft_cf>=96),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax41, label=labels[j], color=colors[j], ls=linestyle)
        
# Decidous forest
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        (dd_velo[j].where((pft_df>=72),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax42, label=labels[j], color=colors[j], ls=linestyle)
        
# Tropical forest
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        (dd_velo[j].where((pft_tf>=100),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax43, label=labels[j], color=colors[j], ls=linestyle)
        
# Cropland
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        (dd_velo[j].where((pft_ac>=70)&(dd_velo[j].lat>=40),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax44, label=labels[j], color=colors[j], ls=linestyle)
        
# Grassland
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        (dd_velo[j].where((pft_gr>=77)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax45, label=labels[j], color=colors[j], ls=linestyle)
                
# Tundra
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        (dd_velo[j].where((pft_tu>=80)&(dd_velo[j].lat>=0),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax46, label=labels[j], color=colors[j], ls=linestyle)
        
# Snow and ice (some other barren ground?)
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        (dd_velo[j].where((pft_is>=100),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax47, label=labels[j], color=colors[j], ls=linestyle)
            
# Ocean
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        (dd_velo[j].where((pft_oc>=100),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax48, label=labels[j], color=colors[j], ls=linestyle)
        
# Desert
for j, (name,linestyle) in enumerate(linestyles.items()):
    if j < len(dd_velo):
        (dd_velo[j].where((pft_de>=100),drop=True).mean(dim='lat').mean(dim='lon')).plot(ax=ax49, label=labels[j], color=colors[j], ls=linestyle)
        
       
ax41.scatter(dd_velo[0].time.data, hardacre_cf, label='Hardacre (2015)', color='red')
ax42.scatter(dd_velo[0].time.data, hardacre_df, label='Hardacre (2015)', color='red')
ax43.scatter(dd_velo[0].time.data, hardacre_tf, label='Hardacre (2015)', color='red')
ax44.scatter(dd_velo[0].time.data, hardacre_ac, label='Hardacre (2015)', color='red')
ax45.scatter(dd_velo[0].time.data, hardacre_gr, label='Hardacre (2015)', color='red')
ax46.scatter(dd_velo[0].time.data, hardacre_tu, label='Hardacre (2015)', color='red')
ax47.scatter(dd_velo[0].time.data, hardacre_is, label='Hardacre (2015)', color='red')
ax48.scatter(dd_velo[0].time.data, hardacre_oc, label='Hardacre (2015)', color='red')
ax49.scatter(dd_velo[0].time.data, hardacre_de, label='Hardacre (2015)', color='red')

ax41.set_title("Coniferous forest", loc='right', x=0.99, y=0.85)
ax42.set_title("Deciduous forest", loc='right', x=0.99, y=0.85)
ax43.set_title("Tropical forest", loc='right', x=0.99, y=0.85)
ax44.set_title("Cropland", loc='right', x=0.99, y=0.85)
ax45.set_title("Grassland", loc='right', x=0.99, y=0.85)
ax46.set_title("Tundra", loc='right', x=0.99, y=0.85)
ax47.set_title("Ice/Snow", loc='right', x=0.99, y=0.85)
ax48.set_title("Ocean", loc='right', x=0.99, y=0.85)
ax49.set_title("Desert", loc='right', x=0.99, y=0.85)

for ax in fig4.axes:
    ticks = ax.get_xticks()
    start = ticks[0]-(ticks[1]-ticks[0])/2.
    end = ticks[-1]
    ax.set_title("")
    ax.set_xticks(np.linspace(start, end, 12))
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_ylim(0,1.5)
for ax in fig4.axes[-3:]:
    ax.set_ylim(0,0.3)
    
ax44.set_ylabel("$v^{O_3}_{DD}$ (%s)" % ("$cm\,s^{-1}$"))
#ax48.set_xlabel("Time (months)")
ax47.set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
ax48.set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
ax49.set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
ax48.legend(bbox_to_anchor=(0.45, -0.55), loc=8, borderaxespad=0.,ncol=5)
#if b_nh:
#    ax42.set_title("Northern hemisphere", weight='bold')
#else:
#    ax42.set_title("Southern hemisphere", weight='bold')

if b_norm:
    # Surface ozone
    fig5 = plt.figure(5, figsize=(16,9))
    ax5 = plt.subplot()
    fig5.canvas.set_window_title("final-norm_diff_surfaceozone")

    for i, (name,linestyle) in enumerate(linestyles.items()):
        if i < len(dd_velo_zonal):
            ((ozone_raw_data[i]-ozone_raw_data[0]).mean(dim='lon').mean(dim='time')/ozone_raw_data[0].mean(dim='lon').mean(dim='time')).plot(label=labels[i], color=colors[i], ls=linestyle)
      
    ax5.set_title('')
    ax5.set_xlabel("Latitude (deg)")
    
    ax5.set_ylabel("$\Delta O_3 / O_3^{ref}$")
    #ax5.scatter(hardacre_lat, hardacre_ddvel/hardacre_ddvel.sum(), label='mmm Hardacre (2015)', color='red')
    ax5.set_ylim(-0.05, 1)
    ax5.legend(ncol=3)

    # Surface ozone split by month
    fig6 = plt.figure(6, figsize=(16,9))
    fig6.canvas.set_window_title("final-norm_diff_surfaceozone_monthly")
    rows, col = 4, 3
    for i in np.arange(1,rows*col+1):
        plt.subplot(rows, col, i)
        ax =  plt.gca()
        for j, (name,linestyle) in enumerate(linestyles.items()):
            if j < len(ozone_raw_data):
                ((ozone_raw_data[j]-ozone_raw_data[0]).mean(dim='lon').isel(time=i-1)/ozone_raw_data[0].mean(dim='lon').isel(time=i-1)).plot(label=labels[j], color=colors[j], ls=linestyle)
                               
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_ylim(-0.1,2)
        ax.set_title("")
        ax.set_title(get_month_name(i), x=0.05, y=0.8, loc='left')
    fig6.axes[-2].set_xlabel("Latitude (deg)")
    fig6.axes[6].set_ylabel("$\Delta O_3 / O_3^{ref}$", y=1)
    
# Global maps
fig7 = plt.figure(7, figsize=(9,16))
fig7.canvas.set_window_title("final-ozone_diff_maps")
ax71 = plt.subplot(311, projection=cp.crs.PlateCarree())
ax72 = plt.subplot(312, projection=cp.crs.PlateCarree())
ax73 = plt.subplot(313, projection=cp.crs.PlateCarree())
levels = np.arange(-1,1.1,0.1)
((ozone_raw_data[-3]-ozone_raw_data[0])/ozone_raw_data[0]).mean(dim='time').plot(ax=ax71, transform=cp.crs.PlateCarree(), vmin=-1, vmax=1, cmap=plt.cm.RdYlBu_r)
((dd_velo[-3]-dd_velo[0])/dd_velo[0]).mean(dim='time').plot(ax=ax72, transform=cp.crs.PlateCarree(), vmin=-1, vmax=1, cmap=plt.cm.RdYlBu_r)
((dry_dep_raw_data[-3].sum(dim='time')-dry_dep_raw_data[0].sum(dim='time'))/dry_dep_raw_data[0].sum(dim='time')).plot(ax=ax73, transform=cp.crs.PlateCarree(), vmin=-1, vmax=1, cmap=plt.cm.RdYlBu_r)

for ax in fig7.axes[:-3]:
    ax.set_global()
    ax.set_aspect('auto')
    ax.set_title('')
    ax.set_xticks(np.arange(-180,181,60), crs=cp.crs.PlateCarree())
    ax.set_yticks(np.arange(-90,91,30), crs=cp.crs.PlateCarree())
    ax.set_xlabel("")
    ax.set_ylabel("")
    #lon_formatter = LongitudeFormatter(number_format='.1f',
    #                                   degree_symbol='',
    #                                   dateline_direction_label=True)
    #lat_formatter = LatitudeFormatter(number_format='.1f',
    #                                  degree_symbol='')
    #ax.xaxis.set_major_formatter(lon_formatter)
    #ax.yaxis.set_major_formatter(lat_formatter)
    ax.coastlines()
    
fig7.axes[0].set_title("(a)")
fig7.axes[1].set_title("(b)")
fig7.axes[2].set_title("(c)")
fig7.axes[2].set_xlabel("Latitude (deg)")
fig7.axes[1].set_ylabel("Longitude (deg)")
fig7.axes[-3].set_ylabel("$<\Delta O_3 / O_3^{ref}>$")
fig7.axes[-2].set_ylabel("$<\Delta v^{O_3}_{DD} / v^{O_3}_{DDref}>$")
fig7.axes[-1].set_ylabel("$<\Delta \Sigma O_3^{DD} / \Sigma O_3^{DDref}>$")
# Show it
plt.show(block=False)

# Print total annual ozone dry deposition
outfile = open("ozone_dry_deposition.txt", 'w')
# Write file line by line
outfile.write("Ocean: NH - SH - Global (Tg a-1)\n")
for each,label in zip(dry_dep_raw_data,labels):
    nh = (each.where((pft_oc>=98)&(each.lat>=0),drop=True)).sum()*1e-9
    sh = (each.where((pft_oc>=98)&(each.lat<0),drop=True)).sum()*1e-9
    gl = each.where(pft_oc>=98).sum()*1e-9
    outfile.write("%3.1f %3.1f %3.1f %s\n" % (nh, sh, gl, label))

outfile.write("Ice: NH - SH - Global (Tg a-1)\n")
for each,label in zip(dry_dep_raw_data,labels):
    nh = (each.where((pft_is>=98)&(each.lat>=0),drop=True)).sum()*1e-9
    sh = (each.where((pft_is>=98)&(each.lat<0),drop=True)).sum()*1e-9
    gl = each.where(pft_is>=98).sum()*1e-9
    outfile.write("%3.1f %3.1f %3.1f %s\n" % (nh, sh, gl, label))

outfile.write("Land: NH - SH - Global (Tg a-1\n)")
for each,label in zip(dry_dep_raw_data,labels):
    nh = (each.where((pft_cf+pft_df+pft_tf+pft_ac+pft_gr+pft_tu+pft_de>=98)&(each.lat>=0),drop=True)).sum()*1e-9
    sh = (each.where((pft_cf+pft_df+pft_tf+pft_ac+pft_gr+pft_tu+pft_de>=98)&(each.lat<0),drop=True)).sum()*1e-9
    gl = each.where(pft_cf+pft_df+pft_tf+pft_ac+pft_gr+pft_tu+pft_de>=98).sum()*1e-9
    outfile.write("%3.1f %3.1f %3.1f %s\n" % (nh, sh, gl, label))
    
outfile.write("Ocean, ice, land: NH - SH - Global (Tg a-1)\n")
for each,label in zip(dry_dep_raw_data,labels):
    nh = (each.where(each.lat>=0,drop=True)).sum()*1e-9
    sh = (each.where(each.lat<0,drop=True)).sum()*1e-9
    gl = each.sum()*1e-9
    outfile.write("%3.1f %3.1f %3.1f %s\n" % (nh, sh, gl, label))

outfile.write("\n\n")
outfile.close()

print("Written global sums to %s" % ("ozone_dry_deposition.txt"))
