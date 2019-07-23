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

b_norm = False
species = 'O_3'
names = {'O_3':'ozone', 'HNO_3':'nitric_acid', 'SO_2':'sulfur_dioxid', 'NH_3':'ammonia', 'SO_4':'sulfate'}
title_name = names[species]
# Data source
b_original_resolution = False
b_ozone = True
if b_original_resolution:
    scav_dir = "scavenging_monthly/"
    mm_dir = "monthly_means/"
else:
    scav_dir = "scavenging_monthly/hardacre_grid/"
    mm_dir = "monthly_means/hardacre_grid/"
experiment = ('C3RUN_default/',
              'C3RUN_mOSaic/',
              'C3RUN_mOSaic_offLight/',
              'C3RUN_mOSaic_offPhen/',
              'C3RUN_mOSaic_SWVL1/',
              'C3RUN_mOSaic_ice/',
              'C3RUN_mOSaic_desert/',
              'C3RUN_mOSaic_emis2014/',
              'C3RUN_mOSaic_hough/'
              #'C3RUN_oDD/',
              #'C3RUN_emep_full/',
              #'C3RUN_emep_offLight/',
              #'C3RUN_emep_offPhen/',
              #'C3RUN_emep_SWVL4/',
              #'C3RUN_emep_ppgs/',
              #'C3RUN_emep_ppgssh/',
              #'C3RUN_emep_ppgssh_ice/',
              #'C3RUN_emep_ppgs_2005/'
)

data_dir = os.environ['DATA']+'/astra_data/ctm_results/' 
    
labels = ('Wesely_type',
          'mOSaic',
          'mOSaic_offLight','mOSaic_offPhen','mOSaic_SWVL1',
          'mOSaic_ice',
          'mOSaic_desert',
          'mOSaic_emis2014',
          'mOSaic_hough'
    #'OsloCTM3: Wesely type',
          #'OsloCTM3: EMEP_swgd',
          #'OsloCTM3: EMEP_full',
          #'OsloCTM3: EMEP_offLight','OsloCTM3: EMEP_offPhen','OsloCTM3: EMEP_SWVL4',
          #'OsloCTM3: EMEP_ppgs','OsloCTM3: EMEP_ppgssh',
          #'OsloCTM3: EMEP_ppgssh_ice',
          #'OsloCTM3: EMEP_ppgs_2005'
)
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
    raw_data = []
    for iexp in experiment:
        subdir = data_dir+iexp+scav_dir+'sum*.nc'
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
            data_list.append(data)
        # Concatenating the list
        raw_data.append(xr.concat(data_list, dim='time'))
        
# Extract drydeposition and rescale units
species_data = [date['dry_%s' % (species.replace('_',''))]/mega for date in raw_data]
for date in species_data:
    date.attrs['unit'] = 'Gg'
    # Monthly zonal total and rescale units
try:
    ozone_zonal = [date.sum(dim='lon')/kilo for date in species_data] #sum(dim='time').
except ValueError:
    ozone_zonal = [date.sum(dim='x')/kilo for date in species_data]
for date in ozone_zonal:
    date.attrs['unit'] = 'Tg'
# Extracting some general information
# WARNING: cdo has summed these also up -> devide by number of days
gridarea = raw_data[0]['gridarea'].isel(time=0)/31. 
tot_surface_area = gridarea.sum()
gridarea_faction = gridarea/tot_surface_area
molarweight = get_molarweight(raw_data[0].isel(time=0))/31.

weights = np.cos(raw_data[0].lat*np.pi/180) #1

# Clean up
plt.close('all')

# Extracted numbers from Hardacre et al. (2015) 
hardacre_data = np.load("Hardacre2015.npz")

# Plot zonal means
plt.close('all')
# General styles
linestyles = OrderedDict(
    [('solid',               (0, ())),
     #('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 2))),
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

hardacre_lat = np.arange(-88.5,89,3, dtype=np.float64)
new_hardacre_lat = (hardacre_lat[4:-4])[::-1]
# Total ozone dry deposition over latitudes
fig1 = plt.figure(1, figsize=(16,9))
if b_norm:
    fig1.canvas.set_window_title("final-norm_total_%s_drydep" % (title_name))
else:
    fig1.canvas.set_window_title("final-total_%s_drydep" % (title_name))
for i, (name,linestyle) in enumerate(linestyles.items()):
    if i < len(ozone_zonal):
        if b_norm:
            (ozone_zonal[i].sum(dim='time')/(ozone_zonal[i]).reduce(np.sum)).plot(label=labels[i], color=colors[i], ls=linestyle)
        else:
            (ozone_zonal[i]).sum(dim='time').plot(label=labels[i], color=colors[i], ls=linestyle)
        
ax1 = plt.gca()
ax1.set_xlabel("Latitude (deg)")
if b_norm:
    ax1.set_ylabel("$%s^{DD} / \Sigma %s^{DD}$" % (species, species))
    if b_ozone:
        ax1.scatter(new_hardacre_lat, hardacre_data['dd_tot'][1]/hardacre_data['dd_tot'][1].sum(), label='mmm Hardacre (2015)', color='red')
    if species=="O_3":
        ax1.set_ylim(0,0.05)
    elif species=="HNO_3":
        ax1.set_ylim(0,0.06)
    elif species=="SO_2":
        ax1.set_ylim(0,0.08)
    elif species=="NH_3":
        ax1.set_ylim(0,0.07)
    elif species=="SO_4":
        ax1.set_ylim(0,0.05)
else:
    ax1.set_ylabel("$%s^{DD}$ (%s lat$^{-1}$ a$^{-1}$)" % (species, ozone_zonal[0].attrs['unit']))
    if b_ozone:
        ax1.scatter(new_hardacre_lat, hardacre_data['dd_tot'][1], label='mmm Hardacre (2015)', color='red')
        #ax1.fill_between(new_hardacre_lat, hardacre_data['dd_tot'][1]*4/5, hardacre_data['dd_tot'][1]*6/5, color='red', alpha=0.15)
    if species=="O_3":
        ax1.set_ylim(0,40)
    elif species=="HNO_3":
        ax1.set_ylim(0,5)
    elif species=="SO_2":
        ax1.set_ylim(0,7)
    elif species=="NH_3":
        ax1.set_ylim(0,3)
    elif species=="SO_4":
        ax1.set_ylim(0,1)
ax1.legend(ncol=2)
# Total ozone dry deposition over latitudes split by month
fig2 = plt.figure(2, figsize=(16,9))
if b_norm:
    fig2.canvas.set_window_title("final-norm_total_%s_drydep_monthly" % (title_name))
else:
    fig2.canvas.set_window_title("final-total_%s_drydep_monthly" % (title_name))

rows, col = 4, 3
for i in np.arange(1,rows*col+1):
    plt.subplot(rows, col, i)
    ax =  plt.gca()
    for j, (name,linestyle) in enumerate(linestyles.items()):
        if j < len(ozone_zonal):
            if b_norm:
                ((ozone_zonal[j]).isel(time=i-1)/(ozone_zonal[j].isel(time=i-1).reduce(np.sum))).plot(label=labels[j], color=colors[j], ls=linestyle)
            else:
                (ozone_zonal[j]).isel(time=i-1).plot(label=labels[j], color=colors[j], ls=linestyle)
                    
    if i==2 and b_ozone:
        if b_norm:
            ax.scatter(new_hardacre_lat, hardacre_data['dd_feb'][1]/hardacre_data['dd_feb'][1].sum(), label='Hardacre (2015)', color='red')
        else:
            ax.scatter(new_hardacre_lat, hardacre_data['dd_feb'][1], label='Hardacre (2015)', color='red')
            #ax.fill_between(new_hardacre_lat, hardacre_data['dd_feb'][1]*4/5, hardacre_data['dd_feb'][1]*6/5, color='red', alpha=0.15)
    elif i==8 and b_ozone:
        if b_norm:
            ax.scatter(new_hardacre_lat, hardacre_data['dd_aug'][1]/hardacre_data['dd_aug'][1].sum(), label='Hardacre (2015)', color='red')
            
        else:
            ax.scatter(new_hardacre_lat, hardacre_data['dd_aug'][1], label='Hardacre (2015)', color='red')
            #ax.fill_between(new_hardacre_lat, hardacre_data['dd_aug'][1]*4/5, hardacre_data['dd_aug'][1]*6/5, color='red', alpha=0.15)
    ax.set_xlabel("")
    ax.set_ylabel("")
    if b_norm:
        if species=="O_3":
            ax.set_ylim(0,0.05)
        elif species=="HNO_3":
            ax.set_ylim(0,0.08)
        elif species=="SO_2":
            ax.set_ylim(0,0.08)
        elif species=="NH_3":
            ax.set_ylim(0,0.08)
        elif species=="SO_4":
            ax.set_ylim(0,0.06)   
    else:
        if species=="O_3":
            ax.set_ylim(0,4)
        elif species=="HNO_3":
            ax.set_ylim(0,0.5)
        elif species=="SO_2":
            ax.set_ylim(0,0.7)
        elif species=="NH_3":
            ax.set_ylim(0,0.6)
        elif species=="SO_4":
            ax.set_ylim(0,0.1)
    ax.set_title("")
    ax.set_title(get_month_name(i), x=0.05, y=0.8, loc='left')
fig2.axes[-2].set_xlabel("Latitude (deg)")
if b_norm:
    fig2.axes[6].set_ylabel("$%s^{DD}/\Sigma %s^{DD}$" % (species, species), y=1)
else:
    fig2.axes[6].set_ylabel("$%s^{DD}$ (%s lat$^{-1}$ month$^{-1}$)" % (species, ozone_zonal[0].attrs['unit']), y=1)
        
# Total ozone dry deposition integrated over hemispheres over months
fig3 = plt.figure(3, figsize=(16,9))
if b_norm:
    fig3.canvas.set_window_title("final-norm_total_%s_drydep_hem" % (title_name))
else:
    fig3.canvas.set_window_title("final-total_%s_drydep_hem" % (title_name))
ax31 = plt.subplot(131)
for j, (name,linestyle) in enumerate(linestyles.items()):
            if j < len(ozone_zonal):
                if b_norm:
                    ((ozone_zonal[j]).sel(lat=slice(30,90)).sum(dim='lat')/(ozone_zonal[j]).sel(lat=slice(30,90)).sum()).plot(label=labels[j], color=colors[j], ls=linestyle)
                else:
                    (ozone_zonal[j]).sel(lat=slice(30,90)).sum(dim='lat').plot(label=labels[j], color=colors[j], ls=linestyle)
                    print(labels[j], ((ozone_zonal[j]-ozone_zonal[0])/ozone_zonal[0]).sel(lat=slice(30,90)).sum(dim='lat').mean(dim='time'))
ax31.set_title("NH (30N-90N)")
if b_norm:
    ax31.set_ylabel("$%s^{DD}/\Sigma %s^{DD}$" % (species, species))
else:
    ax31.set_ylabel("$%s^{DD}$ (%s)" % (species, ozone_zonal[0].attrs['unit']))

ax32 = plt.subplot(132)
for j, (name,linestyle) in enumerate(linestyles.items()):
            if j < len(ozone_zonal):
                if b_norm:
                    ((ozone_zonal[j]).sel(lat=slice(-30,30)).sum(dim='lat')/(ozone_zonal[j]).sel(lat=slice(-30,30)).sum()).plot(label=labels[j], color=colors[j], ls=linestyle)
                else:
                    (ozone_zonal[j]).sel(lat=slice(-30,30)).sum(dim='lat').plot(label=labels[j], color=colors[j], ls=linestyle)
                    print(labels[j], ((ozone_zonal[j]-ozone_zonal[0])/ozone_zonal[0]).sel(lat=slice(-30,30)).sum(dim='lat').mean(dim='time'))
ax32.set_title("TR (30S-30N)")
ax32.set_ylabel("")

ax33 = plt.subplot(133)
for j, (name,linestyle) in enumerate(linestyles.items()):
            if j < len(ozone_zonal):
                if b_norm:
                    ((ozone_zonal[j]).sel(lat=slice(-90,-30)).sum(dim='lat')/(ozone_zonal[j]).sel(lat=slice(-90,-30)).sum()).plot(label=labels[j], color=colors[j], ls=linestyle)
                else:
                    (ozone_zonal[j]).sel(lat=slice(-90,-30)).sum(dim='lat').plot(label=labels[j], color=colors[j], ls=linestyle)
                    print(labels[j], ((ozone_zonal[j]-ozone_zonal[0])/ozone_zonal[0]).sel(lat=slice(-90,-30)).sum(dim='lat').mean(dim='time'))
ax33.set_title("SH (90S-30S)")
ax33.set_ylabel("")

for ax, hardacre in zip(fig3.axes, hardacre_data['dd_hem']):
    ax.set_xlabel("")
    ticks = ax.get_xticks()
    start = ticks[0]-(ticks[1]-ticks[0])/2.
    end = ticks[-1]
    ax.set_xticks(np.linspace(start, end, 12))
    ax.set_xticklabels([get_month_name(imonth, length=1) for imonth in range(1,13)])
    if b_norm:
        ax.set_ylim(0,0.2)
        if b_ozone:
            ax.scatter(ax.get_xticks(), hardacre/hardacre.sum(), label='Hardacre (2015)', color='red')
    else:
        if species=="O_3":
            ax.set_ylim(0,60)
        elif species=="HNO_3":
            ax.set_ylim(0,6)
        elif species=="SO_2":
            ax.set_ylim(0,6)
        elif species=="NH_3":
            ax.set_ylim(0,4)
        elif species=="SO_4":
            ax.set_ylim(0,1)
        if b_ozone:
            ax.scatter(ax.get_xticks(), hardacre, label='Hardacre (2015)', color='red')
        
#ax32.legend(bbox_to_anchor=(0.45, -0.55), loc=8, borderaxespad=0.,ncol=5)    
#ax32.set_xlabel("Time (months)")
ax32.legend(bbox_to_anchor=(0.45, -0.12),loc=8, borderaxespad=0., ncol=5)

if b_norm:
    # Surface ozone
    fig5 = plt.figure(5, figsize=(16,9))
    ax5 = plt.subplot()
    fig5.canvas.set_window_title("final-norm_%s_diff_drydep" % (title_name))

    for i, (name,linestyle) in enumerate(linestyles.items()):
        if i < len(species_data):
            ((species_data[i]-species_data[1]).sum(dim='lon').sum(dim='time')/species_data[1].sum(dim='lon').sum(dim='time')).plot(label=labels[i], color=colors[i], ls=linestyle)
      
    ax5.set_title('')
    ax5.set_xlabel("Latitude (deg)")
    
    ax5.set_ylabel("$\Delta %s / %s^{ref}$" % (species, species))
    #ax5.scatter(hardacre_lat, hardacre_ddvel/hardacre_ddvel.sum(), label='mmm Hardacre (2015)', color='red')
    if species=='O_3':
        #ax5.set_ylim(-1, 0.2)
        ax5.set_ylim(-1, 1)
    elif species=="HNO_3":
        ax5.set_ylim(-0.2, 20)
    elif species=="SO_2":
        ax5.set_ylim(-0.2, 1.2)
    elif species=="NH_3":
        ax5.set_ylim(-1.8, 5)
    elif species=="SO_4":
        ax5.set_ylim(-0.6, 0.2)
        
    ax5.legend(ncol=3)

    # Surface ozone split by month
    fig6 = plt.figure(6, figsize=(16,9))
    fig6.canvas.set_window_title("final-norm_%s_diff_drydep_monthly" % (title_name))
    rows, col = 4, 3
    for i in np.arange(1,rows*col+1):
        plt.subplot(rows, col, i)
        ax =  plt.gca()
        for j, (name,linestyle) in enumerate(linestyles.items()):
            if j < len(species_data):
                ((species_data[j]-species_data[0]).sum(dim='lon').isel(time=i-1)/species_data[0].sum(dim='lon').isel(time=i-1)).plot(label=labels[j], color=colors[j], ls=linestyle)
                               
        ax.set_xlabel("")
        ax.set_ylabel("")
        if species=='O_3':
            ax.set_ylim(-1,0.2)
        elif species=="HNO_3":
            ax.set_ylim(-0.2, 40)
        elif species=="SO_2":
            ax.set_ylim(-0.5, 2.5)
        elif species=="NH_3":
            ax.set_ylim(-2., 5) #1500.
        elif species=="SO_4":
            ax.set_ylim(-1, 0.2)
        ax.set_title("")
        ax.set_title(get_month_name(i), x=0.05, y=0.8, loc='left')
    fig6.axes[-2].set_xlabel("Latitude (deg)")
    fig6.axes[6].set_ylabel("$\Delta %s / %s^{ref}$" % (species, species), y=1)

# Show it
plt.show(block=False)
