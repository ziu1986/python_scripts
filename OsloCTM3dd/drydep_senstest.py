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
# Data source
b_original_resolution = True
if b_original_resolution:
    scav_dir = "scavenging_monthly/"
    mm_dir = "VMR/"
else:
    scav_dir = "scavenging_monthly/regrid_hardacre/"
    mm_dir = "VMR/regrid_hardacre/"
       
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

    
labels = ('OsloCTM3: Wesely type',
          'OsloCTM3: EMEP_full','OsloCTM3: EMEP_offLight','OsloCTM3: EMEP_offPhen','OsloCTM3: EMEP_SWVL4',
          'OsloCTM3: EMEP_ppgs','OsloCTM3: EMEP_ppgssh',
          'OsloCTM3: EMEP_ppgssh_ice',
          'OsloCTM3: EMEP_ppgs_2005')
colors = np.concatenate((('black',),
                         ('blue',),
                         np.repeat('cornflowerblue',3),
                         np.repeat('purple',2),
                         ('grey',),
                         ('orange',)) )
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

# Read the data
try:
    data
except NameError:
    ozone_raw_data = []
    for iexp in experiment:
        subdir = data_dir+iexp+mm_dir+'*.nc'
        print("Reading from path %s" % (os.path.abspath(subdir)))
        data_list = []
        # Open dataset
        for file in sorted(glob.glob(subdir)):
            print("Reading %s" % (os.path.basename(file)))
            data = xr.open_dataset(file)
            # Define new time coordinates and drop the old one
            data.coords['time'] = ([dt.datetime.strptime(os.path.basename(file)[10:18], '%Y%M%d'),])
            if (data.lat.ndim > 1):
                print("Changing coordinates...")
                data.coords['x'] = data.lat.lon[0].data
                data.coords['y'] = data.lat[:,0].data
                data = data.drop(('lon','lat'))
                data = data.rename({'x':'lon','y':'lat'})
            data_list.append(data['O3'])
        # Concatenating the list
        ozone_raw_data.append(xr.concat(data_list, dim='time'))

    #dry_dep_raw_data = []
    #for iexp in experiment:
    #    subdir = data_dir+iexp+scav_dir+'*.nc'
    #    print("Reading from path %s" % (os.path.abspath(subdir)))
    #    data_list = []
        # Open dataset
    #    for file in sorted(glob.glob(subdir)):
    #        print("Reading %s" % (os.path.basename(file)))
    #        data = xr.open_dataset(file)
            # Define new time coordinates and drop the old one
    #        data['time'].reset_coords(drop=True)
    #        data.coords['time'] = ([dt.datetime(data['YEAR'], data['MONTH'], 1),])
    #        if (data.lat.ndim > 1):
    #            print("Changing coordinates...")
    #            data.coords['x'] = data.lat.lon[0].data
    #            data.coords['y'] = data.lat[:,0].data
    #            data = data.drop(('lon','lat'))
    #            data = data.rename({'x':'lon','y':'lat'})
    #        data = data['dry_O3']/mega
    #        data.attrs['unit'] = 'Gg'
    #        data_list.append(data)
        # Concatenating the list
    #    dry_dep_raw_data.append(xr.concat(data_list, dim='time'))
        
    # Compute dry deposition velocity
    #sim = np.array([seconds_in_month(i, np.unique(dry_dep_raw_data[0].time.dt.year)) for i in range(1,13)])
    #dim = sim/(24*60**2)
    #layer_thickness = 8*100 #[1m -> 100cm]
    #dd_velo = []
    #for i in np.arange(len(dry_dep_raw_data)):
    #    tmp_velo = dry_dep_raw_data[i]/ozone_raw_data[i]
    #    tmp_velo_list = []
    #    for imonth in range(0,12):
    #        tmp_velo_list.append(tmp_velo.where(tmp_velo.time.dt.month==imonth+1, drop=True)/sim[imonth]*layer_thickness)
    #    dd_velo.append(xr.concat(tmp_velo_list, dim='time'))

    #dd_velo_zonal = [each.mean(dim='lon') for each in dd_velo]
    #for each in dd_velo_zonal:
    #    each.attrs['unit'] = 'cms-1'
plt.close('all')
surface_ozone = ozone_raw_data
test = xr.concat(surface_ozone[1:], dim='test')
print("Global - SH - NH")
print((test.mean(dim='test')-surface_ozone[0]).mean(dim='lon').mean(), (test.mean(dim='test')-surface_ozone[0]).mean(dim='lon').sel(lat=slice(-90,0)).mean(), (test.mean(dim='test')-surface_ozone[0]).mean(dim='lon').sel(lat=slice(0,90)).mean())
# Surface ozone
surface_ozone = [each.isel(lev=0) for each in ozone_raw_data]
fig1 = plt.figure(1, figsize=(16,9))
ax1 = plt.subplot()
fig1.canvas.set_window_title("final-vmr_surfaceozone")

for i, (name,linestyle) in enumerate(linestyles.items()):
    if i < len(surface_ozone):
        (((surface_ozone[i])*1e9).mean(dim='lon').mean(dim='time')).plot(label=labels[i], color=colors[i], ls=linestyle)
        
ax1.set_title('')
ax1.set_xlabel("Latitude (deg)")

ax1.set_ylabel("$O_3$ (ppb)")
#ax1.scatter(hardacre_lat, hardacre_ddvel/hardacre_ddvel.sum(), label='mmm Hardacre (2015)', color='red')
ax1.set_ylim(0, 60)
ax1.legend(ncol=3)

# Surface ozone split by month
fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("final-vmr_surfaceozone_monthly")
rows, col = 4, 3
for i in np.arange(1,rows*col+1):
    plt.subplot(rows, col, i)
    ax =  plt.gca()
    for j, (name,linestyle) in enumerate(linestyles.items()):
        if j < len(surface_ozone):
            (((surface_ozone[j])*1e9).mean(dim='lon').isel(time=i-1)).plot(label=labels[j], color=colors[j], ls=linestyle)
                
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_ylim(0,60)
    ax.set_title("")
    ax.set_title(get_month_name(i), x=0.05, y=0.8, loc='left')
fig2.axes[-2].set_xlabel("Latitude (deg)")
fig2.axes[6].set_ylabel("$O_3$ (ppb)", y=1)   

# Surface ozone relative
fig3 = plt.figure(3, figsize=(16,9))
ax3 = plt.subplot()
fig3.canvas.set_window_title("final-diff_vmr_surfaceozone")

for i, (name,linestyle) in enumerate(linestyles.items()):
    if i < len(surface_ozone):
        (((surface_ozone[i]-surface_ozone[0])*1e9).mean(dim='lon').mean(dim='time')).plot(label=labels[i], color=colors[i], ls=linestyle)
        
ax3.set_title('')
ax3.set_xlabel("Latitude (deg)")

ax3.set_ylabel("$\Delta O_3$ (ppb)")
ax3.legend(ncol=3)

# Surface ozone split by month
fig4 = plt.figure(4, figsize=(16,9))
fig4.canvas.set_window_title("final-diff_vmr_surfaceozone_monthly")
rows, col = 4, 3
for i in np.arange(1,rows*col+1):
    plt.subplot(rows, col, i)
    ax =  plt.gca()
    for j, (name,linestyle) in enumerate(linestyles.items()):
        if j < len(surface_ozone):
            (((surface_ozone[j]-surface_ozone[0])*1e9).mean(dim='lon').isel(time=i-1)).plot(label=labels[j], color=colors[j], ls=linestyle)
                
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_ylim(-0.1,30)
    ax.set_title("")
    ax.set_title(get_month_name(i), x=0.05, y=0.8, loc='left')
fig4.axes[-2].set_xlabel("Latitude (deg)")
fig4.axes[6].set_ylabel("$\Delta O_3$ (ppb)", y=1)

# Zonal level plot
fig5 = plt.figure(5,figsize=(16,9))
fig5.canvas.set_window_title("final-diff_vmr_zonal_ozone_profile")
ax51 = plt.subplot(231)
ax52 = plt.subplot(232)
ax53 = plt.subplot(233)
#ax54 = plt.subplot(234)
ax55 = plt.subplot(235)
ax56 = plt.subplot(236)

cbaxes = fig5.add_axes([0.325, 0.1, 0.0075, 0.3875])

levels = np.arange(0,80,5)
levels = np.append(levels,(np.arange(80,180,20)))
#levels = np.append(levels,(np.arange(180,1090,100)))

(ozone_raw_data[0]*1e9).mean(dim='lon').mean('time').plot.contourf(ax=ax51, levels=levels, add_colorbar=False)
(ozone_raw_data[-3]*1e9).mean(dim='lon').mean('time').plot.contourf(ax=ax52, levels=levels, add_colorbar=False)
(ozone_raw_data[-2]*1e9).mean(dim='lon').mean('time').plot.contourf(ax=ax53, levels=levels, cbar_kwargs={'label':'O$_3$ (ppb)','fraction':0.046, 'pad':0.04,'aspect':30})

((ozone_raw_data[-3]-ozone_raw_data[0])*1e9).mean(dim='lon').mean('time').plot.contourf(ax=ax55, levels=np.arange(-50,51,2), cmap=plt.cm.RdBu_r, add_colorbar=False)
((ozone_raw_data[-2]-ozone_raw_data[0])*1e9).mean(dim='lon').mean('time').plot.contourf(ax=ax56, levels=np.arange(-50,51,2), cmap=plt.cm.RdBu_r, cbar_kwargs={'label':'$\Delta O_3$ (ppb)','fraction':0.046, 'pad':0.04,'aspect':30})

(np.fabs(ozone_raw_data[-3]-ozone_raw_data[0])/ozone_raw_data[0]).mean(dim='lon').mean('time').plot.contourf(ax=ax55, levels=np.arange(0,1.1,0.25), hatches=[ None, '...','---','///','**' ], colors='none', add_colorbar=False)
(np.fabs(ozone_raw_data[-2]-ozone_raw_data[0])/ozone_raw_data[0]).mean(dim='lon').mean('time').plot.contourf(ax=ax56, levels=np.arange(0,1.1,0.25), hatches=[ None, '...', '---', '///' ,'**'], colors='none', add_colorbar=True, cbar_kwargs={'label':'$\Delta O_3 / O_3^{ref}$','cax':cbaxes})

#np.arange(0,1.1,0.25)

axes = (ax51,ax52,ax53,ax55,ax56)
for ax in axes:
    ax.invert_yaxis()
    ax.set_ylim(USstdP(0)/hecto, USstdP(15*kilo)/hecto)
    ax.set_xlabel("Latitude (deg)")
    ax.set_ylabel("")
    
ax51.set_ylabel("Pressure (hPa)")
ax55.set_ylabel("Pressure (hPa)")
ax52.set_yticklabels("")
ax53.set_yticklabels("")
ax56.set_yticklabels("")
#ax51.set_title("%s" % (labels[0]))
#ax52.set_title("%s" % (labels[-3]))
#ax53.set_title("%s" % (labels[-2]))
ax51.set_title("%s" % ("(a)"))
ax52.set_title("%s" % ("(b)"))
ax53.set_title("%s" % ("(c)"))
ax55.set_title("%s" % ("(d)"))
ax56.set_title("%s" % ("(e)"))
ax52.set_xlabel("")
ax53.set_xlabel("")

#cbaxes.set_xlabel('$\Delta O_3 / O_3^{ref}$')
cbaxes.yaxis.set_label_position('left')
cbaxes.yaxis.set_ticks_position('left')

# Ozone relative
#fig6 = plt.figure(6, figsize=(16,9))
#fig6.canvas.set_window_title("final-diff_ozone")
#ax61 = plt.subplot()
#ax62 = plt.subplot(312)
#ax63 = plt.subplot(313)

#for i, (name,linestyle) in enumerate(linestyles.items()):
#    if i < len(surface_ozone):
#        (((surface_ozone[i]-surface_ozone[0])*1e9).mean(dim='lon').mean(dim='time')).plot(ax=ax61, label=labels[i], color=colors[i], ls=linestyle)
        
#for i, (name,linestyle) in enumerate(linestyles.items()):
#    if i < len(surface_ozone):
#        ((dd_velo[i]-dd_velo[0]).mean(dim='lon').mean(dim='time')).isel(lev=0).plot(ax=ax62, label=labels[i], color=colors[i], ls=linestyle)

#for i, (name,linestyle) in enumerate(linestyles.items()):
#    if i < len(surface_ozone):
#        ((dry_dep_raw_data[i]-dry_dep_raw_data[0]).mean(dim='lon').mean(dim='time')).plot(ax=ax63, label=labels[i], color=colors[i], ls=linestyle)      

#for ax in fig6.axes:
#    ax.set_title('')
#    ax.set_xlabel('')
#ax61.set_ylabel("$\Delta O_3$ (ppb)")
#ax62.set_ylabel("$\Delta v_{O_3}^{DD} (cm s^{-1})$")
#ax63.set_ylabel("$\Delta O_3^{DD}$ (Gg)")   
#ax61.set_xlabel("Latitude (deg)")


#ax61.legend(ncol=3)
# Show it
plt.show(block=False)
