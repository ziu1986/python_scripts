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

# Hardacre
ulborg_mmm = 20/(403.13-311.92)*(np.array((317.46,319.40,321.76,326.06,330.63,333.13,332.29,328.69,324.95,321.62,319.96,317.05))-311.92)
ulborg_obs = 20/(403.13-311.92)*(np.array((337.56,351.98,346.16,355.45,361.68,368.75,363.76,360.44,357.66,354.06,342.55,334.51))-311.92)
castel_mmm = 20/(403.13-311.92)*(np.array((320.05,323.84,329.94,334.24,337.56,341.86,340.34,337.42,332.43,327.30,321.62,319.68))-311.92)
castel_obs = 20/(403.13-311.92)*(np.array((317.46,319.26,324.39,330.49,325.50,328.14,322.87,320.93,319.40,322.45,316.07,313.58))-311.92)
auchen_mmm = 20/(303.46-212.10)*(np.array((220.98,223.33,226.80,231.23,234.98,235.81,233.04,229.47,227.91,225.97,224.03,220.98))-212.10)
auchen_obs = 20/(303.46-212.10)*(np.array((230.40,238.86,240.52,250.09,248.84,241.77,241.77,243.29,234.42,239.83,232.62,232.20))-212.10)
hyytia_mmm = 20/(303.46-212.10)*(np.array((215.43,215.98,216.68,222.08,227.91,232.48,234.70,229.15,224.30,219.45,217.09,215.43))-212.10)
hyytia_obs = 20/(303.46-212.10)*(np.array((214.74,215.15,214.46,215.43,221.39,224.03,225.41,221.95,219.87,217.65,216.68,214.88))-212.10)
havard_mmm = 20/(204.06-112.57)*(np.array((116.59,119.64,123.38,132.11,142.10,147.78,146.81,146.38,137.66,128.65,121.30,117.97))-112.57)
havard_obs = 20/(204.06-112.57)*(np.array((120.33,118.81,121.86,126.85,128.51,137.24,142.10,136.13,129.62,124.63,120.47,119.08))-112.57)
califo_mmm = 20/(204.06-112.57)*(np.array((124.77,127.95,132.95,137.38,140.15,138.63,139.32,137.80,135.30,130.73,126.85,124.21))-112.57)
califo_obs = 20/(204.06-112.57)*(np.array((121.44,122.41,126.98,134.89,134.89,135.03,135.86,133.50,136.27,133.08,125.04,122.27))-112.57)
blodge_mmm = 20/(99.95-8.73)*(np.array((19.27,21.49,26.62,30.36,32.30,31.75,31.61,30.50,27.59,25.09,22.04,19.13))-8.73)
blodge_obs = 20/(99.95-8.73)*(np.array((19.27,21.49,27.86,35.49,50.60,58.22,64.32,59.19,36.88,32.30,29.25,23.84))-8.73)

hardacre_mmm = (ulborg_mmm, castel_mmm, auchen_mmm, hyytia_mmm, havard_mmm, califo_mmm, blodge_mmm)
hardacre_obs = (ulborg_obs, castel_obs, auchen_obs, hyytia_obs, havard_obs, califo_obs, blodge_obs)
hardacre_labels = ('Ulborg','Castel Porziano','Auchencorth Moss','Hyytiala','Harvard Forest','Citrus Orchard','Blodgett Forest')
station_lat = (56,41,55,62,42,36,39)
station_lon = (8,12,360-3,24,360-77,360-120,360-120)

# Data source
b_original_resolution = True
if b_original_resolution:
    scav_dir = "scavenging_monthly/"
    mm_dir = "monthly_means/"
else:
    scav_dir = "scavenging_monthly/regrid_hardacre/"
    mm_dir = "monthly_means/regrid_hardacre/"
experiment = ('C3RUN_default/',
              'C3RUN_mOSaic/',
              'C3RUN_mOSaic_offLight/',
              'C3RUN_mOSaic_offPhen/',
              'C3RUN_mOSaic_SWVL1/',
              'C3RUN_mOSaic_ice/',
              'C3RUN_mOSaic_desert/',
              'C3RUN_mOSaic_emis2014/',
              'C3RUN_mOSaic_hough/'
            )
data_dir = os.environ['DATA']+'/astra_data/ctm_results/' 
    
labels = ('Wesely_type',
          'mOSaic',
          'mOSaic_offLight','mOSaic_offPhen','mOSaic_SWVL1',
          'mOSaic_ice',
          'mOSaic_desert',
          'mOSaic_emis2014',
          'mOSaic_hough'
         )
colors = np.concatenate((('black',),
                         ('blue',),
                         np.repeat('cornflowerblue',3),
                         np.repeat('purple',2),
                         ('grey',),
                         ('orange',)) )

# Data source
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
        
    # Extract ozone drydeposition and rescale units
    ozone_data = [date['dry_O3']/mega for date in raw_data]
    for date in ozone_data:
        date.attrs['unit'] = 'Gg'
    # Monthly zonal total and rescale units
    try:
        ozone_zonal = [date.sum(dim='lon')/kilo for date in ozone_data] #sum(dim='time').
    except ValueError:
        ozone_zonal = [date.sum(dim='x')/kilo for date in ozone_data]
    for date in ozone_zonal:
        date.attrs['unit'] = 'Tg'
    # Extracting some general information
    # WARNING: cdo has summed these also up -> devide by number of days
    gridarea = raw_data[0]['gridarea'].isel(time=0)/31. 
    tot_surface_area = gridarea.sum()
    gridarea_faction = gridarea/tot_surface_area
    molarweight = get_molarweight(raw_data[0].isel(time=0))/31.
    
# Clean up
plt.close('all')
# Unit convertion
grampy = np.array([gridarea.data*seconds_in_month(imonth,2005)*molarweight.sel(name='O3').data for imonth in range(1,13)])*nano/giga
# Mean and standard deviation from sensitivity tests
mean_dd = ozone_data[1:]
cases = xr.concat(ozone_data[1:], dim='test')/grampy
cases_mean = cases.mean(dim='test')
cases_std = cases.std(dim='test')

# Plotting
'''
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("final-annual_total_ozone_drydep_zonal_%s" % (b_original_resolution))
ax11 = plt.subplot(131)
ax12 = plt.subplot(132)
ax13 = plt.subplot(133)
(ozone_data[0]/grampy).sel(lat=slice(56,80)).sel(lon=24, method='nearest').transpose().plot.contourf(ax=ax11, levels=np.arange(0,10.1,0.1),extend='max', cbar_kwargs={'label':'%s (%s)' %
                                                                                                                                                                      ('O$_3^{DD}$', 'nmol m$^{-2}$ s$^{-1}$'),'fraction':0.046, 'pad':0.04,'aspect':30,'orientation':'horizontal'})
(cases_mean).sel(lat=slice(56,80)).sel(lon=24, method='nearest').transpose().plot.contourf(ax=ax12, levels=np.arange(0,6.1,0.1),extend='max', cbar_kwargs={'label':'','fraction':0.046, 'pad':0.04,'aspect':30,'orientation':'horizontal'})
((cases_mean-ozone_data[0]/grampy)).sel(lat=slice(56,80)).sel(lon=24, method='nearest').transpose().plot.contourf(ax=ax13, levels=np.arange(-6.0,6.1,0.1),extend='max', cmap=plt.cm.seismic, cbar_kwargs={'label':'%s (%s)' %
                                                                                                                                                                     ('$\Delta O_3^{DD}$', 'nmol m$^{-2}$ s$^{-1}$'),'fraction':0.046, 'pad':0.04,'aspect':30,'orientation':'horizontal'})
for ax in fig1.axes[:3]:
    ax.set_xticklabels([get_month_name(imonth, length=1) for imonth in range(1,13)])
    ax.set_xlabel('')
    ax.set_ylabel('')
    
ax11.set_ylabel('Latitude (deg)')
ax12.set_yticklabels('')
ax13.set_yticklabels('')
'''
fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("final-monthly_ozone_drydep_stations_%s" % (b_original_resolution))

for i in range(1,len(station_lat)):
    ax = plt.subplot(2,3,i)
    (ozone_data[0]/grampy).sel(lat=station_lat[i],lon=station_lon[i],method='nearest').plot(ax=ax, label=labels[0])
    (cases_mean).sel(lat=station_lat[i],lon=station_lon[i],method='nearest').plot(ax=ax, label="Mean sensetivity tests", color='orange')
    ax.fill_between(cases_mean.time.data, ((cases_mean+cases_std)).sel(lat=station_lat[i],lon=station_lon[i],method='nearest'), ((cases_mean-cases_std)).sel(lat=station_lat[i],lon=station_lon[i],method='nearest'), color='orange', alpha=0.45)
    ax.scatter(ozone_data[0].time.data, hardacre_mmm[i], color='red', label='mmm (Hardacre et al., 2015)')
    ax.fill_between(ozone_data[0].time.data, hardacre_mmm[i]*4/3., hardacre_mmm[i]*2/3., color='red', alpha=0.1)
    ax.scatter(ozone_data[0].time.data, hardacre_obs[i], color='black', marker='v', label='obs (Hardacre et al., 2015)')
    ax.set_title(hardacre_labels[i])

for ax in fig2.axes:
    ticks = ax.get_xticks()
    start = ticks[0]-(ticks[1]-ticks[0])/2.
    end = ticks[-1]
    ax.set_xticks(np.linspace(start, end, 12))
    ax.set_yticks(np.arange(0,21,4))
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_ylim(0,20)
    ax.set_xticklabels('')

fig2.axes[3].set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
fig2.axes[4].set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
fig2.axes[5].set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
fig2.axes[0].set_ylabel('')
fig2.axes[3].set_ylabel('%s (%s)' %
                  ('O$_3^{DD}$', 'nmol m$^{-2}$ s$^{-1}$'), y=1)
fig2.axes[4].legend(bbox_to_anchor=(0.35, -0.25), loc=8, borderaxespad=0.,ncol=4)

'''
fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("stations_%s" % (b_original_resolution))
ax31 = plt.subplot(projection=cp.crs.PlateCarree())
delta_ozone = cases_mean-ozone_data[0]/grampy
delta_ozone.attrs['unit'] = 'Gg'

plot_ozone_drydep(ax31, delta_ozone, transform=cp.crs.PlateCarree(), mode='sum', glob=True, cmap=plt.cm.seismic,
                  cbar_kwargs={'label':'%s (%s/a)' %
                               ('O$_3^{DD}$', ozone_data[0].attrs['unit']),'orientation':'horizontal','fraction':0.046, 'pad':0.04,'aspect':30})
ax31.set_title("WESELY-<mOSaic>")
for lat,lon in zip(station_lat, station_lon):
    plt.plot((lon),(lat), marker='*', markersize='15', transform=cp.crs.Geodetic(), color='purple')
'''
fig4 = plt.figure(4, figsize=(16,9))
fig4.canvas.set_window_title("final-monthly_drydep_divergence_%s" % (b_original_resolution))

for i in range(1,len(station_lat)):
    ax = plt.subplot(2,3,i)
    ((ozone_data[0]/grampy).sel(lat=station_lat[i],lon=station_lon[i], method='nearest')-hardacre_obs[i]).plot(ax=ax, zorder=2, ls='--',
                                                                                                          label='$\mathrm{OsloCTM3_{Wesely}}-\mathrm{obs}$', color='blue')
    ((cases_mean).sel(lat=station_lat[i],lon=station_lon[i], method='nearest')-hardacre_obs[i]).plot(ax=ax, zorder=2, ls='--',
                                                                             label='$\mathrm{OsloCTM3_{<mOSaic>}}-\mathrm{obs}$', color='orange')
    #((ozone_data[0]/grampy).sel(lat=station_lat[i],lon=station_lon[i], method='nearest')-hardacre_mm[i]).plot(ax=ax, zorder=1, alpha=0.5,
    #                                                                         label='$\mathrm{OsloCTM3_{Wesely}}$-$\mathrm{Hardacre_{mmm}}$', color='blue')
    #((cases_mean).sel(lat=station_lat[i],lon=station_lon[i], method='nearest')-hardacre_mmm[i]).plot(ax=ax, zorder=1, alpha=0.5,
    #                                                                         label='$\mathrm{OsloCTM3_{<mOSaic>}}$-$\mathrm{Hardacre_{mmm}}$', color='orange')
    ax.set_title(hardacre_labels[i])
    chi2_wesely = ((((ozone_data[0]/grampy).sel(lat=station_lat[i],lon=station_lon[i], method='nearest')-hardacre_obs[i])**2).sum(dim='time')/12).data
    #chi2_mOSaic = (((((cases_mean).sel(lat=station_lat[i],lon=station_lon[i], method='nearest')-hardacre_obs[i])**2/(((cases_std).sel(lat=station_lat[i],lon=station_lon[i], method='nearest')-1)**2)).sum(dim='time'))/12).data
    chi2_mOSaic = ((((cases_mean).sel(lat=station_lat[i],lon=station_lon[i], method='nearest')-hardacre_obs[i])**2).sum(dim='time')/12).data
    ax.text(1,0.1,'$\chi^2/NDF=$%1.2f' % (chi2_wesely),
          verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, color='blue')
    ax.text(1,0,'$\chi^2/NDF=$%1.2f' % (chi2_mOSaic),
          verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, color='orange')

for ax in fig4.axes:
    ticks = ax.get_xticks()
    start = ticks[0]-(ticks[1]-ticks[0])/2.
    end = ticks[-1]
    ax.set_xticks(np.linspace(start, end, 12))
    ax.set_yticks(np.arange(-12,13,4))
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_ylim(-12,12)
    ax.set_xticklabels('')

fig4.axes[3].set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
fig4.axes[4].set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
fig4.axes[5].set_xticklabels([get_month_name(i, length=1) for i in range(1,13)])
fig4.axes[0].set_ylabel('')
fig4.axes[3].set_ylabel('%s (%s)' %
                        ('$\Delta O_3^{DD}$', 'nmol m$^{-2}$ s$^{-1}$'), y=1)
fig4.axes[4].legend(bbox_to_anchor=(0.35, -0.25), loc=8, borderaxespad=0.,ncol=4)

# Show it
plt.show(block=False)
