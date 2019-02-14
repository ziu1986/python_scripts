import os, glob
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import * # Get physics constants
import matplotlib.pyplot as plt
import matplotlib.dates as mdates # handling of dates in plotting
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime as dt  # Python standard library datetime module
from mytools.med_tools import *
from mytools.netcdf_tools import * # ncdump implementation for python
from local_solar_time import local_solar_time
from plot_BrO_col import resample_lst


# Compute space masked simulation data
def space_masked(input_field, lat_mask, hemisphere='nh'):
    output_field = []
    if hemisphere=='nh':
        output_field.append(input_field.sel(time=slice('2000-01-01', '2000-01-31')).where(input_field.lat<=lat_mask[0], drop=True))
        output_field.append(input_field.sel(time=slice('2000-02-01', '2000-02-29')).where(input_field.lat<=lat_mask[1], drop=True))
        output_field.append(input_field.sel(time=slice('2000-03-01', '2000-03-31')).where(input_field.lat<=lat_mask[2], drop=True))
        output_field.append(input_field.sel(time=slice('2000-04-01', '2000-04-30')).where(input_field.lat<=lat_mask[3], drop=True))
        output_field.append(input_field.sel(time=slice('2000-05-01', '2000-05-31')).where(input_field.lat<=lat_mask[4], drop=True))
        output_field.append(input_field.sel(time=slice('2000-06-01', '2000-06-30')).where(input_field.lat<=lat_mask[5], drop=True))
        output_field.append(input_field.sel(time=slice('2000-07-01', '2000-07-31')).where(input_field.lat<=lat_mask[6], drop=True))
        output_field.append(input_field.sel(time=slice('2000-08-01', '2000-08-31')).where(input_field.lat<=lat_mask[7], drop=True))
        output_field.append(input_field.sel(time=slice('2000-09-01', '2000-09-30')).where(input_field.lat<=lat_mask[8], drop=True))
        output_field.append(input_field.sel(time=slice('2000-10-01', '2000-10-31')).where(input_field.lat<=lat_mask[9], drop=True))
        output_field.append(input_field.sel(time=slice('2000-11-01', '2000-11-30')).where(input_field.lat<=lat_mask[10], drop=True))
        output_field.append(input_field.sel(time=slice('2000-12-01', '2000-12-31')).where(input_field.lat<=lat_mask[11], drop=True))
    else:
        output_field.append(input_field.sel(time=slice('2000-01-01', '2000-01-31')).where(input_field.lat>=lat_mask[0], drop=True))
        output_field.append(input_field.sel(time=slice('2000-02-01', '2000-02-29')).where(input_field.lat>=lat_mask[1], drop=True))
        output_field.append(input_field.sel(time=slice('2000-03-01', '2000-03-31')).where(input_field.lat>=lat_mask[2], drop=True))
        output_field.append(input_field.sel(time=slice('2000-04-01', '2000-04-30')).where(input_field.lat>=lat_mask[3], drop=True))
        output_field.append(input_field.sel(time=slice('2000-05-01', '2000-05-31')).where(input_field.lat>=lat_mask[4], drop=True))
        output_field.append(input_field.sel(time=slice('2000-06-01', '2000-06-30')).where(input_field.lat>=lat_mask[5], drop=True))
        output_field.append(input_field.sel(time=slice('2000-07-01', '2000-07-31')).where(input_field.lat>=lat_mask[6], drop=True))
        output_field.append(input_field.sel(time=slice('2000-08-01', '2000-08-31')).where(input_field.lat>=lat_mask[7], drop=True))
        output_field.append(input_field.sel(time=slice('2000-09-01', '2000-09-30')).where(input_field.lat>=lat_mask[8], drop=True))
        output_field.append(input_field.sel(time=slice('2000-10-01', '2000-10-31')).where(input_field.lat>=lat_mask[9], drop=True))
        output_field.append(input_field.sel(time=slice('2000-11-01', '2000-11-30')).where(input_field.lat>=lat_mask[10], drop=True))
        output_field.append(input_field.sel(time=slice('2000-12-01', '2000-12-31')).where(input_field.lat>=lat_mask[11], drop=True))
    return xr.concat(output_field, dim='time')

# Clean up
plt.close('all')
nc_src = os.environ['DATA']
# EMAC BrO column density
subd_emac = '/processed_data/BrXplo/'
src_emac = 'BrO_col_2000*_BrXplo_mysic_corr.nc'
src_emac_ref = 'BrO_col_2000*_BrXplo_ref.nc'
src_emac_fysic = 'BrO_col_2000*_BrXplo.nc'
src_emac_mysic_rs = 'BrO_col_2000*_BrXplo_mysic_rsnow.nc'
# GOME BrO column density
subd_gome = '/BrO/gome_2000/'
src_gome = 'gome_*.nc'

month_name = ("January", "February", "March", "April", "May", "June", 
              "July", "August", "September", "October", "November", "December")
lat_mask_nh = (60,70,80,90,90,90,90,90,82,75,65,58)
lat_mask_sh = (-86,-86,-78,-68,-58,-53,-55,-65,-76,-86,-86,-86)
try:
    EMAC_data_list
except NameError:
    # Read EMAC and re-sample data
    EMAC_data_list = []
    EMAC_data_list_ref = []
    EMAC_data_list_fysic = []
    EMAC_data_list_mysic_rs = []
    # Open dataset
    for file in sorted(glob.glob(nc_src+subd_emac+src_emac)):
        EMAC_data_list.append(resample_lst(xr.open_dataset(file)))
    for file in sorted(glob.glob(nc_src+subd_emac+src_emac_ref)):
        EMAC_data_list_ref.append(resample_lst(xr.open_dataset(file)))
    for file in sorted(glob.glob(nc_src+subd_emac+src_emac_fysic)):
        EMAC_data_list_fysic.append(resample_lst(xr.open_dataset(file)))
    for file in sorted(glob.glob(nc_src+subd_emac+src_emac_mysic_rs)):
        EMAC_data_list_mysic_rs.append(resample_lst(xr.open_dataset(file)))
    # Read GOME data
    GOME_data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src+subd_gome+src_gome)):
        GOME_data_list.append(xr.open_dataset(file))
    # Concatenate the datasets
    EMAC_data_list = xr.concat(EMAC_data_list, dim='time')
    EMAC_data_list_ref = xr.concat(EMAC_data_list_ref, dim='time')
    EMAC_data_list_fysic = xr.concat(EMAC_data_list_fysic, dim='time')
    EMAC_data_list_mysic_rs = xr.concat(EMAC_data_list_mysic_rs, dim='time')
    GOME_data_list = xr.concat(GOME_data_list, dim='time')
    # Compute temporal means
    EMAC_time_mean = EMAC_data_list.groupby('time.month').mean(dim='time')
    EMAC_time_mean_ref = EMAC_data_list_ref.groupby('time.month').mean(dim='time')
    EMAC_time_mean_fysic = EMAC_data_list_fysic.groupby('time.month').mean(dim='time')
    EMAC_time_mean_mysic_rs = EMAC_data_list_mysic_rs.groupby('time.month').mean(dim='time')
    GOME_time_mean = GOME_data_list.groupby('time.month').mean(dim='time')

# Plot histogramms
for month in range(1,1):
    fig1 = plt.figure(month)
    fig1.canvas.set_window_title("hist_compare_BrO_column_%s" % (month_name[month-1]))
    ax11 = plt.subplot(211)
    ax11.set_title("Northern Hemisphere", x=0.65, y=0.85)
    GOME_hist_nh = ((GOME_time_mean['BrO']).sel(month=month)).where(GOME_time_mean.lat>45).plot.hist(bins=400, normed=True, range=(0,8), histtype='step', label="GOME", color='red')
    EMAC_hist_nh_ref = ((EMAC_time_mean_ref['BrO_column']*1e-13).sel(month=month)).where((EMAC_time_mean_ref.lat>45)&(EMAC_time_mean.lat<=lat_mask_nh[month-1])).plot.hist(bins=400, normed=True, range=(0,8), histtype='step', label="EMAC_ref", color='black')
    EMAC_hist_nh_fysic = ((EMAC_time_mean_fysic['BrO_column']*1e-13).sel(month=month)).where((EMAC_time_mean_fysic.lat>45)&(EMAC_time_mean.lat<=lat_mask_nh[month-1])).plot.hist(bins=400, normed=True, range=(0,8), histtype='step', label="EMAC_fysic", color='cornflowerblue')
    EMAC_hist_nh = (((EMAC_time_mean['BrO_column'])*1e-13).sel(month=month)).where((EMAC_time_mean.lat>45)&(EMAC_time_mean.lat<=lat_mask_nh[month-1])).plot.hist(bins=400, normed=True, range=(0,8), histtype='step', label="EMAC_mysic", color='blue')
    plt.legend(loc='best', frameon=False)
    ax12 = plt.subplot(212)
    ax12.set_title("Southern Hemisphere", x=0.65, y=0.85)
    GOME_hist_sh = ((GOME_time_mean['BrO']).sel(month=month)).where(GOME_time_mean.lat<-45).plot.hist(bins=400, normed=True, range=(0,8), histtype='step', label="GOME", color='red')
    EMAC_hist_sh_ref = ((EMAC_time_mean_ref['BrO_column']*1e-13).sel(month=month)).where((EMAC_time_mean.lat<-45)&(EMAC_time_mean.lat>=lat_mask_sh[month-1])).plot.hist(bins=400, normed=True, range=(0,8), histtype='step', label="EMAC_ref", color='black')
    EMAC_hist_sh_fysic = ((EMAC_time_mean_fysic['BrO_column']*1e-13).sel(month=month)).where((EMAC_time_mean.lat<-45)&(EMAC_time_mean.lat>=lat_mask_sh[month-1])).plot.hist(bins=400, normed=True, range=(0,8), histtype='step', label="EMAC_fysic", color='cornflowerblue')
    EMAC_hist_sh = (((EMAC_time_mean['BrO_column'])*1e-13).sel(month=month)).where((EMAC_time_mean.lat<-45)&(EMAC_time_mean.lat>=lat_mask_sh[month-1])).plot.hist(bins=400, normed=True, range=(0,8), histtype='step', label="EMAC_mysic", color='blue')

    for ax in fig1.axes:
        ax.set_xlabel("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'))
        ax.set_ylabel("Probability Density")
        ax.set_ylim(0,1)

fig13 = plt.figure(13)
fig13.canvas.set_window_title("BrO_col_zonal_mean_nh")
ax131 = plt.subplot()
fysic_nh = EMAC_data_list_fysic.where(EMAC_data_list_fysic.lat>45, drop=True)
fysic_nh_masked = space_masked(fysic_nh, lat_mask_nh, hemisphere='nh')
mysic_nh = EMAC_data_list.where(EMAC_data_list.lat>45, drop=True)
mysic_nh_masked = space_masked(mysic_nh, lat_mask_nh, hemisphere='nh')
ref_nh = EMAC_data_list_ref.where(EMAC_data_list_ref.lat>45, drop=True)
ref_nh_masked = space_masked(ref_nh, lat_mask_nh, hemisphere='nh')
mysic_rs_nh = EMAC_data_list_mysic_rs.where(EMAC_data_list.lat>45, drop=True)
mysic_rs_nh_masked = space_masked(mysic_rs_nh, lat_mask_nh, hemisphere='nh')
gome_nh = GOME_data_list.where(GOME_data_list.lat>45, drop=True)

average(fysic_nh.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(fysic_nh.lat*np.pi/180)).plot(zorder=3, color='cornflowerblue', ls='--', label='fysic')
average(fysic_nh_masked.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(fysic_nh_masked.lat*np.pi/180)).plot(zorder=3, color='cornflowerblue', ls=':', label='fysic_masked')
average(mysic_nh.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(mysic_nh.lat*np.pi/180)).plot(zorder=3, color='blue', label='mysic')
average(mysic_nh_masked.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(mysic_nh_masked.lat*np.pi/180)).plot(zorder=3, color='blue', ls=':', label='mysic_masked')
average(ref_nh.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(ref_nh.lat*np.pi/180)).plot(zorder=3, color='black', label='ref')
average(ref_nh_masked.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(ref_nh_masked.lat*np.pi/180)).plot(zorder=3, color='black', ls=':', label='ref_masked')
average(mysic_rs_nh.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(mysic_rs_nh.lat*np.pi/180)).plot(zorder=3, color='cyan', label='mysic_rs')
average(mysic_rs_nh_masked.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(mysic_rs_nh_masked.lat*np.pi/180)).plot(zorder=3, color='cyan', ls=':', label='mysic_rs_masked')
average(gome_nh.where(gome_nh.BrO>=0).mean(dim='lon')['BrO'], 
        dim='lat', weights=np.cos(gome_nh.lat*np.pi/180)).plot(zorder=2, ls='none', marker='x', color='red', label='obs')
ax131.set_ylabel("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'))
ax131.set_xlabel("Time")
ax131.set_ylim(0,6.5)
ax131.set_title("Northern Hemisphere", x=0.2, y=0.05)
plt.legend(frameon=False, ncol=5, loc='best')

fig14 = plt.figure(14)
fig14.canvas.set_window_title("BrO_col_zonal_mean_sh")
ax141 = plt.subplot()
fysic_sh = EMAC_data_list_fysic.where(EMAC_data_list_fysic.lat<-45, drop=True)
fysic_sh_masked = space_masked(fysic_sh, lat_mask_sh, hemisphere='sh')
mysic_sh = EMAC_data_list.where(EMAC_data_list.lat<-45, drop=True)
mysic_sh_masked = space_masked(mysic_sh, lat_mask_sh, hemisphere='sh')
ref_sh = EMAC_data_list_ref.where(EMAC_data_list_ref.lat<-45, drop=True)
ref_sh_masked = space_masked(ref_sh, lat_mask_sh, hemisphere='sh')
mysic_rs_sh = EMAC_data_list_mysic_rs.where(EMAC_data_list.lat<-45, drop=True)
mysic_rs_sh_masked = space_masked(mysic_rs_sh, lat_mask_sh, hemisphere='sh')
gome_sh = GOME_data_list.where(GOME_data_list.lat<-45, drop=True)
average(fysic_sh.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(fysic_sh.lat*np.pi/180)).plot(zorder=3, color='cornflowerblue', ls='--', label='fysic')
average(fysic_sh_masked.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(fysic_sh_masked.lat*np.pi/180)).plot(zorder=3, color='cornflowerblue', ls=':', label='fysic_masked')
average(mysic_sh.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(mysic_sh.lat*np.pi/180)).plot(zorder=3, color='blue', label='mysic')
average(mysic_sh_masked.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(mysic_sh_masked.lat*np.pi/180)).plot(zorder=3, color='blue', ls=':', label='mysic_masked')
average(ref_sh.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(ref_sh.lat*np.pi/180)).plot(zorder=3, color='black', label='ref')
average(ref_sh_masked.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(ref_sh_masked.lat*np.pi/180)).plot(zorder=3, color='black', ls=':', label='ref_masked')
average(mysic_rs_sh.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(mysic_rs_sh.lat*np.pi/180)).plot(zorder=3, color='cyan', label='mysic_rs')
average(mysic_rs_sh_masked.mean(dim='lon')['BrO_column']*1e-13, 
        dim='lat', weights=np.cos(mysic_rs_sh_masked.lat*np.pi/180)).plot(zorder=3, color='cyan', ls=':', label='mysic_rs_masked')
average(gome_sh.where(gome_sh.BrO>=0).mean(dim='lon')['BrO'], 
        dim='lat', weights=np.cos(gome_sh.lat*np.pi/180)).plot(zorder=2, ls='none', marker='x', color='red', label='obs')
ax141.set_ylabel("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'))
ax141.set_xlabel("Time")
ax141.set_ylim(0,6.5)
ax141.set_title("Southern Hemisphere", x=0.2, y=0.05)
plt.legend(frameon=False, ncol=5, loc='best')

fig15 = plt.figure(15)
fig15.canvas.set_window_title("BrO_col_zonal_means_gome")
ax151 = plt.subplot(211)
color_range = np.arange(0,5.1,0.5)
gome_nh_zonal_mean = gome_nh.where(gome_nh.BrO>=0).mean(dim='lon')
x_time = gome_nh_zonal_mean['time']
csp151 = xr.plot.contourf(gome_nh_zonal_mean['BrO'], x='time', y='lat', levels=color_range, extend='both')
ax152 = plt.subplot(212)
gome_sh_zonal_mean = gome_sh.where(gome_sh.BrO>=0).mean(dim='lon')
csp152 = xr.plot.contourf(gome_sh_zonal_mean['BrO'], x='time', y='lat', levels=color_range, extend='both')
#cbar = fig15.colorbar(csp151, ax=(ax151,ax152), orientation='vertical', fraction=0.046, pad=0.04, aspect=30)
fig15.autofmt_xdate()

fig16 = plt.figure(16)
fig16.canvas.set_window_title("BrO_col_zonal_means_diff")
ax161 = plt.subplot(211)
color_range = np.arange(0,5.1,0.5)
mysic_nh_zonal_mean = (mysic_nh-ref_nh).mean(dim='lon')*1e-13
x_time = mysic_nh_zonal_mean['time']
csp161 = xr.plot.contourf(mysic_nh_zonal_mean['BrO_column'], x='time', y='lat', levels=color_range, extend='both', cmap=plt.cm.viridis)
ax162 = plt.subplot(212)
mysic_sh_zonal_mean = (mysic_sh-ref_sh).mean(dim='lon')*1e-13
csp162 =xr.plot.contourf(mysic_sh_zonal_mean['BrO_column'], x='time', y='lat', levels=color_range, extend='both', cmap=plt.cm.viridis)
#cbar = fig16.colorbar(csp161, ax=(ax161,ax162), orientation='vertical', fraction=0.046, pad=0.04, aspect=30)
fig16.autofmt_xdate()

fig17 = plt.figure(17, figsize=(16,9))
month_name = ("January", "February", "March", "April", "May", "June", 
              "July", "August", "September", "October", "November", "December")
fig17.canvas.set_window_title("BrO_col_zonal_mean_monthly")
for imonth in range(1,13):
    ax = plt.subplot(3,4,imonth)
    mysic_select = EMAC_time_mean.mean(dim='lon').sel(month=imonth)
    mysic_rs_select = EMAC_time_mean_mysic_rs.mean(dim='lon').sel(month=imonth)
    fysic_select = EMAC_time_mean_fysic.mean(dim='lon').sel(month=imonth)
    ref_select = EMAC_time_mean_ref.mean(dim='lon').sel(month=imonth)
    gome_select = GOME_time_mean.mean(dim='lon').sel(month=imonth)
    ((mysic_select['BrO_column']-ref_select['BrO_column'])*1e-13).plot(label='mysic-ref')
    ((fysic_select['BrO_column']-ref_select['BrO_column'])*1e-13).plot(color='cornflowerblue', ls='--', label='fysic-ref')
    ((mysic_rs_select['BrO_column']-ref_select['BrO_column'])*1e-13).plot(label='mysic_rs-ref', color='cyan', ls=':')
    #(ref_select['BrO_column']*1e-13).plot(color='black', label='ref')
    gome_select['BrO'].plot(color='red', ls='none', marker='x', label='obs')
    ax.set_title(month_name[imonth-1], ha='left', y=0.85, x=0.02)

for ax in fig17.axes:
    ax.set_xticks(np.arange(-90,91,45))
    ax.set_xlabel("")
    ax.set_ylim(0,6)
    ax.set_ylabel("")
    
ax.legend(loc='best', frameon=False)
fig17.axes[4].set_ylabel("Latitude (deg)")
fig17.axes[10].set_xlabel("Longitude (deg)", x=-0.25)
# Show it
plt.show(block=False)
