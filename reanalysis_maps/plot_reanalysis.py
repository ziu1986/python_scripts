import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import xarray as xr
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from matplotlib.dates import date2num # Convert dates to matplotlib axis coords
from matplotlib import dates
import cartopy as cart
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#import cartopy.feature as cfeature
#import cartopy.io.img_tiles as cimgt

from mytools.ozone_tools import *
from mytools.plot_tools import print_all

def sample_to_da(clim_sample, lat, lon, time, **karg):
    b_camsaq = karg.pop('CAMSAQ', False)
    
    clim_3d = [np.meshgrid(lon, lat, each)[2] for each in clim_sample]
    clim_3d_squeezed = np.array(clim_3d).squeeze()
    time = time.resample(time='D').first().data
    if b_camsaq:
        da_clim_3d = xr.DataArray(
            data=clim_3d_squeezed,
            dims=["time", 'lat', 'lon'],
            coords=dict(
                time=time,
                lat=lat,
                lon=lon,
                #reference_time=reference_time,
            ),
            attrs=dict(
                description="Fennoscandic ozone climatology.",
                units="ppb",
            ),
        )
    else:
        da_clim_3d = xr.DataArray(
            data=clim_3d_squeezed,
            dims=["time", 'latitude', 'longitude'],
            coords=dict(
                time=time,
                latitude=lat,
                longitude=lon,
                #reference_time=reference_time,
            ),
            attrs=dict(
                description="Fennoscandic ozone climatology.",
                units="ppb",
            ),
        )

    return(da_clim_3d)

def weighted_seasonal_means(data, **karg):
    month_length = data.time.dt.days_in_month
    # Calculate the weights by grouping by 'time.season'.
    weights = month_length.groupby('time.season') / month_length.groupby('time.season').sum().astype(float)
    #print(weights)
    # Calculate the weighted average
    ds_weighted = (data * weights).groupby('time.season').sum(dim='time')

    return(ds_weighted)
    
def plot_map(data, **karg):
    # This is the map projection we want to plot *onto*
    map_proj = ccrs.PlateCarree(central_longitude=18.5)

    levels = karg.pop('levels', np.arange(0,61,2.5))
    clabel = karg.pop('clabel', '$[O_3]$ (ppb)')
    title = karg.pop('title', None)
    #cmap = karg.pop('cmap', )

    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()

    fig2 = data.plot(levels=levels, transform=ccrs.PlateCarree(), robust=True, col='season', col_wrap=2,figsize=(16,8), cbar_kwargs={'label': clabel, 'shrink':.7, 'anchor':(0.5,0.5), 'drawedges':True}, subplot_kws={'projection': map_proj})

    for ax in fig2.axes.flat:
        ax.coastlines(resolution='10m')
        ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', color='lightgrey')
        ax.set_extent([13.75, 33.5, 64.75, 71.75])
        ax.set_xticks(np.arange(14,34,2), crs=ccrs.PlateCarree())
        ax.set_yticks(np.arange(65, 72,  2), crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.set_xlabel("")
        ax.set_ylabel("")


        fig2.axes.flat[-2].set_xlabel("Longitude ($^\circ E$)", x=1)
        fig2.axes.flat[-2].set_ylabel("Latitude ($^\circ N$)", y=1)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        fig2.fig.canvas.set_window_title(title)
        #plt.tight_layout()

    


plt.close('all')

dataset = 'CAMSAQ'

src = {'CAMS': os.environ['DATA']+'/nird_data/reanalysis/ECMWF/CAMS_reanalysis/netcdf/VMR/vmr_cams_r_o3_ml60_climatology.nc',
       'MACC': os.environ['DATA']+'/nird_data/reanalysis/ECMWF/MACC_reanalysis/netcdf/VMR/vmr_macc_r_o3_ml60_3h_climatology.nc',
       'CAMSAQ': os.environ['DATA']+'/nird_data/reanalysis/Copernicus/ensemble_ozone/SCA_ENSa.yearlyrea.climatology.nc'}

var = {'CAMS': 'go3',
       'MACC': 'go3',
       'CAMSAQ': 'O3'}



data = xr.open_dataset(src[dataset])

if dataset == 'CAMSAQ':
    selection = data.sel(lat=slice(65,71.5), lon=slice(14,33))*0.5
else:
    selection = data.sel(latitude=slice(71.5,65), longitude=slice(14,33))*1e9
    
date_freq = pd.date_range('2012-01-01 00:00:00', '2013', freq='1H')[:-1].size/selection.time.size
new_time = pd.date_range('2012-01-01 00:00:00', '2013', freq='%sH' % date_freq)[:-1]
selection['time'] = new_time

selection_seasonal_means = weighted_seasonal_means(selection)
selection_annual_mean = selection.mean(dim='time')

import pickle
with open( os.environ['PY_SCRIPTS']+'/ozone_metrics/ozone_anomalies/obs_climatologies.pkl', 'rb') as input:
    climatology_nord = pickle.load(input)
    climatology_svanvik = pickle.load(input)
    
sample_spl_nord = climatology_nord(np.unique(selection.time.dt.dayofyear))

if dataset == 'CAMSAQ':
    da_clim_3d = sample_to_da(sample_spl_nord, selection.lat.data, selection.lon.data, selection.time, CAMSAQ=True)
    
else:
    da_clim_3d = sample_to_da(sample_spl_nord, selection.latitude, selection.longitude, selection.time)
    
diff = (selection.resample(time='D').mean())-da_clim_3d
diff_sig = diff/1

#diff['time'] = selection.time[1:].resample(time='D').first().data


diff_seasonal_mean = weighted_seasonal_means(diff)
diff_sig_seasonal_mean = weighted_seasonal_means(diff_sig)
    
plot_map(selection_seasonal_means[var[dataset]], title="ozone_seasonal_average_%s" % dataset)
plot_map(diff_seasonal_mean[var[dataset]], title="ozone_seasonal_diff_%s" % dataset, clabel='$\Delta[O_3]$ (ppb)', levels=np.arange(-20, 21, 1))
plot_map(diff_sig_seasonal_mean[var[dataset]], title="ozone_seasonal_diff_sig_%s" % dataset, clabel='$\Delta[O_3]$ $(\sigma_{clim})$', levels=np.arange(-20, 21, 1))
plt.show(block=False)
