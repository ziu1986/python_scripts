import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *
from collections import OrderedDict


# Closing plots from previous runs
plt.close('all')

#experiment = ('C3RUN_oDD/',
#              'C3RUN_emep_full/',
#              'C3RUN_emep_offLight/',
#              'C3RUN_emep_offPhen/',
#              'C3RUN_emep_SWVL4/',
#              'C3RUN_emep_ppgs/',
#              'C3RUN_emep_ppgssh/',
#              'C3RUN_emep_ppgssh_ice/',
#              'C3RUN_emep_ppgs_2005/')

#labels = ('OsloCTM3: Wesely type',
#          'OsloCTM3: EMEP_full',
#          'OsloCTM3: EMEP_offLight','OsloCTM3: EMEP_offPhen','OsloCTM3: EMEP_SWVL4',
#          'OsloCTM3: EMEP_ppgs','OsloCTM3: EMEP_ppgssh',
#          'OsloCTM3: EMEP_ppgssh_ice',
#          'OsloCTM3: EMEP_ppgs_2005')

#experiment = ('C3RUN_emep_ppgssh_desert_*/',)
#data_dir = os.environ['DATA']+'/abel/'
#labels = ('OsloCTM3: EMEP_ppgssh_desert_2005',)

data_dir = os.environ['DATA']+'/astra_data/ctm_results/'
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
labels = ('Wesely_type',
          'mOSaic',
          'mOSaic_offLight','mOSaic_offPhen','mOSaic_SWVL1',
          'mOSaic_ice',
          'mOSaic_desert',
          'mOSaic_emis2014',
          'mOSaic_hough'
)

macc_ozone_dir = os.environ['DATA']+'/astra_data/ECMWF/MACC_reanalysis/netcdf/monthly_mean/regrid_macc/mm_vmr_macc_r_o3_ml60_2005*.nc'

# Read the data
try:
    data
except NameError:
    ozone_data = []
    for iexp in experiment:
        subdir = data_dir+iexp+'VMR/avg_ozone_*.nc'
        data = []
        # Open dataset
        for file in sorted(glob.glob(subdir)):
            print("Reading %s" % (os.path.basename(file)))
            data.append(xr.open_dataset(file, decode_times=False))
        ozone_data.append(xr.concat(data, dim='time'))

    # Read ECWMF MACC reanalysis data
    macc_data = []
    for file in sorted(glob.glob(macc_ozone_dir)):
        data = xr.open_dataset(file)
        if (data.lat.ndim > 1):
            print("Changing coordinates...")
            data.coords['x'] = data.lat.lon[0].data
            data.coords['y'] = data.lat[:,0].data
            data = data.drop(('lon','lat'))
            data = data.rename({'x':'lon','y':'lat'})
        macc_data.append(data)
    macc_data = xr.concat(macc_data, dim='time')

# Add cyclic point
ozone_data_cyclic_da = [addcyclicpoint(data['O3'].isel(lev=0)*1e9, data['lon']) for data in (ozone_data)]
for data in ozone_data_cyclic_da:
    data.attrs['units'] = ozone_data[0]['O3'].attrs['units']

macc_ozone_data_cyclic_da = addcyclicpoint(macc_data['go3']*1e9, macc_data['lon'])
macc_ozone_data_cyclic_da.attrs['units'] = macc_data['go3'].attrs['units']
    
# Plot it
levels = np.arange(12,61,4)
levels = np.insert(levels,0,0)
levels = np.append(levels,140)
diff_levels = np.arange(-0.2,0.21,0.02)

for i in np.arange(len(experiment)): 
    fig = plt.figure(figsize=(9,16))
    fig.canvas.set_window_title("annual_mean_surface_ozone_2005_%s" % (labels[i]))
    ax11 = plt.subplot(311,projection=cp.crs.PlateCarree())
    ax12 = plt.subplot(312,projection=cp.crs.PlateCarree())
    ax13 = plt.subplot(313,projection=cp.crs.PlateCarree())

    cf11 = macc_ozone_data_cyclic_da.mean(dim='time').plot.contourf(ax=ax11, levels=levels, cmap=plt.cm.Spectral_r, transform=cp.crs.PlateCarree(), cbar_kwargs={'ticks': levels[1:-1], 'extend':'neither'})
    #cf11.cmap.set_under('yellow')
    #cf11.cmap.set_over('cyan')
    cf11.changed()
    ax11.set_title("ECWMF MACC-reanalyis")

    cf12 = ozone_data_cyclic_da[i].mean(dim='time').plot.contourf(ax=ax12, levels=levels, cmap=plt.cm.Spectral_r, transform=cp.crs.PlateCarree(), cbar_kwargs={'ticks': levels[1:-1]})
    #cf12.cmap.set_under('yellow')
    #cf12.cmap.set_over('cyan')
    cf12.changed()
    ax12.set_title("%s" % (labels[i]))

    diff_data = xr.DataArray((ozone_data_cyclic_da[i].values-macc_ozone_data_cyclic_da.values)/macc_ozone_data_cyclic_da.values, coords=[ozone_data_cyclic_da[i].time, ozone_data_cyclic_da[i].lat, ozone_data_cyclic_da[i].lon], dims=['time', 'lat', 'lon'])
    diff_data.mean(dim='time').plot.contourf(ax=ax13, levels=diff_levels, cmap=plt.cm.RdYlBu_r, transform=cp.crs.PlateCarree(), cbar_kwargs={'ticks': diff_levels[::2], 'extend':'neither'})
    ax13.set_title("%s - %s" % (labels[i],"ECWMF MACC-reanalyis"))

    for ax in fig.axes[:3]:
        ax.set_global()
        ax.set_aspect('auto')
        ax.coastlines()
        ax.set_xticks(np.arange(-180, 181, 60), crs=cp.crs.PlateCarree())
        ax.set_yticks(np.arange(-90, 91, 30), crs=cp.crs.PlateCarree())
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.set_xlabel("")
        ax.set_ylabel("")
    for ax in fig.axes[3:]:
        ax.set_ylabel("$[O_3]$ (ppb)")

    fig.axes[-1].set_ylabel("$\Delta_{exp-MACC} [O_3] /[O_3]_{MACC}$")

# Show it
plt.show(block=False)

