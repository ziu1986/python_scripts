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
          #'OsloCTM3: EMEP/MEGAN_corr','OsloCTM3: EMEP','OsloCTM3: EMEP_swgd',
          'OsloCTM3: EMEP_full',
          'OsloCTM3: EMEP_offLight','OsloCTM3: EMEP_offPhen','OsloCTM3: EMEP_SWVL4',
          'OsloCTM3: EMEP_ppgs','OsloCTM3: EMEP_ppgssh',
          'OsloCTM3: EMEP_ppgssh_ice',
          'OsloCTM3: EMEP_ppgs_2005')

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

# Add cyclic point
ozone_data_cyclic_da = [addcyclicpoint(data['O3'].isel(lev=0)*1e9, data['lon']) for data in (ozone_data)]
for data in ozone_data_cyclic_da:
    data.attrs['units'] = ozone_data[0]['O3'].attrs['units']
    
# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("annual_mean_surface_ozone_2005")
ax11 = plt.subplot(projection=cp.crs.PlateCarree())
levels = np.arange(12,61,4)
ozone_data_cyclic_da[-1].mean(dim='time').plot.contourf(ax=ax11, levels=levels, cmap=plt.cm.Spectral_r, transform=cp.crs.PlateCarree(), cbar_kwargs={'ticks': levels})
for ax in fig1.axes[::2]:
    ax.set_title("%s" % (labels[-1]))
    ax.set_global()
    ax.set_aspect('auto')
    ax.coastlines()
    ax.set_xticks(np.arange(-180, 181, 30), crs=cp.crs.PlateCarree())
    ax.set_yticks(np.arange(-90, 91, 30), crs=cp.crs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_xlabel("")
    ax.set_ylabel("")
for ax in fig1.axes[1::2]:
    ax.set_ylabel("$[O_3]$ (ppb)")

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("annual_mean_surface_ozone_2005-diff_emis")
ax21 = plt.subplot(projection=cp.crs.PlateCarree())
levels = np.arange(-10,11,1)
((ozone_data_cyclic_da[5]-ozone_data_cyclic_da[-1])/ozone_data_cyclic_da[5]*100).mean(dim='time').plot.contourf(ax=ax21, levels=levels, cmap=plt.cm.seismic, transform=cp.crs.PlateCarree(), cbar_kwargs={'ticks': levels[::2]})
for ax in fig2.axes[::2]:
    ax.set_title("(%s $-$ %s)/ %s" % (labels[5][10:], labels[-1][10:], labels[5][10:]))
    ax.set_global()
    ax.set_aspect('auto')
    ax.coastlines()
    ax.set_xticks(np.arange(-180, 181, 30), crs=cp.crs.PlateCarree())
    ax.set_yticks(np.arange(-90, 91, 30), crs=cp.crs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_xlabel("")
    ax.set_ylabel("")
for ax in fig2.axes[1::2]:
    ax.set_ylabel("$\Delta [O_3]_{2014-2005}/[O_3]_{2014}$ (%)")
# Show it
plt.show(block=False)

