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
          #'OsloCTM3: EMEP/MEGAN_corr','OsloCTM3: EMEP','OsloCTM3: EMEP_swgd',
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
    o3_burden = []
    for iexp in experiment:
        subdir = data_dir+iexp+'dobson_nmet.nc'
        # Open dataset
        print("Reading %s" % (os.path.basename(subdir)))
        data = xr.open_dataset(subdir, decode_times=False)
        o3_burden.append(data['o3col_trp'])
    data = xr.open_dataset(data_dir+iexp+"scavenging_monthly/sum_scavenging_200501_2d.nc")
    gridarea = data['gridarea'].isel(time=0)/31. 

du_conv = 2.687*1e20/Avogadro*48*gridarea

fig1 = plt.figure(3, figsize=(16,11))
fig1.canvas.set_window_title("mean_trop_ozone_burden")
for i in range(9):
    ax = plt.subplot(3,3,i+1,projection=cp.crs.PlateCarree())
    o3_burden[i].mean(dim='time').plot(ax=ax, levels=np.arange(0,70,5), transform=cp.crs.PlateCarree())
    ax.set_title("%s" % labels[i])
    print("%s (%3.0f +/- %2.0f) Tg" % (labels[i], (o3_burden[i].mean(dim='time')*du_conv).sum()*1e-12, ((o3_burden[i]*du_conv).sum(dim='lat').sum(dim='lon')*1e-12).std(dim='time')))

for ax in fig1.axes[::2]:
    ax.set_global()
    ax.set_aspect('auto')
    ax.coastlines()
    ax.set_xticks(np.arange(-180, 181, 90), crs=cp.crs.PlateCarree())
    ax.set_yticks(np.arange(-90, 91, 45), crs=cp.crs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_xlabel("")
    ax.set_ylabel("")
for ax in fig1.axes[1::2]:
    ax.set_ylabel("$O_3$ (DU)")
    
  
# Show it
plt.show(block=False)
