import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *

month = '12'
# Data source
nc_src_corr = os.environ['DATA']+'/abel/C3RUN_nitrate_corr_2005/scavenging_daily/scavenging_daily_2d_2005'+month+'*.nc'
nc_src_emep =  os.environ['DATA']+'/abel/C3RUN_emep.090318.9433/scavenging_daily/scavenging_daily_2d_2005'+month+'*.nc'

try:
    data
except NameError:
    data_list = []
    for subdir in (nc_src_corr,nc_src_emep):
        raw_data = []
        print("Reading from path %s" % (os.path.abspath(subdir)))
        # Open dataset
        for file in sorted(glob.glob(subdir)):
            print("Reading %s" % (os.path.basename(file)))
            data = xr.open_dataset(file)
            # Defining new time coordinates
            data.coords['time'] = dt.datetime(data['YEAR'], data['MONTH'], data['DAY'])
            raw_data.append(data['dry_O3'])
        # Concatenating the list
        data_list.append(xr.concat(raw_data, dim='time'))

ozone_data = [data/mega for data in data_list]
for data in ozone_data :
    data.attrs['unit'] = 'Gg'
# Plot it
plt.close('all')
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("compare-zonal_average_ozone_drydep-timeline-%s" % month)
ax11 = plt.subplot(131)
ax12 = plt.subplot(132)
ax13 = plt.subplot(133)
ozone_data[0].mean(dim='lon').plot(ax=ax11, add_colorbar=False)
ozone_data[1].mean(dim='lon').plot(ax=ax12,
                                  cbar_kwargs={'label':'%s Dry Dep (%s)' %
                                               ('O$_3$', ozone_data[1].attrs['unit'])})

(100*(ozone_data[1].mean(dim='lon')-ozone_data[0].mean(dim='lon'))/ozone_data[0].mean(dim='lon')).plot(ax=ax13,
                                            cbar_kwargs={'label':'%s Dry Dep (%s)' %
                                                         ('$\Delta O_3$', '%')})
for ax in fig1.axes:
    ax.set_xlabel('Latitude (deg)')
ax12.set_ylabel('')
ax12.set_yticklabels('')
ax13.set_ylabel('')
ax13.set_yticklabels('')

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("compare-zonal_average_ozone_drydep-%s" % month)
ax21 = plt.subplot(121)
ax22 = plt.subplot(122)
ozone_data[0].sum(dim='time').mean(dim='lon').plot(ax=ax21, label='base')
ozone_data[1].sum(dim='time').mean(dim='lon').plot(ax=ax21, label='emep')
#for itime in np.arange(len(ozone_data[1].time)):
#    ozone_data[0].isel(time=itime).mean(dim='lon').plot(ax=ax22, ls=':',color='grey')
#    ozone_data[1].isel(time=itime).mean(dim='lon').plot(ax=ax22, ls='--',color='grey')
ozone_data[0].mean(dim='time').mean(dim='lon').plot(ax=ax22, label='base')
ax22.fill_between(ozone_data[0].lat,
                  ozone_data[0].mean(dim='time').mean(dim='lon')-
                  ozone_data[0].std(dim='time').mean(dim='lon')/np.sqrt(len(ozone_data[1].time)),
                  ozone_data[0].mean(dim='time').mean(dim='lon')+
                  ozone_data[0].std(dim='time').mean(dim='lon')/np.sqrt(len(ozone_data[1].time)),
                  alpha=0.5)
ozone_data[1].mean(dim='time').mean(dim='lon').plot(ax=ax22, label='emep')
ax22.fill_between(ozone_data[0].lat,
                  ozone_data[1].mean(dim='time').mean(dim='lon')-
                  ozone_data[1].std(dim='time').mean(dim='lon')/np.sqrt(len(ozone_data[1].time)),
                  ozone_data[1].mean(dim='time').mean(dim='lon')+
                  ozone_data[1].std(dim='time').mean(dim='lon')/np.sqrt(len(ozone_data[1].time)),
                  alpha=0.5)
ax21.set_title("Zonal average of monthly sums")
ax21.legend()
ax22.set_title("Zonal average and daily variation")
ax22.legend()
for ax in fig2.axes:
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('$O_3^{drydep}$ (%s)' % ozone_data[0].unit)
# Show it
plt.show(block=False)
