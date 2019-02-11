import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

hour = 60*60
# Clean up
plt.close('all')
# Data source
nc_src = os.environ['DATA']+'/input/metdataOpenIFS/cy38r1nc4_2005/T159N80L60/01/*'
nc_src_mean = os.environ['DATA']+'/input/Indata_CTM3/DRYDEP/par30T_climmean_1997-2010.nc'
try:
    data
except NameError:
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src)):
        print(file)
        data = xr.open_dataset(file)
        data.coords['time'] = dt.datetime.strptime(str(data['data_date'].data),"%Y%m%d%H")
        data_list.append(data['PAR']/(3*hour))
    # Longitude has 540. instead of 180.
    # Intend to correct it with this
    
try:
    data_mean
except NameError:
    data_mean = xr.open_dataset(nc_src_mean)
# Concatinating the list
data = xr.concat(data_list, dim='time')
maxi = []
maximum = []
for itime in range(len(data.time)):
    maximum.append(data.isel(time=itime).max())
    maxi.append(data.isel(time=itime)/data.isel(time=itime).max())
maxi = xr.concat(maxi, dim='time')
maximum = xr.concat(maximum, dim='time')

ucum_max = []
for itime in range(1,8):
    ucum_max.append(((data.isel(time=itime)-data.isel(time=itime-1))).max())

# Plot it
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("PAR_january")
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)

(data).sel(lat=60, method='nearest').sel(lon=18).plot(ax=ax11, label='60N18E')
(data).sel(lat=60, method='nearest').sel(lon=189).plot(ax=ax11, label='60N189E', color=ax11.lines[-1].get_color(), ls='--')
(data).sel(lat=0, method='nearest').sel(lon=18).plot(ax=ax11, label='018E')
(data).sel(lat=0, method='nearest').sel(lon=189).plot(ax=ax11, label='0189E', color=ax11.lines[-1].get_color(), ls='--')
(data).sel(lat=-60, method='nearest').sel(lon=18).plot(ax=ax11, label='60S18E')
(data).sel(lat=-60, method='nearest').sel(lon=189).plot(ax=ax11, label='60S189E', color=ax11.lines[-1].get_color(), ls='--')

#(maxi).sel(lat=60, method='nearest').sel(lon=18).plot(ax=ax11, label='60N', ls='--')
#(maxi).sel(lat=0, method='nearest').sel(lon=18).plot(ax=ax11, label='0',ls='--')
#(maxi).sel(lat=-60, method='nearest').sel(lon=18).plot(ax=ax11, label='60S',ls='--')

data.sel(lat=60, method='nearest').sel(lon=18).sel(time=slice('2005-01-01')).plot(ax=ax12, marker='x', label='60N18E')
data.sel(lat=60, method='nearest').sel(lon=189.0).sel(time=slice('2005-01-01')).plot(ax=ax12, color=ax12.lines[-1].get_color(), ls ='--', marker='+', label='60N189E')
data.sel(lat=60, method='nearest').sel(lon=351.0).sel(time=slice('2005-01-01')).plot(ax=ax12, color=ax12.lines[-1].get_color(), ls =':', marker='+', label='60N351E')
data.sel(lat=0, method='nearest').sel(lon=18).sel(time=slice('2005-01-01')).plot(ax=ax12, marker='x', label='018E')
data.sel(lat=0, method='nearest').sel(lon=189.0).sel(time=slice('2005-01-01')).plot(ax=ax12, color=ax12.lines[-1].get_color(), ls ='--', marker='+', label='0189E')
data.sel(lat=0, method='nearest').sel(lon=351.0).sel(time=slice('2005-01-01')).plot(ax=ax12, color=ax12.lines[-1].get_color(), ls =':', marker='+', label='0351E')
data.sel(lat=-60, method='nearest').sel(lon=18).sel(time=slice('2005-01-01')).plot(ax=ax12, marker='x', label='60S18E')
data.sel(lat=-60, method='nearest').sel(lon=189.0).sel(time=slice('2005-01-01')).plot(ax=ax12, color=ax12.lines[-1].get_color(), ls ='--', marker='+', label='60S189E')
data.sel(lat=-60, method='nearest').sel(lon=351.0).sel(time=slice('2005-01-01')).plot(ax=ax12, color=ax12.lines[-1].get_color(), ls =':', marker='+', label='60S351E')
plt.legend(ncol=3)
ax11.set_ylabel("")
ax12.set_ylabel("PAR ($W\,m^{-2}s^{-1}$)", y=1)
for ax in fig1.axes:
    ax.set_ylim(0,1200)

fig2 = plt.figure(2)
fig2.canvas.set_window_title("PAR_as_is")
ax21 = plt.subplot(421)
ax22 = plt.subplot(422)
ax23 = plt.subplot(423)
ax24 = plt.subplot(424)
ax25 = plt.subplot(425)
ax26 = plt.subplot(426)
ax27 = plt.subplot(427)
ax28 = plt.subplot(428)
(data.sel(time='2005-01-01T00')).plot(ax=ax21,vmin=0,vmax=1600)
(data.sel(time='2005-01-01T03')).plot(ax=ax22,vmin=0,vmax=1600)
(data.sel(time='2005-01-01T06')).plot(ax=ax23,vmin=0,vmax=1600)
(data.sel(time='2005-01-01T09')).plot(ax=ax24,vmin=0,vmax=1600)
(data.sel(time='2005-01-01T12')).plot(ax=ax25,vmin=0,vmax=1600)
(data.sel(time='2005-01-01T15')).plot(ax=ax26,vmin=0,vmax=1600)
(data.sel(time='2005-01-01T18')).plot(ax=ax27,vmin=0,vmax=1600)
(data.sel(time='2005-01-01T21')).plot(ax=ax28,vmin=0,vmax=1600)
time_counter = 0

for ax in fig2.axes[:-len(fig2.axes)/2]:
    sunpos = 200
    sunpos -= time_counter*15
    if sunpos < 0:
        sunpos += 360
    ax.axvline(sunpos,color='grey')
    ax.set_title('%sh' % (time_counter), x=sunpos/360.)
    time_counter += 3 
for ax in fig2.axes[:-len(fig2.axes)/2-2]:
    ax.set_xticklabels("")
    ax.set_xlabel("")
    ax.set_ylabel('')
ax27.set_xlabel('')
ax27.set_ylabel('Latitude (deg)', y=2.25)
ax28.set_ylabel('')
ax28.set_xlabel('Longitude (deg)', x=-0.25)


fig3 = plt.figure(3)
fig3.canvas.set_window_title("PAR_partly_deaccumulated")
ax31 = plt.subplot(421)
ax32 = plt.subplot(422)
ax33 = plt.subplot(423)
ax34 = plt.subplot(424)
ax35 = plt.subplot(425)
ax36 = plt.subplot(426)
ax37 = plt.subplot(427)
ax38 = plt.subplot(428)
((data.sel(time='2005-01-01T00'))/3).plot(ax=ax31,vmin=0,vmax=300)
((data.sel(time='2005-01-01T03')-data.sel(time='2005-01-01T00'))).plot(ax=ax32,vmin=0,vmax=300)
((data.sel(time='2005-01-01T06')-data.sel(time='2005-01-01T03'))).plot(ax=ax33,vmin=0,vmax=300)
((data.sel(time='2005-01-01T09')-data.sel(time='2005-01-01T06'))).plot(ax=ax34,vmin=0,vmax=300)
((data.sel(time='2005-01-01T12')-data.sel(time='2005-01-01T09'))).plot(ax=ax35,vmin=0,vmax=300)
((data.sel(time='2005-01-01T15')-data.sel(time='2005-01-01T12'))).plot(ax=ax36,vmin=0,vmax=300)
((data.sel(time='2005-01-01T18')-data.sel(time='2005-01-01T15'))).plot(ax=ax37,vmin=0,vmax=300)
((data.sel(time='2005-01-01T21')-data.sel(time='2005-01-01T18'))).plot(ax=ax38,vmin=0,vmax=300)

time_counter = 0

for ax in fig3.axes[:-len(fig3.axes)/2]:
    sunpos = 200
    sunpos -= time_counter*15
    if sunpos < 0:
        sunpos += 360
    ax.axvline(sunpos,color='grey')
    ax.set_title('%sh' % (time_counter), x=sunpos/360.)
    time_counter += 3 
for ax in fig3.axes[:-len(fig3.axes)/2-2]:
    ax.set_xticklabels("")
    ax.set_xlabel("")
    ax.set_ylabel('')
ax37.set_xlabel('')
ax37.set_ylabel('Latitude (deg)', y=2.25)
ax38.set_ylabel('')
ax38.set_xlabel('Longitude (deg)', x=-0.25)

fig4 = plt.figure(4)
fig4.canvas.set_window_title("PAR_deaccumulating_00h")
ax41 = plt.subplot(311)
ax42 = plt.subplot(312)
ax43 = plt.subplot(313)
estimate = (data.sel(time='2005-01-02T00')-(data.sel(time='2005-01-01T21')-data.sel(time='2005-01-01T12')))
data.sel(time='2005-01-02T00').plot(ax=ax41)
(data.sel(time='2005-01-01T21')-data.sel(time='2005-01-01T12')).plot(ax=ax42)
estimate.plot(ax=ax43,vmin=0,vmax=300)
for ax in fig4.axes[:2]:
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels('')
ax42.set_title("2005-01-01T21-2005-01-01T12")
ax43.set_xlabel('Longitude (deg)')
ax43.set_ylabel('Latitude (deg)', y=1.75)
ax43.set_title("2005-01-02T00-(2005-01-01T21-2005-01-01T12)")

# Create a new dataset
def deaccumulate_par(data):
    '''
    Deaccumulate PAR and return PPFD.
    '''
    test = []
    for i in np.arange(1,len(data.time)):
        if np.mod(i,8) > 0:
            test.append((data.isel(time=i)-data.isel(time=i-1)).copy(deep=True).data)
        else:
            test.append((data.isel(time=i)-(data.isel(time=i-1)-data.isel(time=i-4))).copy(deep=True).data)
    test_lon = [i*(data.lon[1]-data.lon[0]) for i in range(len(data.lon))]
    test_lat = data.lat.copy(deep=True).data
    test_time = data.time.isel(time=slice(1,len(data.time))).copy(deep=True).data
    test = np.array(test)
    test[np.where(test < 0)] = 0
    test_data = xr.Dataset({'PPFD':(['time','lat','lon'],test)},
                           coords={'lon':test_lon,
                                   'lat':test_lat,
                                   'time': test_time})
    test_data.attrs['unit'] = 'W m-2 s-1'
    return(test_data)
    
    
test = np.array([(data.sel(time='2005-01-01T03')-data.sel(time='2005-01-01T00')).copy(deep=True).data,
        (data.sel(time='2005-01-01T06')-data.sel(time='2005-01-01T03')).copy(deep=True).data,
        (data.sel(time='2005-01-01T09')-data.sel(time='2005-01-01T06')).copy(deep=True).data,
        (data.sel(time='2005-01-01T12')-data.sel(time='2005-01-01T09')).copy(deep=True).data,
        (data.sel(time='2005-01-01T15')-data.sel(time='2005-01-01T12')).copy(deep=True).data,
        (data.sel(time='2005-01-01T18')-data.sel(time='2005-01-01T15')).copy(deep=True).data,
        (data.sel(time='2005-01-01T21')-data.sel(time='2005-01-01T18')).copy(deep=True).data,
        (data.sel(time='2005-01-02T00')-(data.sel(time='2005-01-01T21')-data.sel(time='2005-01-01T12'))).copy(deep=True).data])
test_lon = data.lon.copy(deep=True).data
#test_lon[np.where(test_lon==540.)] = 180.
test_lat = data.lat.copy(deep=True).data
test_time = data.time.sel(time=slice('2005-01-01T03','2005-01-02T00')).copy(deep=True).data
test_data = xr.Dataset({'PPFD':(['time','lat','lon'],test)},
                       coords={'lon':test_lon,
                               'lat':test_lat,
                               'time': test_time})
test_data.attrs['unit'] = 'W m^{-2} s^{-1}'
test = deaccumulate_par(data)

fig5 = plt.figure(5,figsize=(16,9))
fig5.canvas.set_window_title("PPFD_january_deaccumulated")
ax52 = plt.subplot(211)
ax51 = plt.subplot(212)
test_data.sel(lat=60, method='nearest').sel(lon=18)['PPFD'].plot(ax=ax51, marker='x', label='60N18E')
test_data.sel(lat=60, method='nearest').sel(lon=189.0)['PPFD'].plot(ax=ax51, color=ax51.lines[-1].get_color(), ls ='--', marker='+', label='60N189E')
test_data.sel(lat=60, method='nearest').sel(lon=351.0)['PPFD'].plot(ax=ax51, color=ax51.lines[-1].get_color(), ls =':', marker='+', label='60N351E')
test_data.sel(lat=0, method='nearest').sel(lon=18)['PPFD'].plot(ax=ax51, marker='x', label='018E')
test_data.sel(lat=0, method='nearest').sel(lon=189.0)['PPFD'].plot(ax=ax51, color=ax51.lines[-1].get_color(), ls ='--', marker='+', label='0189E')
test_data.sel(lat=0, method='nearest').sel(lon=351.0)['PPFD'].plot(ax=ax51, color=ax51.lines[-1].get_color(), ls =':', marker='+', label='0351E')
test_data.sel(lat=-60, method='nearest').sel(lon=18)['PPFD'].plot(ax=ax51, marker='x', label='60S18E')
test_data.sel(lat=-60, method='nearest').sel(lon=189.0)['PPFD'].plot(ax=ax51, color=ax51.lines[-1].get_color(), ls ='--', marker='+', label='60S189E')
test_data.sel(lat=-60, method='nearest').sel(lon=351.0)['PPFD'].plot(ax=ax51, color=ax51.lines[-1].get_color(), ls =':', marker='+', label='60S351E')

test['PPFD'].sel(lat=60, method='nearest').sel(lon=18).plot(ax=ax52, label='60N18E')
test['PPFD'].sel(lat=60, method='nearest').sel(lon=189).plot(ax=ax52, label='60N189E', color=ax52.lines[-1].get_color(), ls='--')
test['PPFD'].sel(lat=0, method='nearest').sel(lon=18).plot(ax=ax52, label='018E')
test['PPFD'].sel(lat=0, method='nearest').sel(lon=189).plot(ax=ax52, label='0189E', color=ax52.lines[-1].get_color(), ls='--')
test['PPFD'].sel(lat=-60, method='nearest').sel(lon=18).plot(ax=ax52, label='60S18E')
test['PPFD'].sel(lat=-60, method='nearest').sel(lon=189).plot(ax=ax52, label='60S189E', color=ax52.lines[-1].get_color(), ls='--')

ax51.legend(ncol=3)
ax51.set_ylabel("PPFD ($W\,m^{-2}\,s^{-1}$)",y=1.)
ax52.set_ylabel('')
for ax in fig5.axes:
    ax.set_ylim(0,300)

#data.sel(time=slice('2005-01-01')).plot(col='time', col_wrap=2)
#data_mean['PAR'].plot(col='time', col_wrap=4)

# Show it
plt.show(block=False)
