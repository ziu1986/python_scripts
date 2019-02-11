import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *
import matplotlib.pyplot as plt

# Show data
b_show = False

# Input/output data
nc_src = "tmp_aot40_2005.nc"

# Read input
data = read_data(nc_src)
# Generate phased longitude [0->180;-180->0]
lon_c = data.lon.data[:len(data.lon)/2]
lon_c = np.array([lon_c,-lon_c[::-1]]).ravel()
# Generate simple longitude based localtime
timedelta_c = lon_c * 24./360
# Get localtime datetimes
local_time = []
for itime in data.time.data:
    ltime = [itime+np.timedelta64(int(itdelta*60**2),'s') for itdelta in timedelta_c]
    local_time.append(ltime)
local_time = np.array(local_time)
# Save localtime variable to data
data['localtime'] = (('time','lon'), local_time)
# Select data
data_selection_nh = data['AOT40'].where((data.lat > 23) &
                                        (data.localtime.dt.month > 5) &
                                        (data.localtime.dt.month <= 8),
                                        drop=True)
data_selection_tr = data['AOT40'].where((data.lat >= -23) & (data.lat <= 23),
                                        drop=True)
data_selection_sh = data['AOT40'].where((data.lat < -23) &
                                        ~((data.localtime.dt.month > 2) &
                                        (data.localtime.dt.month <= 11)),
                                        drop=True)
data_selection = xr.concat((data_selection_sh,data_selection_tr,data_selection_nh),dim='lat')
# Resample data
aot40 = data_selection.resample('d','time',np.nansum)*1e-3
aot40.attrs['units'] = 'ppm'
aot40_820 = data_selection.where((data.localtime.dt.hour >= 8) &
                                 (data.localtime.dt.hour <= 20), drop=True).resample('d','time',np.nansum)*1e-3
aot40_820.attrs['units'] = 'ppm'

sum40 = data_selection.resample('Y','time',np.nansum)*1e-3
sum40.attrs['units'] = 'ppm'
sum40_820 = data_selection.where((data.localtime.dt.hour >= 8) &
                                 (data.localtime.dt.hour <= 20), drop=True).resample('Y','time',np.nansum)*1e-3
sum40_820.attrs['units'] = 'ppm'
# Write to file
dataset1 = xr.Dataset({'AOT40_8-20':aot40_820, 'AOT40_0-24':aot40})
dataset1.to_netcdf(nc_src[4:])
dataset2 =  xr.Dataset({'SUM40_8-20':sum40_820, 'SUM40:0-24':sum40})
dataset2.to_netcdf("sum"+nc_src[7:])

if b_show:
    # Show result
    plt.close('all')
    fig1 = plt.figure(1)
    ax11 = plt.subplot(221)
    sum40.plot()
    ax12 = plt.subplot(222)
    (data_selection_nh.resample('Y','time',np.nansum)*1e-3).plot()
    ax13 = plt.subplot(223)
    (data_selection_tr.resample('Y','time',np.nansum)*1e-3).plot()
    ax14 = plt.subplot(224)
    (data_selection_sh.resample('Y','time',np.nansum)*1e-3).plot()
    
    plt.show(block=False)

