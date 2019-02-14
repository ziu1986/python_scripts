import os, glob
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
from scipy.constants import * # Get physics constants
from local_solar_time import local_solar_time

def resample_lst(data):
    '''
    Compute local solar time and resample data accordingly.
    '''
    time = np.arange(-1,23) # shifted time, since 23 is not valid
    lon_shifted = np.concatenate((data.lon.data[np.where(data.lon.data<180)],data.lon.data[np.where(data.lon.data>=180)]-360))
    lst = local_solar_time(time, lon_shifted) # get local solar-time wrt UTC and longitude
    x1, y1 = np.where((lst<11) & (lst>10)) # GOME equator crossing time 10.30 am = local time
    lst_lon = data.lon[x1] # note longitude and time
    lst_time = time[y1] # each time five-fold! Take only one later on.
    lst_time[np.where(lst_time==-1)] = 23 # "shift" time back
    data_sample = [] # resampling
    for iday in range(int(len(data.time)/24.)):
        data_test = []
        for i in range(len(lst_time[::5])):
            data_test.append(data.sel(time=data.time[iday*24+lst_time[::5][i]]).where(data.lon==lst_lon[5*i:5*(i+1)])) # select data accordingly
        data_test = xr.concat((data_test), dim='lon')
        data_test.coords['time'] = dt.datetime(data_test.coords['time.year'][0].data,data_test.coords['time.month'][0].data,iday+1,10)
        data_sample.append(data_test)
    data_sample = xr.concat((data_sample), dim='time')
    return data_sample

nc_src = os.environ['DATA']
subd = '/BrXplo/EMAC_total_BrO/'
src = 'BrO_col_2000*_BrXplo_ref.nc'

for file in sorted(glob.glob(nc_src+subd+src)):
    print("Reading " + file)
    data = xr.open_dataset(file)
    data_sample = resample_lst(data)
    data_sample = data_sample.groupby("time.month").mean(dim='time')
    data_sample = data_sample.mean(dim='lon')
    data_sample.to_netcdf(file[file.rfind("/")+1:])
    print("Wrote " + file[file.rfind("/")+1:])

