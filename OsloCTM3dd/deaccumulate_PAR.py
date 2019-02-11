import os, glob, sys
import numpy as np
import xarray as xr
import datetime as dt
from scipy.constants import *     # Get physics constants
#from mytools.met_tools import *

# Create a new dataset
def deaccumulate_par(data, start12, start21):
    '''
    Deaccumulate PAR and return PPFD.
    '''
    test = [(data.isel(time=0)-(start21-start12)).copy(deep=True).data/4.*3]
    
    for i in np.arange(1,len(data.time)):
        if np.mod(i,8) > 0:
            test.append((data.isel(time=i)-data.isel(time=i-1)).copy(deep=True).data)
        else:
            test.append((data.isel(time=i)-(data.isel(time=i-1)-data.isel(time=i-4))).copy(deep=True).data)
    test_lon = [i*(data.lon[1]-data.lon[0]) for i in range(len(data.lon))]
    test_lat = data.lat.copy(deep=True).data
    test_time = data.time.isel(time=slice(0,len(data.time))).copy(deep=True).data
    test = np.array(test)
    test[np.where(test < 0)] = 0
    test_data = xr.Dataset({'PPFD':(['time','lat','lon'],test)},
                           coords={'lon':test_lon,
                                   'lat':test_lat,
                                   'time': test_time})
    test_data.attrs['units'] = 'W/m2/s, deaccumulated'
    test_data.attrs['long_name'] = 'Photosynthetic Photon Flux Density'
    test_data.attrs['history'] = 'Deaccumulated PAR using deaccumulate_PAR.py locally - S.Falk March 2018'
    
    return(test_data)

hour = 60*60
year = 2004
leapyear = (1904, 1908, 1912, 1916, 1920, 1924, 1928, 1932, 1936, 1940, 1944,
            1948, 1952, 1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988,
            1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2024, 2028, 2032,
            2036, 2040, 2044, 2048, 2052, 2056, 2060, 2064, 2068, 2072, 2076,
            2080, 2084, 2088, 2092, 2096, 2104)
months = np.arange(1,13)
imonth = 2
days = np.array((31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
if year in leapyear:
    days[1] = 29
    
# Data source
# Month that shall be processed
nc_src = os.environ['DATA']+'/CTM3_input_data/metdataOpenIFS/cy38r1nc4_'
nc_src_s1 = nc_src + str(year) + '/T159N80L60/'+ str(months[imonth]).zfill(2) +'/*'
# 21 and 12 hour data of the last day of previous month
# needed to deaccumulate the first time step at day one of the current month
nc_src_s21 = nc_src + str(year) + '/T159N80L60/'+ str(months[imonth-1]).zfill(2) +'/*d' + str(days[imonth-1]) + 'h21*'
nc_src_s12 = nc_src + str(year) + '/T159N80L60/'+ str(months[imonth-1]).zfill(2) +'/*d' + str(days[imonth-1]) + 'h12*'

try:
    data
except NameError:
    print(nc_src_s1, nc_src_s21, nc_src_s12)
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src_s1)):
        #print(file)
        data = xr.open_dataset(file)
        data.coords['time'] = dt.datetime.strptime(str(data['data_date'].data),"%Y%m%d%H")
        data.attrs['unit'] = 'W/m2/s'
        data_list.append(data['PAR']/(3*hour))
    # Longitude has 540. instead of 180.
    # Intend to correct it with this
# Concatinating the list
data = xr.concat(data_list, dim='time')
try:
    data_s1
except NameError:
    data_s21 = xr.open_dataset(glob.glob(nc_src_s21)[0])
    data_s21.coords['time'] = dt.datetime.strptime(str(data_s21['data_date'].data),"%Y%m%d%H")
    data_s21 = (data_s21['PAR']/(3*hour))
    data_s21.attrs['units'] = 'W/m2/s'
    data_s12 = xr.open_dataset(glob.glob(nc_src_s12)[0])
    data_s12.coords['time'] = dt.datetime.strptime(str(data_s12['data_date'].data),"%Y%m%d%H")
    data_s12 = (data_s12['PAR']/(3*hour))
    data_s12.attrs['units'] = 'W/m2/s'

data_new = deaccumulate_par(data, data_s12, data_s21)
data_new.to_netcdf("PPFD_%s-%s.nc" % (year, str(months[imonth]).zfill(2)))

