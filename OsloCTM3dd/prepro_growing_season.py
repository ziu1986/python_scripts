import os, glob, sys
import numpy as np
import xarray as xr
import datetime as dt
#from growing_season import growing_season_stadyn
execfile('growing_season.py')

# Read the data
def open_data(src_path):
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(src_path)):
        print(file)
        data_raw = xr.open_dataset(file)
        data_raw.coords['time'] = dt.datetime.strptime(str(data_raw['data_date'].data),"%Y%m%d%H")
        data_list.append(data_raw['T2M'])
              
    # Concatinating the list
    print("+++>Concatinating data_list<+++")
    data = xr.concat(data_list, dim='time')
    daily_ave = []
    for itime in range(len(data.time)/8):
        # Construct the daily averages
        daily_ave.append(data.isel(time=slice(itime*8,(itime+1)*8)).mean(dim='time'))
    daily_ave = xr.concat(daily_ave, dim='time')
    daily_ave.coords['time'] = data.time[4::8]

    return(data, daily_ave)
    
# Data source
nc_src = os.environ['DATA']+'/CTM3_input_data/metdataOpenIFS/cy38r1nc4_2004/T159N80L60/0[1,2,3,4,]/*'
nc_src_2 = os.environ['DATA']+'/CTM3_input_data/metdataOpenIFS/cy38r1nc4_2004/T159N80L60/0[5,6,7,8,]/*'
nc_src_3 = os.environ['DATA']+'/CTM3_input_data/metdataOpenIFS/cy38r1nc4_2004/T159N80L60/09/*'
nc_src_4 = os.environ['DATA']+'/CTM3_input_data/metdataOpenIFS/cy38r1nc4_2004/T159N80L60/1[0,1,2,]/*'

try:
    data
except NameError:
    data_1, data_ave_1 = open_data(nc_src)
    data_2, data_ave_2 = open_data(nc_src_2)
    data_3, data_ave_3 = open_data(nc_src_3)
    data_4, data_ave_4 = open_data(nc_src_4)
   
    data = xr.concat((data_1, data_2, data_3, data_4),dim='time')
    daily_ave = xr.concat((data_ave_1, data_ave_2, data_ave_3, data_ave_4),dim='time')


# Growing season
#latitude = 61
#longitude = (162)*320/360-1
#test_location = (daily_ave-273.15).sel(lat=latitude, method='nearest').sel(lon=longitude)
#test_location_2 = (daily_ave-273.15).sel(lat=-latitude, method='nearest').sel(lon=longitude)

if(True):
    gday_ilon = []
    glen_ilon = []
    for ilon in daily_ave.lon:
        gday_ilat = []
        glen_ilat = []
        for ilat in daily_ave.lat:
            test_location = (daily_ave-273.15).sel(lat=ilat).sel(lon=ilon)
            gs = growing_season_stadyn(test_location)
            gday_ilat.append(gs[0])
            glen_ilat.append(gs[1])
        
        test_glen = xr.DataArray(glen_ilat,[('lat', daily_ave.lat.data)])
        test_gday = xr.concat(gday_ilat,dim='lat')
        test_gday.coords['lat'] = daily_ave.lat.data
        gday_ilon.append(test_gday)
        glen_ilon.append(test_glen)

gday = xr.concat(gday_ilon,dim='lon')
gday.coords['lon'] = np.arange(0,360,daily_ave.lon.data[1]-daily_ave.lon.data[0])
gday.attrs['description'] = "Days of growing."
gday.attrs['unit'] = "days"
glen = xr.concat(glen_ilon,dim='lon')
glen.coords['lon'] = np.arange(0,360,daily_ave.lon.data[1]-daily_ave.lon.data[0])
glen.attrs['description'] = "Length of the growing season."
glen.attrs['unit'] = "days"

output = xr.Dataset({'GDAY':gday,'GLEN':glen})
output.attrs['history'] = "Generated global input fields of GDAY and GLEN for OsloCTM3 from metdataOpenIFS/cy38r1nc4_2004/T159N80L60. Use fixed growing season definition for all region south of 45deg and north of 85deg. Inbetween of this use 5 days in a row above/below 5degC estimate. - S.Falk August 2018"
output.to_netcdf("GROWING_SEASON_2004_she.nc")
