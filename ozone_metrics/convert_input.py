import os, glob, sys
import numpy as np
import xarray as xr
from mytools.netcdf_tools import *
for iyear in np.arange(2010,2012):
    nc_src = os.environ['DATA']+"/CTM3_Oivind_data/CIXPAG_%s/ctm3_icbc_glo/gas*%s*.nc" % (iyear,iyear)
    data = read_data(nc_src,var='NO2',datatype='osloctm_ozone')
    
    if len(np.where(data.lon.data==540)[0])>0:
        test_lon = [i*(data.lon[1]-data.lon[0]) for i in range(len(data.lon))]
        test_lat = data.lat.copy(deep=True).data
        test_lev = data.lev.copy(deep=True).data
        test_time = data.time.copy(deep=True).data
        test_data = xr.Dataset({'NO2':(['time','lev','lat','lon'],data.copy(deep=True).data)},
                               coords={'lon':test_lon,
                                       'lat':test_lat,
                                       'lev':test_lev,
                                       'time': test_time})
        test_data.attrs['units'] = data.attrs['units']
        data = test_data

    # Save it to new file
    new_file = "osloctm_no2_%s.nc" % (iyear)
    print('Save to file...%s' % new_file)
    data.to_netcdf(new_file)
