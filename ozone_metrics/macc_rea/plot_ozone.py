import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.station_info import station_location

def read_date(input_files):
    data_list = []
    print("Reading")
    for file_name in sorted(glob.glob(input_files)):
        print(file_name)
        data = xr.open_dataset(file_name)
        data_list.append(data)
        

    # Concat list
    data = xr.concat(data_list, dim='time')
    return(data)

def select(ag_data, lat,lon):
    ag_data_sel = {}
    for each in ag_data:
        ag_data_sel[each] = ag_data[each].sel(latitude=lat,longitude=lon, method='nearest')
    return(ag_data_sel)

try:
    data
except NameError:
    data = []
    for year in range(2003,2013):
        nc_src = os.environ['DATA']+"/astra_data/ECMWF/MACC_reanalysis/netcdf/*.nc" 
        data.append(select(read_date(nc_src),station_location['Esrange'].lat,station_location['Esrange'].lon))

