import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
from matplotlib.dates import date2num
import datetime as dt
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from scipy.stats import linregress     # Get linearregression stats
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Get the read-in routine for EBAS data
def read_station_data(infile, **kargs):
    '''
    Conversion of NASA aimes files from nilo.no data sets.
    '''
    import nappy as nap # Read and write NASA data files
    # Open file (read header only)
    nas_file = nap.openNAFile(infile)
    # Read actual data
    nas_file.readData()
    # Access the ozone data
    data_raw = np.array(nas_file.getNADict()['V'])
    # Close the nas-file
    nas_file.close()
    # Keywords
    conversion = kargs.pop('conversion', True)
    tracer = kargs.pop('tracer', ('O3',))
    verbose = kargs.pop('v', False)
    if (data_raw.shape[0]-1)/2 < len(tracer):
        for each in tracer:
            if nas_file.getNADict()['NCOM'][-1].find(each) > 0:
                return(read_station_data(infile, tracer=(each,)))
    # Conversion = 1/air_dens(25degC)/mass_fraction
    # 1/2 for O3 (0.5)
    # ca. 0.38 for SO2 
    # ca. 0.81 for NO
    M_air = 28.949                        # [g/mol] 
    air_dens_25 = 1.1839                  # [kg/m3]
    mass_fraction = {'O3':3*15.9994/M_air,
                     'SO2':32.065/M_air,# ug(S)/m3
                     'SO4':32.065/M_air,
                     'NO':14.0067/M_air,# ug(N)/m3
                     'NO2':14.0067/M_air}
       
    nas_date = nas_file.getNADict()['DATE']
    start_date = dt.datetime.strptime("%s-0%s-0%s 00:00:00" % (nas_date[0], nas_date[1], nas_date[2]), '%Y-%m-%d %H:%M:%S')
    # Filter the data (FLAG: 0 - valid, >0 - invalid)
    # Add day fraction at stop time to start date 
    x_time_station = np.ma.masked_where(data_raw[2]>0, 
                                        datetime_from_time(start_date, data_raw[0])[0])
    data_dic = {'time':x_time_station}
    if verbose:
        for each in tracer:
            print(each, data_raw.shape,
                  (nas_file.getNADict()['NCOM'][-1]).split().index(each),
                  (nas_file.getNADict()['NCOM'][-1]).split()[(nas_file.getNADict()['NCOM'][-1]).split().index(each)])
    #i = 1
    for each in tracer:
        # Find column in which tracer appears
        # Split table header and look for it
        # Table header has start/end time while
        # nas_file.getNADict()['V'] only gives one time field
        i = (nas_file.getNADict()['NCOM'][-1]).split().index(each)-1
        if conversion:
            data_dic[each] = np.ma.masked_where(data_raw[i+1]>0, data_raw[i]*1/air_dens_25/mass_fraction[each])
        else:
            data_dic[each] = np.ma.masked_where(data_raw[i+1]>0, data_raw[i])
        #i += 1
    
    return(data_dic)



