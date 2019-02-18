import numpy as np
import datetime as dt  # Python standard library datetime module
from mytools.met_tools import *
from mytools.netcdf_tools import * # ncdump implementation for python

def read_station_data(infile):
    '''
    Conversion of NASA aimes files from nilo.no data sets.
    '''
    import nappy as nap # Read and write NASA data files
    # Open file (read header only)
    nas_file = nap.openNAFile(infile)
    # Read actual data
    nas_file.readData()
    # Access the ozone data
    ozone_data_raw = np.array(nas_file.getNADict()['V'])
    # Close the nas-file
    nas_file.close()
    # Conversion
    #M_O = 15.9994                         # [g/mol] 
    #M_air = 28.949                        # [g/mol] 
    #mass_fraction = (3*M_O/M_air)
    #air_dens = np.array((1.4224, 1.1455)) # estimate for standard air density (-25 degC, 35 degC)
    #conversion = 1/np.mean(air_dens)/mass_fraction
    conversion = 1/2.
    nas_date = nas_file.getNADict()['DATE']
    start_date = dt.datetime.strptime("%s-0%s-0%s 00:00:00" % (nas_date[0], nas_date[1], nas_date[2]), '%Y-%m-%d %H:%M:%S')
    # Filter the data (0 - valid, >0 - invalid)
    x_time_station = np.ma.masked_where(ozone_data_raw[2]>0, 
                                        datetime_from_time(start_date, ozone_data_raw[0])[0])
    ozone_data = np.ma.masked_where(ozone_data_raw[2]>0, ozone_data_raw[1]*conversion)
    
    #yerr = ((ozone_data/np.mean(air_dens)-ozone_data/air_dens[0])/mass_fraction, 
    #        (ozone_data/air_dens[1]-ozone_data/np.mean(air_dens))/mass_fraction)
   
    data_dic = {'time':x_time_station, 'ozone':ozone_data} #'yerr':yerr
    return data_dic

def read_station_data_neumeyer(infile, **kargs):
    '''
    Conversion of Neumeyer station data. Pseudo NASA aimes format.
    '''
    import codecs
    conversion = kargs.pop('conversion', 1.9953)
    lines = codecs.open(infile, "r", 'utf-8').readlines()
    start_data = int(lines[0].split()[0])
    date = lines[6].split()
    time = []
    ozone_data_raw = []
    valid = []
    for each in lines[start_data+1:]:
        col = each.split()
        time.append((float(col[0])+float(col[1]))/2.)
        ozone_data_raw.append(float(col[2]))
        valid.append(float(col[-1]))

    start_date = dt.datetime.strptime("%s-%s-%s 00:00:00" % (date[0], date[1], date[2]), '%Y-%m-%d %H:%M:%S')
    x_time_station = np.ma.masked_where(np.array(valid)>0, 
                                        datetime_from_time(start_date, np.array(time))[0])
    ozone_data = np.ma.masked_where(np.array(valid)>0, np.array(ozone_data_raw)/conversion)
    data_dic = {'time':x_time_station, 'ozone':ozone_data}
    return data_dic

def read_station_data_noaa(infile, **kargs):
    '''
    Import of NOAA station data.
    **kargs utc - timeshift wrt UTC
            start_data - index where data block starts
            column - column in which data is found
    '''
    import codecs
    # Key words
    start_data = kargs.pop('start_data', 0)
    utc = kargs.pop('utc', 0)
    column = kargs.pop('column', 2)
    station = kargs.pop('station', 'standard')
    # Reading data
    lines = codecs.open(infile, "r", 'utf-8').read().splitlines()
    # Compute the timeshift
    timeshift = dt.timedelta(hours=utc)
    # Data structure
    time = []
    ozone_data_raw = []
    for each in lines[start_data:]:
        # Split the lines
        col = each.split()
        # Seperate ozone data
        ozone_data_raw.append(float(col[column]))
        # Generate timestamp
        if (station=='southpole') | (station=='other'):
            if  int(col[4]) == 24:
                temp_time = dt.datetime.strptime("%s-%s-%s %s" % (col[1], col[2], col[3], '00'), '%Y-%m-%d %H')+dt.timedelta(days=1)
            else:
                temp_time = dt.datetime.strptime("%s-%s-%s %s" % (col[1], col[2], col[3], col[4]), '%Y-%m-%d %H')
        else:
            if int(col[1][:2]) == 24:
                temp_time = dt.datetime.strptime("%s" % (col[0]+" 00:00"), '%Y-%m-%d %H:%M')+dt.timedelta(days=1)
            else:
                temp_time = dt.datetime.strptime("%s" % (col[0]+" "+col[1]), '%Y-%m-%d %H:%M')
        time.append(temp_time-timeshift)
    # Mask the data (negative ozone values)
    ozone_data = np.ma.masked_where(np.array(ozone_data_raw)<0, np.array(ozone_data_raw))
    x_time = np.ma.masked_where(np.array(ozone_data_raw)<0, np.array(time))
    # Output dictonary
    data_dic = {'time':x_time, 'ozone':ozone_data}
    return data_dic
