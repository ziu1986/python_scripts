import os, glob # Access environment variables
import numpy as np 
import nappy as nap # Read and write NASA data files
import matplotlib.pyplot as plt
import matplotlib.dates as mdates # handling of dates in plotting
import datetime as dt  # Python standard library datetime module
from mytools.med_tools import *
from mytools.netcdf_tools import * # ncdump implementation for python

# Standalone plot
b_stand = True
# Clean up
if b_stand:
    plt.close('all')

# Access environment variable for directory
nas_data = os.environ['DATA']
# Data directory
nas_subd = '/Ebas_Ozone'
nas_src = '/Zeppelin_Mountain/NO0042G.20000101000000.20130101000000.uv_abs.ozone.air.1y.1h.NO01L_uv_abs_uk_0042.NO01L_uv_abs..nas'
nas_src_3 = '/Neumeyer/DE0060G.20000101000000.20170201090710.uv_abs.ozone.air.1y.1h.DE06L_O3Neumayer2.DE06L_uv_ab.lev2.nas'

# Open file (read header only)
nas_file = nap.openNAFile(nas_data+nas_subd+nas_src_3)
# Read actual data
nas_file.readData()
# Access the ozone data
ozone_data_raw = np.array(nas_file.getNADict()['V'])
# Filter the data (0 - valid, >0 - invalid)
ozone_data = ozone_data_raw[:,np.where(ozone_data_raw[2]==0)[0]]
# Close the nas-file
nas_file.close()

station_name = nas_file.getNADict()['NCOM'][13][30:-14]
air_dens = np.array((1.4224, 1.1839)) # estimate for standard air density (-25 degC, 35 degC)
M_O = 15.9994    # [g/mol] 
M_air = 28.949   # [g/mol] 
mass_fraction = (3*M_O/M_air)
# Computations
nas_date = nas_file.getNADict()['DATE']
start_date = dt.datetime.strptime("%s-0%s-0%s 00:00:00" % (nas_date[0], nas_date[1], nas_date[2]), '%Y-%m-%d %H:%M:%S')
x_time_station, delta_time_station = datetime_from_time(start_date, ozone_data[0])
yerr = (ozone_data[1]/np.mean(air_dens)-ozone_data[1]/air_dens[0], 
        ozone_data[1]/air_dens[1]-ozone_data[1]/np.mean(air_dens))
if b_stand:
    # Plot it
    fig1 = plt.figure(1, figsize=(16,9))
    fig1.canvas.set_window_title("station_data_%s" % station_name)
    ax11 = plt.subplot()  
    ax11.set_title(station_name)
    ax11.errorbar(x_time_station, ozone_data[1]/np.mean(air_dens)/mass_fraction, yerr=np.array(yerr)/np.mean(air_dens)/mass_fraction, marker='x', color='red', ls='none', alpha=0.5)
    ax11.set_ylabel("Ozone (ppb)")
    #ax11.set_xlabel("Time")

    # Autoscale
    fig1.autofmt_xdate()

    plt.show(block=False)



