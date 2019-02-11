import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation

# Read the data
#nc_src = os.environ['DATA']+'/SN_d_tot_V2.0.txt'

def read_sunspots(file_name):
    '''
    Function to read sunspot data.
    Returns a pandas DataFrame object.
    '''
    # Open the file and read all lines at once
    data_sun_spots = open(file_name).readlines()
    date_list = []
    data_raw_list = []
    # Loop through the lines
    for line in data_sun_spots:
        data = line.split()
        date_list.append(dt.datetime(int(data[0]), int(data[1]), int(data[2])))
        data_raw_list.append((int(data[4]), float(data[5])))
    # Set up the DataFrame
    data_sunspots = pd.DataFrame(data_raw_list, index=date_list, columns=['Ntot', 'std'])
    # Return the masked values
    return(data_sunspots.where(data_sunspots[['Ntot']]!=-1))
