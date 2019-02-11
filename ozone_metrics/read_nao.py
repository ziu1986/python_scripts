import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt   # Plotting data
# Read the Northern Atlantic Oscillation index and display it
'''
The daily NAO index is constructed by projecting the daily (00Z) 500mb height anomalies over the Northern Hemsiphere onto the loading pattern of the NAO.  Please note that all year-round monthly mean anomaly data has been used to obtain the loading pattern of the NAO (Methodology).  Since the NAO has the largest variability during the cold season, the loading pattern primarily captures characteristics of the cold season NAO pattern.

The daily NAO index and its forecasts using MRF and Ensemble mean forecast data are shown for the previous 120 days as indicated and they are normalized by standard deviation of the monthly NAO index from 1950 to 2000.

'''

# Close the previous plots
plt.close('all')
# Directories of data
src_nao = os.environ['DATA']+'/processed_data/norm.nao.monthly.b5001.current.ascii'

def read_nao(infile):
    # Reading data
    lines = open(infile).readlines()
    data_raw = []
    time = []
    for each in lines:
        # Split the lines
        col = each.split()
        # Seperate ozone data
        data_raw.append(float(col[2]))
        time.append(dt.datetime.strptime("%s-%s" % (col[0], col[1]), '%Y-%m'))
    data_raw = pd.Series(data_raw,index=time)
    return(data_raw)

data = read_nao(src_nao)

# Plot the data
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot()
data_selection = data['2000':'2019']
data_selection.plot(color='black')
ax11.fill_between(data_selection.where(data_selection>0).index, data_selection.where(data_selection>=0), 0, color='red')
ax11.fill_between(data_selection.where(data_selection>0).index, data_selection.where(data_selection<=0), 0, color='blue')
ax11.set_xlabel('Time (years)')
ax11.set_ylabel('NAO index')
ax11.axhline(0, color='black', ls='--')
ax11.set_ylim(-4,4)
# Show it
plt.show(block=False)
    
