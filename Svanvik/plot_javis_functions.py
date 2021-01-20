import os, glob, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
from scipy.constants import *     # Get physics constants
from javis_model.javis_model import *
from mytools.plot_tools import *
from mytools.ozone_tools import VPD

# Data sources
src_svanvik = os.environ['DATA']+'/astra_data/observations/metdata_svanvik/Svanvik_temp_relhum_wind_*.csv'
src_svanvik_rad = os.environ['DATA']+'/astra_data/observations/metdata_svanvik/svanvik_glob_rad_*.csv'


def import_data(src):
    data_list = []
    # Read data
    for file in sorted(glob.glob(src)):
        print(file)
        data_svanvik = pd.read_csv(file)
        data_svanvik.index = pd.date_range("%s" % data_svanvik['Time measured'][0][:-3], "%s" % data_svanvik['Time measured'].iloc[-1][:-3], freq='H')
        data_svanvik = data_svanvik.drop(columns=["Time measured"])

        data_list.append(data_svanvik)
        
    data_svanvik = pd.concat(data_list)

    return(data_svanvik)

#data_svanvik_temp = (data_svanvik['Temp_degC'].where((data_svanvik['Temp_degC']>-50) & (data_svanvik['Temp_degC']<50)).dropna()
# main

# Evergreen
#Tmin=0, Tmax=200, Topt=20, fmin=0.1, Dmax=0.8, Dmin=2.8, alpha=0.006
# Birch
#Tmin=5, Tmax=200, Topt=20, fmin=0.1, Dmax=0.5, Dmin=2.7, alpha=0.0042
# Grassland
#Tmin=10, Tmax=36, Topt=24, fmin=0.1, Dmax=1.75, Dmin=4.5, alpha=0.011

data_temp = import_data(src_svanvik)
data_rad = import_data(src_svanvik_rad)

f_temp = data_temp.iloc[:,0].apply(lambda x: f_temp(x, Tmin=0, Tmax=200, Topt=20))

vpd = VPD(data_temp.iloc[:,1], data_temp.iloc[:,0])/kilo
f_vpd = vpd.apply(lambda x: f_vpd(x, fmin=0.1, Dmin=2.8,Dmax=0.8))

f_light = data_rad.iloc[:,0].apply(lambda x: f_light(x, alpha=0.006))

# Plot it
# Clean up
plt.close('all')

fig1 = plt.figure(1, figsize=(10,12))
fig1.canvas.set_window_title("javis_funcs_evergreen")
ax11 = plt.subplot(311)
ax12 = plt.subplot(312)
ax13 = plt.subplot(313)

f_temp['2019-05':'2019-08'].plot(ax=ax11)
f_vpd['2019-05':'2019-08'].plot(ax=ax12)
f_light['2019-05':'2019-08'].plot(ax=ax13)

ax11.set_ylabel("f_temp")
ax12.set_ylabel("f_vpd")
ax13.set_ylabel("f_light")
ax13.set_xlabel("Time (months)")

for ax in fig1.axes:
    ax.set_ylim(0,1)

# Show it
plt.show(block=False)
