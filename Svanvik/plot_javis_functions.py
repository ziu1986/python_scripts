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

# Set up the different species
# Evergreen
evergreen = JavisModel('evergreen', Tmin=0, Tmax=200, Topt=20, fmin=0.1, Dmax=0.8, Dmin=2.8, alpha=0.006, gmax=125)
# Birch
birch = JavisModel('birch', Tmin=5, Tmax=200, Topt=20, fmin=0.1, Dmax=0.5, Dmin=2.7, alpha=0.0042, gmax=240)
# Grassland
grassland = JavisModel('grassland', Tmin=10, Tmax=36, Topt=24, fmin=0.1, Dmax=1.75, Dmin=4.5, alpha=0.011, gmax=190)


# Load data
data_temp = import_data(src_svanvik)
data_rad = import_data(src_svanvik_rad)

# Clean up
plt.close('all')


# Loop through species
for species, i in zip((evergreen, birch, grassland), np.arange(1,4)):
                      
    f_temp = data_temp.iloc[:,0].apply(lambda x: species.f_temp(x))

    vpd = VPD(data_temp.iloc[:,1], data_temp.iloc[:,0])/kilo
    f_vpd = vpd.apply(lambda x: species.f_vpd(x))

    f_light = data_rad.iloc[:,0].apply(lambda x: species.f_light(x))

    # Plot it
    fig = plt.figure(i, figsize=(10,12))
    fig.canvas.set_window_title("javis_funcs_%s" % species.name)
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)

    f_temp['2019-05':'2019-08'].plot(ax=ax1)
    f_vpd['2019-05':'2019-08'].plot(ax=ax2)
    f_light['2019-05':'2019-08'].plot(ax=ax3)

    ax1.set_ylabel("f_temp")
    ax2.set_ylabel("f_vpd")
    ax3.set_ylabel("f_light")
    ax3.set_xlabel("Time (months)")

    for ax in fig.axes:
        ax.set_ylim(0,1)

# Show it
plt.show(block=False)
