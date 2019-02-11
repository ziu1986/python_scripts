import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Data source
nc_src_corr = os.environ['DATA']+'/abel/C3RUN_nitrate_corr_2005/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src = os.environ['DATA']+'/abel/C3RUN_test_ncorr.010318.r206/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src_r212 = os.environ['DATA']+'/abel/C3RUN_test.010318.r212/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src_r213 = os.environ['DATA']+'/abel/C3RUN_test.010318.r213/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src_r214 = os.environ['DATA']+'/abel/C3RUN_test.010318.r214_meanVegH/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src_rc1 = os.environ['DATA']+'/abel/C3RUN_test.020318.9107/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src_rc2 = os.environ['DATA']+'/abel/C3RUN_test.020318.30931/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src_rc3 = os.environ['DATA']+'/abel/C3RUN_test.050318.10547/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src_rc4 = os.environ['DATA']+'/abel/C3RUN_test.050318.28205/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src_rc5 =  os.environ['DATA']+'/abel/C3RUN_test.060318.23060/scavenging_daily/scavenging_daily_2d_20050201.nc'
nc_src_emep =  os.environ['DATA']+'/abel/C3RUN_emep.080318/scavenging_daily_2d_20050101.nc'
nc_src_corr =  os.environ['DATA']+'/abel/C3RUN_test.090318.31846/scavenging_daily/scavenging_daily_2d_20050101.nc'
#nc_src_corr2 =  os.environ['DATA']+'/abel/C3RUN_test.090318.25072/scavenging_daily/scavenging_daily_2d_20050101.nc'
nc_src_corr2 =  os.environ['DATA']+'/abel/C3RUN_emep.090318.9433/scavenging_daily/scavenging_daily_2d_20050101.nc'
#labels = ('corr','r206','r212','r213','r214','wc1','wc_wrong_gsto','wc_condXres','wc4_PAR','growseason')
labels = ('base', 'emep', 'corr', 'LAI')
try:
    data
except NameError:
    data_list = []
    for subdir in (nc_src_corr,nc_src_emep,nc_src_corr,nc_src_corr2):#(nc_src_corr,nc_src,nc_src_r212,nc_src_r213,nc_src_r214,nc_src_rc1,nc_src_rc2,nc_src_rc3,nc_src_rc4,nc_src_rc5):
        print("Reading from path %s" % (os.path.abspath(subdir)))
        # Open dataset
        data = xr.open_dataset(subdir)
        data_list.append(data)

# Extract ozone drydeposition and rescale units
ozone_data = [data['dry_O3']/mega for data in data_list]
for data in ozone_data:
    data.attrs['unit'] = 'Gg'

# Plot it
plt.close('all')
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(141)
ax12 = plt.subplot(142)
ax13 = plt.subplot(143)
ax14 = plt.subplot(144)
#ax15 = plt.subplot(335)
#ax16 = plt.subplot(336)
#ax17 = plt.subplot(337)
#ax18 = plt.subplot(338)
#(ozone_data[1]-ozone_data[0]).plot(ax=ax11, vmin=-5., vmax=5., cmap=plt.cm.seismic)
#(ozone_data[2]-ozone_data[0]).plot(ax=ax11, cmap=plt.cm.seismic)
(ozone_data[0]).plot(ax=ax11, vmin=0., vmax=5.)
(ozone_data[1]-ozone_data[0]).plot(ax=ax12, vmin=-.05, vmax=.05, cmap=plt.cm.seismic)
(ozone_data[2]-ozone_data[0]).plot(ax=ax13, vmin=-.05, vmax=.05, cmap=plt.cm.seismic)
(ozone_data[3]-ozone_data[0]).plot(ax=ax14, vmin=-.05, vmax=.05, cmap=plt.cm.seismic)
#(ozone_data[6]-ozone_data[0]).plot(ax=ax16, vmin=-2., vmax=2., cmap=plt.cm.seismic)
#(ozone_data[7]-ozone_data[5]).plot(ax=ax15, cmap=plt.cm.seismic)
#(ozone_data[8]-ozone_data[7]).plot(ax=ax16, cmap=plt.cm.seismic)
#(ozone_data[9]-ozone_data[8]).plot(ax=ax17, cmap=plt.cm.seismic)

fig2 = plt.figure(2, figsize=(16,9))
ax21 = plt.subplot(211)
for i in np.arange(len(ozone_data)):
    ozone_data[i].mean(dim='lon').plot(ax=ax21, label=labels[i])
ax21.legend()
ax22 = plt.subplot(212)
for i in np.arange(1,len(ozone_data)):
    (ozone_data[i]-ozone_data[0]).mean(dim='lon').plot(ax=ax22, label=labels[i])
ax22.legend()
# Show it
plt.show(block=False)
