import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
sys.path.append("../mytools/") 
from med_tools import *

# Data
src = 'BrO_col_2000*_BrXplo_mysic.nc'
src_ref = 'BrO_col_2000*_BrXplo_ref.nc'
month = np.arange(1,13)
month_name = ("January", "February", "March", "April", "May", "June", 
              "July", "August", "September", "October", "November", "December")
# Read data
data_list = []
for file in sorted(glob.glob(src)):
    data = xr.open_dataset(file)
    data_list.append(data)
data_list = xr.concat(data_list, dim="month")

data_list_ref = []
for file in sorted(glob.glob(src_ref)):
    data = xr.open_dataset(file)
    data_list_ref.append(data)
data_list_ref = xr.concat(data_list_ref, dim="month")    
# Clean up
plt.close('all')
fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("EMAC_BrO_tot_zonal_mean")
ax21 = plt.subplot(221)
ax21.set_title("JFM")
for imonth in (1,2,3):
    ax21.plot(data_list.lat, data_list.sel(month=imonth).BrO_column*1e-13, label=month_name[imonth-1])
    ax21.plot(data_list_ref.lat, data_list_ref.sel(month=imonth).BrO_column*1e-13,
              ls=':', color=ax21.lines[-1].get_color(), label='')
ax22 = plt.subplot(222)
ax22.set_title("AMJ")
for imonth in (4,5,6):
    ax22.plot(data_list.lat, data_list.sel(month=imonth).BrO_column*1e-13, label=month_name[imonth-1])
    ax22.plot(data_list_ref.lat, data_list_ref.sel(month=imonth).BrO_column*1e-13,
              ls=':', color=ax22.lines[-1].get_color(), label='')
ax23 = plt.subplot(223)
ax23.set_title("JAS")
for imonth in (7,8,9):
    ax23.plot(data_list.lat, data_list.sel(month=imonth).BrO_column*1e-13, label=month_name[imonth-1])
    ax23.plot(data_list_ref.lat, data_list_ref.sel(month=imonth).BrO_column*1e-13,
              ls=':', color=ax23.lines[-1].get_color(), label='')
ax24 = plt.subplot(224)
ax24.set_title("OND")
for imonth in (10,11,12):
    ax24.plot(data_list.lat, data_list.sel(month=imonth).BrO_column*1e-13, label=month_name[imonth-1])
    ax24.plot(data_list_ref.lat, data_list_ref.sel(month=imonth).BrO_column*1e-13,
              ls=':', color=ax24.lines[-1].get_color(), label='')
for ax in fig2.axes:
    ax.set_xticks(np.arange(-90, 91, 15))
    ax.set_xlim(-90,90)
    ax.set_ylim(0,7)
    ax.axvline(-45, color='grey', ls='--')
    ax.axvline(45, color='grey', ls='--')
    ax.legend(frameon=False,loc='upper right')
ax24.set_xlabel("Latitude (deg)", x=-0.1)
ax23.set_ylabel("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'), y=1.1)

ax21.text(-75, 6.5, "BrXplo_ref")
ax21.plot((-86,-78), (6.6,6.6), ls=':', color='grey')
ax21.text(-75, 6, "BrXplo_mysic")
ax21.plot((-86,-78), (6.1,6.1), color='grey')
plt.show(block=False)
