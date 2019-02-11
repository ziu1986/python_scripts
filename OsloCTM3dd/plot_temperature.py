import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *

execfile('growing_season.py')


def open_data(src_path):
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(src_path)):
        print(file)
        data_raw = xr.open_dataset(file)
        data_raw.coords['time'] = dt.datetime.strptime(str(data_raw['data_date'].data),"%Y%m%d%H")
        data_list.append(data_raw['T2M'])
              
    # Concatinating the list
    print("+++>Concatinating data_list<+++")
    data = xr.concat(data_list, dim='time')
    daily_ave = []
    for itime in range(len(data.time)/8):
        # Construct the daily averages
        daily_ave.append(data.isel(time=slice(itime*8,(itime+1)*8)).mean(dim='time'))
    daily_ave = xr.concat(daily_ave, dim='time')
    daily_ave.coords['time'] = data.time[4::8]

    return(data, daily_ave)
    
# Clean up
plt.close('all')
# Data source
nc_src = os.environ['DATA']+'/CTM3_input_data/metdataOpenIFS/cy38r1nc4_2005/T159N80L60/0[1,2,3,4,]/*'
nc_src_2 = os.environ['DATA']+'/CTM3_input_data/metdataOpenIFS/cy38r1nc4_2005/T159N80L60/0[5,6,7,8,]/*'
nc_src_3 = os.environ['DATA']+'/CTM3_input_data/metdataOpenIFS/cy38r1nc4_2005/T159N80L60/09/*'
nc_src_4 = os.environ['DATA']+'/CTM3_input_data/metdataOpenIFS/cy38r1nc4_2005/T159N80L60/1[0,1,2,]/*'

try:
    data
except NameError:
    data_1, data_ave_1 = open_data(nc_src)
    data_2, data_ave_2 = open_data(nc_src_2)
    data_3, data_ave_3 = open_data(nc_src_3)
    data_4, data_ave_4 = open_data(nc_src_4)
   
    data = xr.concat((data_1, data_2, data_3, data_4),dim='time')
    daily_ave = xr.concat((data_ave_1, data_ave_2, data_ave_3, data_ave_4),dim='time')

# Growing season
latitude = 64
longitude = (360-69)*320/360-1
test_location = (daily_ave-273.15).sel(lat=latitude, method='nearest').isel(lon=longitude)
test_location_2 = (daily_ave-273.15).sel(lat=-latitude, method='nearest').isel(lon=longitude)
# Classical
gs = growing_season(test_location)
sgs = gs['sgs']
egs = gs['egs']

#Sigmoid
s_shift = 182
growth = growing_season_sigmoid(test_location, criterium=(5,5))
sgs_2 = np.argmax(np.isclose(growth,[1],rtol=1e-5)==True)
sgs_2p = np.argmax(np.isclose(growth,[1],rtol=1e-7)==True)
sgs_2m = np.argmax(np.isclose(growth,[1],rtol=1e-3)==True)
fall = falling_season_sigmoid(test_location, criterium=(5,5),
                              s_shift=max(182,sgs_2), sgs=sgs_2)
egs_2 = np.argmax(np.isclose(fall,[0.5],rtol=1e-5)==True)
egs_2p = np.argmax(np.isclose(fall,[0.],rtol=1e-3)==True)
egs_2m = np.argmax(np.isclose(fall,[0.5],rtol=1e-3)==True)

#Fixed
sgs_fixed = int(start_growing_season_fixed(test_location.lat.data)[0])
egs_fixed = int(end_growing_season_fixed(test_location.lat.data)[0])

#Moving sum
#gs_ms = growing_season_moving(test_location, verbose=False)

# Plot it
fig1 = plt.figure(1, figsize=(10,8))
fig1.canvas.set_window_title("Relative_temperature_%s" % (latitude))
ax11 = plt.subplot()

(data-273.15).sel(lat=latitude, method='nearest').isel(lon=longitude).plot(label='3-hourly')
test_location.plot(ls='None', marker='x', label='daily avg')


ax11.axhline(5,color='grey',ls='--')
ax11.set_ylabel("T$_{2m}$ ($^\circ$C)")
ax11.set_xlabel("Time")
plt.legend(loc='upper left')

# Mark begin of growing season
text_ypos = ax11.get_ybound()[0]
# Sigmoid
if sgs_2 > 0:
    ax11.axvspan(daily_ave.time[sgs_2m].data,daily_ave.time[sgs_2p].data,color='orange',alpha=0.25)
    ax11.axvline(daily_ave.time[sgs_2].data, color='orange', linestyle=':', lw=2)
    ax11.text(daily_ave.time[sgs_2].data,
              text_ypos+3, "sigmoid sgs", color='orange', size=14)
if egs_2>0:
    ax11.axvspan(daily_ave.time[egs_2m].data,daily_ave.time[egs_2p].data,color='orange',alpha=0.25)
    ax11.axvline(daily_ave.time[egs_2].data, color='orange', linestyle=':', lw=2)
    ax11.text(daily_ave.time[egs_2].data,
              text_ypos+3, "sigmoid egs", color='orange', size=14)
else:
    ax11.axvline(daily_ave.time[-1].data, color='orange', linestyle=':', lw=2)
    ax11.text(daily_ave.time[-1].data,
              text_ypos+3, "sigmoid egs", color='orange', size=14)
#Classical
if sgs[1]:
    ax11.axvline(daily_ave.time[sgs[0]].data, color='darkcyan', linestyle=':', lw=2)
    ax11.text(daily_ave.time[sgs[0]].data, text_ypos+1, "classical sgs", color='darkcyan', size=14)
if egs[1]:
    ax11.axvline(daily_ave.time[egs[0]].data, color='darkcyan', linestyle=':', lw=2)
    ax11.text(daily_ave.time[egs[0]].data, text_ypos+1, "classical egs", color='darkcyan', size=14)
else:
    ax11.axvline(daily_ave.time[-1].data, color='darkcyan', linestyle=':', lw=2)
    ax11.text(daily_ave.time[-1].data, text_ypos+1, "classical egs", color='darkcyan', size=14)

#Fixed
if sgs_fixed > 0:
    ax11.axvline(daily_ave.time[sgs_fixed-1].data,
                 color='grey', linestyle=':', lw=2)
    ax11.text(daily_ave.time[sgs_fixed-1].data,
              text_ypos+2, "fixed sgs", color='grey', size=14)
if egs_fixed < 365:
    ax11.axvline(daily_ave.time[egs_fixed-1].data,
                 color='grey', linestyle=':', lw=2)
    ax11.text(daily_ave.time[egs_fixed-1].data,
              text_ypos+2, "fixed egs", color='grey', size=14)
else:
    ax11.axvline(daily_ave.time[-1].data,
                 color='grey', linestyle=':', lw=2)
    ax11.text(daily_ave.time[-1].data,
              text_ypos+2, "fixed egs", color='grey', size=14)

#Moving sum
#if gs_ms[0]:
#    ax11.axvline(daily_ave.time[gs_ms[0]].data,
#                 color='purple', linestyle=':', lw=2)
#    ax11.text(daily_ave.time[gs_ms[0]].data,
#              text_ypos+4, "ms sgs", color='purple', size=14)
#if gs_ms[1]:
#    ax11.axvline(daily_ave.time[gs_ms[1]].data,
#                 color='purple', linestyle=':', lw=2)
#    ax11.text(daily_ave.time[gs_ms[1]].data,
#              text_ypos+4, "ms egs", color='purple', size=14)
#else:
#    ax11.axvline(daily_ave.time[-1].data,
#                 color='purple', linestyle=':', lw=2)
#    ax11.text(daily_ave.time[-1].data,
#              text_ypos+4, "ms egs", color='purple', size=14)

    
fig2 = plt.figure(2, figsize=(10,8))
fig2.canvas.set_window_title("SGS_sigmoid_test_%s_%d" % (latitude,data.lon[longitude]))
ax21 = plt.subplot()

ax21.plot(daily_ave.time.data, growth)
ax21.axvline(daily_ave.time.data[sgs_2-1], color='blue', linestyle=':', lw=2)
ax21.set_xlabel("Time (days)")
ax21.set_ylabel("f(x) = 1/(1+exp(-(x-5)))")
ax21.plot(daily_ave.time.data, fall)
ax21.axvline(daily_ave.time.data[egs_2-1], color='red', linestyle=':', lw=2)


#fig3 = plt.figure(3)
#fig3.canvas.set_window_title("Sigmoid")
#ax31 = plt.subplot()
#x = np.linspace(-10,10)
#sig1 = 1/(1+np.exp(-2*x))
#sig2 = 1/(1+np.exp(-2*(x-5)))
#sig3 = np.exp(-2*(x-5))/(1+np.exp(-2*(x-5)))

#ax31.plot(x, sig1)
#ax31.plot(x, sig2)
#ax31.plot(x, sig3)

fig4 = plt.figure(4)
ax41 = plt.subplot()

gs_ms.plot()
ax41.axhline(25, color='purple',ls='--')
ax41.axvline(gs_ms.time.data[(np.where(gs_ms >= 25)[0])[0]], color='purple',ls='--')

fig5 = plt.figure(5)
fig5.canvas.set_window_title("Relative_temperature_%s" % (-latitude))
ax51 = plt.subplot()
(data-273.15).sel(lat=-latitude, method='nearest').isel(lon=longitude).plot(label='3-hourly')
test_location_2.plot(ls='None', marker='x', color='red', label='daily avg')
ax51.plot(test_location_2.time.data, test_location_2.roll(time=182).data, ls='None', marker='+', color='black', label='daily avg')

# Mark begin of growing season
text_ypos = ax51.get_ybound()[0]
#Fixed
sgs_fixed = int(start_growing_season_fixed(test_location_2.lat.data)[0])
egs_fixed = int(end_growing_season_fixed(test_location_2.lat.data)[1])
# Classical
gs = growing_season(test_location_2, sh=True)
sgs = gs['sgs']
egs = gs['egs']
#Fixed
if sgs_fixed > 0:
    ax51.axvline(daily_ave.time[sgs_fixed-1].data,
                 color='grey', linestyle=':', lw=2)
    ax51.text(daily_ave.time[sgs_fixed-1].data,
              text_ypos+2, "fixed sgs", color='grey', size=14)
if egs_fixed < 365:
    ax51.axvline(daily_ave.time[egs_fixed-1].data,
                 color='grey', linestyle=':', lw=2)
    ax51.text(daily_ave.time[egs_fixed-1].data,
              text_ypos+2, "fixed egs", color='grey', size=14)
else:
    ax51.axvline(daily_ave.time[-1].data,
                 color='grey', linestyle=':', lw=2)
    ax51.text(daily_ave.time[-1].data,
              text_ypos+2, "fixed egs", color='grey', size=14)

#Classical
if sgs[1]:
    ax51.axvline(daily_ave.time[sgs[0]].data, color='darkcyan', linestyle=':', lw=2)
    ax51.text(daily_ave.time[sgs[0]].data, text_ypos+1, "classical sgs", color='darkcyan', size=14)
if egs[1]:
    ax51.axvline(daily_ave.time[egs[0]].data, color='darkcyan', linestyle=':', lw=2)
    ax51.text(daily_ave.time[egs[0]].data, text_ypos+1, "classical egs", color='darkcyan', size=14)
else:
    ax51.axvline(daily_ave.time[-1].data, color='darkcyan', linestyle=':', lw=2)
    ax51.text(daily_ave.time[-1].data, text_ypos+1, "classical egs", color='darkcyan', size=14)

ax51.axhline(5,color='grey',ls='--')
ax51.set_ylabel("T$_{2m}$ ($^\circ$C)")
ax51.set_xlabel("Time")
plt.legend(loc='upper left')
# Show it
plt.show(block=False)
