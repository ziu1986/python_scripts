import os, glob # Access environment variables
import numpy as np
import pandas as pd
from scipy.constants import * # Get physics constants
from scipy import stats # Linear regressions etc.
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt
import matplotlib.dates as mdates # handling of dates in plotting
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime as dt  # Python standard library datetime module
from mytools.med_tools import *
from mytools.netcdf_tools import * # ncdump implementation for python

#from read_station_data import * # Various data formats :-|
execfile('read_station_data.py')
# Clean up
plt.close('all')
b_stat = False
# Access environment variable for directory
nc_data = os.environ['DATA']
# Data directory
nc_subd = '/BrXplo'
nc_src = '/BrXplo_fysic/BrXplo*2000*tr_Ox_HOx.nc'
nc_src_ref = '/BrXplo_ref/BrXplo_ref*2000*tr_Ox_HOx.nc'
nc_src_mysic = '/BrXplo_mysic/BrXplo_mysic*2000*tr_Ox_HOx.nc'
nc_src_mysic_rsnow = '/BrXplo_mysic_rsnow/BrXplo_mysic*2000*tr_Ox_HOx.nc'
# Surface Ozone data
nas_subd = '/Ebas_Ozone'
nas_src = '/Zeppelin_Mountain/NO0042G.20000101000000.20130101000000.uv_abs.ozone.air.1y.1h.NO01L_uv_abs_uk_0042.NO01L_uv_abs..nas'
nas_src_2 = '/Alert/CA0420G.20000101000000.20140501000000.uv_abs.ozone.air.1y.1h.CA01L_O3mon_420.CA01L_UV_absorption..nas'
nas_src_3 = '/Neumeyer/DE0060G.20000101000000.20170201090710.uv_abs.ozone.air.1y.1h.DE06L_O3Neumayer2.DE06L_uv_ab.lev2.nas'
barrow_src = '/Barrow/BRW_Ozone_Hourly_2000'
summit_src = '/Summit/sum_ozone_hourly_2000.dat'
southpole_src = '/South_Pole/sptcl2000'

# Steering of x-axis
xShift = (dt.date(2000, 1, 1), dt.date(2001, 1, 1))
xShift2 = xShift#(dt.date(2000, 8, 1), dt.date(2000, 12, 1))#
b_rs = True
try:
    data_ref
except NameError:
    print("Reading data...")
    print("AirSnow...")
    print("Reference...")
    data_ref = fetch_data(nc_data+nc_subd+nc_src_ref, 'O3', level=-1)
    print("Multi-year SIC...")
    data_mysic = fetch_data(nc_data+nc_subd+nc_src_mysic, 'O3', level=-1, coordinates=False)
    print("Multi-year SIC rsnow...")
    data_mysic_rs = fetch_data(nc_data+nc_subd+nc_src_mysic_rsnow, 'O3', level=-1, coordinates=False)
    print("First-year SIC...")
    data = fetch_data(nc_data+nc_subd+nc_src, 'O3', level=-1, coordinates=False)
    
    print("Read station data...")
    data_zeppelin = read_station_data(nc_data+nas_subd+nas_src)
    data_alert = read_station_data(nc_data+nas_subd+nas_src_2)
    data_neumeyer = read_station_data_neumeyer(nc_data+nas_subd+nas_src_3)
    data_barrow = read_station_data_noaa(nc_data+nas_subd+barrow_src, start_data=26, utc=-9)
    data_summit = read_station_data_noaa(nc_data+nas_subd+summit_src, start_data=32, column=4)
    data_southpole = read_station_data_noaa(nc_data+nas_subd+southpole_src, start_data=0, column=5, utc=-12, station='southpole')
# Computations
start_date = dt.datetime.strptime(data_ref[1]['time'][1][10:], '%Y-%m-%d %H:%M:%S')
x_time, delta_time = datetime_from_time(start_date, data_ref[0]['time'])

unfin = len(data[0]['O3'])
# Plot it
fig1 = plt.figure(1, figsize=(16,12))
if b_rs:
    fig1.canvas.set_window_title("rs_surface_ozone")
else:
    fig1.canvas.set_window_title("surface_ozone")
ax11 = plt.subplot(411)
ax11.set_title("Alert (Canada): 82.50 N, 62.30 W, 210 m a.s.l.")
ax11.errorbar(data_alert['time'].compressed(), data_alert['ozone'].compressed(),
              yerr=(data_alert['yerr'][0].compressed(),data_alert['yerr'][1].compressed()),
              zorder=2, marker='x', color='red', ls='none', alpha=0.75, label='obs')
if not b_rs:
    ax11.plot(x_time, data_ref[0]['O3'][:,find_nearest(82.5, data_ref[0]['lat'], index=True),
                                        find_nearest(360-62.3, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='black')
    ax11.plot(x_time, data[0]['O3'][:,find_nearest(82.5, data_ref[0]['lat'], index=True),
                                    find_nearest(360-62.3, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='cornflowerblue', ls='--')
ax11.plot(x_time, data_mysic[0]['O3'][:,find_nearest(82.5, data_ref[0]['lat'], index=True),
                                      find_nearest(360-62.3, data_ref[0]['lon'], index=True)]*1e9,
          zorder=3, color='blue')
if b_rs:
    ax11.plot(x_time, data_mysic_rs[0]['O3'][:,find_nearest(82.5, data_ref[0]['lat'], index=True),
                                             find_nearest(360-62.3, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='c', ls=':')


ax12 = plt.subplot(412, sharex=ax11)
ax12.set_title('Barrow (Alaska): 71.32 N, 156.61 W, 8 m a.s.l.')
ax12.errorbar(data_barrow['time'], data_barrow['ozone'], 
              zorder=2, marker='x', ls='none', color='red', alpha=0.75)
if not b_rs:
    ax12.plot(x_time, data_ref[0]['O3'][:,find_nearest(71.32, data_ref[0]['lat'], index=True),
                                        find_nearest(360-156.61, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='black')
    ax12.plot(x_time, data[0]['O3'][:,find_nearest(71.32, data_ref[0]['lat'], index=True),
                                    find_nearest(360-156.61, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='cornflowerblue', ls='--')
ax12.plot(x_time, data_mysic[0]['O3'][:,find_nearest(71.32, data_ref[0]['lat'], index=True),
                                      find_nearest(360-156.61, data_ref[0]['lon'], index=True)]*1e9,
          zorder=3, color='blue')
if b_rs:
    ax12.plot(x_time, data_mysic_rs[0]['O3'][:,find_nearest(71.32, data_ref[0]['lat'], index=True),
                                             find_nearest(360-156.61, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='c', ls=':')

ax13 = plt.subplot(413, sharex=ax11)
ax13.set_title('Zeppelin Mountain (Norway): 78.90 N, 11.88 E, 474 m a.s.l.')
ax13.errorbar(data_zeppelin['time'].compressed(), data_zeppelin['ozone'].compressed(),
              yerr=(data_zeppelin['yerr'][0].compressed(),data_zeppelin['yerr'][1].compressed()),
              zorder=2, marker='x', color='red', ls='none', alpha=0.75, label='obs')
if not b_rs:
    ax13.plot(x_time, data_ref[0]['O3'][:,find_nearest(78.9, data_ref[0]['lat'], index=True),
                                        find_nearest(11.88, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='black')
    ax13.plot(x_time, data[0]['O3'][:,find_nearest(78.9, data_ref[0]['lat'], index=True),
                                    find_nearest(11.88, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='cornflowerblue', ls='--')
ax13.plot(x_time, data_mysic[0]['O3'][:,find_nearest(78.9, data_ref[0]['lat'], index=True),
                                       find_nearest(11.88, data_ref[0]['lon'], index=True)]*1e9,
          zorder=3, color='blue')
if b_rs:
    ax13.plot(x_time, data_mysic_rs[0]['O3'][:,find_nearest(78.9, data_ref[0]['lat'], index=True),
                                             find_nearest(11.88, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='c', ls=':')

ax14 = plt.subplot(414, sharex=ax11)
ax14.set_title('Summit (Greenland): 72.54 N, 38.48 W, 3238 m a.s.l.')
ax14.errorbar(data_summit['time'], data_summit['ozone'], 
              zorder=2, marker='x', ls='none', color='red', alpha=0.75, label='obs')
if not b_rs:
    ax14.plot(x_time, data_ref[0]['O3'][:,find_nearest(72.57, data_ref[0]['lat'], index=True),
                                        find_nearest(360-38.48, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='black', label='ref')
    ax14.plot(x_time, data[0]['O3'][:,find_nearest(72.57, data_ref[0]['lat'], index=True),
                                    find_nearest(360-38.48, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='cornflowerblue', ls='--', label='fysic')
ax14.plot(x_time, data_mysic[0]['O3'][:,find_nearest(72.57, data_ref[0]['lat'], index=True),
                                    find_nearest(360-38.48, data_ref[0]['lon'], index=True)]*1e9,
          zorder=3, color='blue', label='mysic')
if b_rs:
        ax14.plot(x_time, data_mysic_rs[0]['O3'][:,find_nearest(72.57, data_ref[0]['lat'], index=True),
                                                 find_nearest(360-38.48, data_ref[0]['lon'], index=True)]*1e9,
                  zorder=3, color='c', ls=':', label='mysic_rs')

if b_rs:
    ax14.legend(frameon=False, ncol=5, loc='best')
else:
    ax14.legend(frameon=False, ncol=4, loc='best')

for ax in fig1.axes:
    ax.set_ylim(0,60)
    ax.set_xlim(xShift)
ax13.set_ylabel("%s (%s)" % ('O$_3$', 'ppb'), y=1.25)

fig1.autofmt_xdate()

fig2 = plt.figure(2, figsize=(16,12))
if b_rs:
    fig2.canvas.set_window_title("rs_surface_ozone_sh")
else:
    fig2.canvas.set_window_title("surface_ozone_sh") #_08-11
ax21 = plt.subplot(411)
ax21.set_title("Palmer Station (Antarctica, USA): 64.77 S, 64.05 W, 21 m a.s.l.")
if not b_rs:
    ax21.plot(x_time, data_ref[0]['O3'][:,find_nearest(-64.77, data_ref[0]['lat'], index=True),
                                        find_nearest(360-64.05, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='black')
    ax21.plot(x_time, data[0]['O3'][:,find_nearest(-64.77, data_ref[0]['lat'], index=True),
                                    find_nearest(360-64.05, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='cornflowerblue', ls='--')
ax21.plot(x_time, data_mysic[0]['O3'][:,find_nearest(-64.77, data_ref[0]['lat'], index=True),
                                      find_nearest(360-64.05, data_ref[0]['lon'], index=True)]*1e9,
          zorder=3, color='blue')
if b_rs:
    ax21.plot(x_time, data_mysic_rs[0]['O3'][:,find_nearest(-64.77, data_ref[0]['lat'], index=True),
                                             find_nearest(360-64.05, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='c', ls=':')

ax22 = plt.subplot(412, sharex=ax21)
ax22.set_title('Neumeyer Station (Antarctica, Germany): 70.68 S, 8.26 W, 43 m a.s.l.')
ax22.errorbar(data_neumeyer['time'], data_neumeyer['ozone'], 
              zorder=2, marker='x', ls='none', color='red', alpha=0.75)
if not b_rs:
    ax22.plot(x_time, data_ref[0]['O3'][:,find_nearest(-70.68, data_ref[0]['lat'], index=True),
                                        find_nearest(360-8.26, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='black')
    ax22.plot(x_time, data[0]['O3'][:,find_nearest(-70.68, data_ref[0]['lat'], index=True),
                                    find_nearest(360-8.26, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='cornflowerblue', ls='--')
ax22.plot(x_time, data_mysic[0]['O3'][:,find_nearest(-70.68, data_ref[0]['lat'], index=True),
                                      find_nearest(360-8.26, data_ref[0]['lon'], index=True)]*1e9,
          zorder=3, color='blue')
if b_rs:
    ax22.plot(x_time, data_mysic_rs[0]['O3'][:,find_nearest(-70.68, data_ref[0]['lat'], index=True),
                                             find_nearest(360-8.26, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='c', ls=':')

ax23 = plt.subplot(413, sharex=ax21)
ax23.set_title('Arrival Heights (Antarctica, New Zealand): 77.85 S, 166.78 E, 22 m a.s.l.')
if not b_rs:
    ax23.plot(x_time, data_ref[0]['O3'][:,find_nearest(-77.85, data_ref[0]['lat'], index=True),
                                        find_nearest(166.78, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='black')
    ax23.plot(x_time, data[0]['O3'][:,find_nearest(-77.85, data_ref[0]['lat'], index=True),
                                    find_nearest(166.78, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='cornflowerblue', ls='--')
ax23.plot(x_time, data_mysic[0]['O3'][:,find_nearest(-77.85, data_ref[0]['lat'], index=True),
                                    find_nearest(166.78, data_ref[0]['lon'], index=True)]*1e9,
          zorder=3, color='blue')
if b_rs:
    ax23.plot(x_time, data_mysic_rs[0]['O3'][:,find_nearest(-77.85, data_ref[0]['lat'], index=True),
                                             find_nearest(166.78, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='c', ls=':')

ax24 = plt.subplot(414, sharex=ax21)
ax24.set_title('South Pole Station (Antarctica, USA): 89.98 S, 24.8 W, 2810 m a.s.l.')
ax24.errorbar(data_southpole['time'], data_southpole['ozone'], 
              zorder=2, marker='x', ls='none', color='red', alpha=0.75, label='obs')
if not b_rs:
    ax24.plot(x_time, data_ref[0]['O3'][:,find_nearest(-90, data_ref[0]['lat'], index=True),
                                        find_nearest(360-24.8, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='black', label='ref')
    ax24.plot(x_time, data[0]['O3'][:,find_nearest(-90, data_ref[0]['lat'], index=True),
                                    find_nearest(360-24.8, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='cornflowerblue', ls='--', label='fysic')
ax24.plot(x_time, data_mysic[0]['O3'][:,find_nearest(-90, data_ref[0]['lat'], index=True),
                                      find_nearest(360-24.8, data_ref[0]['lon'], index=True)]*1e9,
          zorder=3, color='blue', label='mysic')
if b_rs:
    ax24.plot(x_time, data_mysic_rs[0]['O3'][:,find_nearest(-90, data_ref[0]['lat'], index=True),
                                             find_nearest(360-24.8, data_ref[0]['lon'], index=True)]*1e9,
              zorder=3, color='c', ls=':', label='mysic_rs')

if b_rs:
    ax24.legend(frameon=False, ncol=5, loc='best')
else:
    ax24.legend(frameon=False, ncol=4, loc='best')

for ax in fig2.axes:
    ax.set_ylim(0,60)
    ax.set_xlim(xShift2)
ax23.set_ylabel("%s (%s)" % ('O$_3$', 'ppb'), y=1.25)

fig2.autofmt_xdate()

fig3 = plt.figure(3)
fig3.canvas.set_window_title("surface_ozone_corr_mysic_rsnow")
ax31 = plt.subplot(221)
ax31.set_title("Alert", ha='left', y=0.85, x=0.02)
ax31.hist2d(data_alert['ozone'][:-24],
            data_mysic_rs[0]['O3'][:,find_nearest(82.5, data_ref[0]['lat'], index=True),
                                find_nearest(360-62.3, data_ref[0]['lon'], index=True)]*1e9,
            bins=40, range=((0,40),(0,40)), cmap=plt.cm.Reds)
ax31.set_ylabel("Ozone EMAC mysic_rs (ppb)", y=-0.25)
ax32 = plt.subplot(222)
ax32.set_title("Barrow", ha='left', y=0.85, x=0.02)
ax32.hist2d(data_barrow['ozone'][:-24],
            data_mysic_rs[0]['O3'][:,find_nearest(71.32, data_ref[0]['lat'], index=True),
                                      find_nearest(360-156.61, data_ref[0]['lon'], index=True)]*1e9,
            bins=40, range=((0,40),(0,40)), cmap=plt.cm.Reds)
ax33 = plt.subplot(223)
ax33.set_title("Zeppelin Mountain", ha='left', y=0.85, x=0.02)
ax33.hist2d(data_zeppelin['ozone'][:-24],
            data_mysic_rs[0]['O3'][:,find_nearest(78.9, data_ref[0]['lat'], index=True),
                                       find_nearest(11.88, data_ref[0]['lon'], index=True)]*1e9,
            bins=40, range=((0,40),(0,40)), cmap=plt.cm.Reds)
ax34 = plt.subplot(224)
ax34.set_title("Summit", ha='left', y=0.85, x=0.02)
ax34.hist2d(data_summit['ozone'][1776:-24],
            data_mysic_rs[0]['O3'][24*30*7+24*16:,find_nearest(72.57, data_ref[0]['lat'], index=True),
                                    find_nearest(360-38.48, data_ref[0]['lon'], index=True)]*1e9,
            bins=40, range=((0,40),(0,40)), cmap=plt.cm.Reds)
ax34.set_xlabel("Ozone in situ observation (ppb)", x=-0.25)

for ax in fig3.axes:
    ax.plot((0,40),(0,40), ls='--', color='black')

fig4 = plt.figure(4)
fig4.canvas.set_window_title("surface_ozone_corr_mysic")
ax41 = plt.subplot(221)
ax41.set_title("Alert", ha='left', y=0.85, x=0.02)
ax41.hist2d(data_alert['ozone'][:-24],
            data_mysic[0]['O3'][:,find_nearest(82.5, data_ref[0]['lat'], index=True),
                                find_nearest(360-62.3, data_ref[0]['lon'], index=True)]*1e9,
            bins=40, range=((0,40),(0,40)), cmap=plt.cm.Reds)
ax41.set_ylabel("Ozone EMAC mysic (ppb)", y=-0.25)
ax42 = plt.subplot(222)
ax42.set_title("Barrow", ha='left', y=0.85, x=0.02)
ax42.hist2d(data_barrow['ozone'][:-24],
            data_mysic[0]['O3'][:,find_nearest(71.32, data_ref[0]['lat'], index=True),
                                      find_nearest(360-156.61, data_ref[0]['lon'], index=True)]*1e9,
            bins=40, range=((0,40),(0,40)), cmap=plt.cm.Reds)
ax43 = plt.subplot(223)
ax43.set_title("Zeppelin Mountain", ha='left', y=0.85, x=0.02)
ax43.hist2d(data_zeppelin['ozone'][:-24],
            data_mysic[0]['O3'][:,find_nearest(78.9, data_ref[0]['lat'], index=True),
                                       find_nearest(11.88, data_ref[0]['lon'], index=True)]*1e9,
            bins=40, range=((0,40),(0,40)), cmap=plt.cm.Reds)
ax44 = plt.subplot(224)
ax44.set_title("Summit", ha='left', y=0.85, x=0.02)
ax44.hist2d(data_summit['ozone'][1776:-24],
            data_mysic[0]['O3'][24*30*7+24*16:,find_nearest(72.57, data_ref[0]['lat'], index=True),
                                    find_nearest(360-38.48, data_ref[0]['lon'], index=True)]*1e9,
            bins=40, range=((0,40),(0,40)), cmap=plt.cm.Reds)
ax44.set_xlabel("Ozone in situ observation (ppb)", x=-0.25)

for ax in fig4.axes:
    ax.plot((0,40),(0,40), ls='--', color='black')

if b_stat:
    fig4 = plt.figure(4)
    fig4.canvas.set_window_title("surface_ozone_dist")
    ax41 = plt.subplot(221)
    ax41.set_title("Alert", ha='left', y=0.85, x=0.02)
    ax41.hist(data_mysic[0]['O3'][:,find_nearest(82.5, data_ref[0]['lat'], index=True),
                                  find_nearest(360-62.3, data_ref[0]['lon'], index=True)]*1e9-data_alert['ozone'][:-24],
              bins=40, range=(-40,40))
    ax41.hist(data[0]['O3'][:,find_nearest(82.5, data_ref[0]['lat'], index=True),
                            find_nearest(360-62.3, data_ref[0]['lon'], index=True)]*1e9-data_alert['ozone'][:-24],
              bins=40, range=(-40,40), histtype='step', color='cornflowerblue', ls='--', lw=2)
    ax41.hist(data_ref[0]['O3'][:,find_nearest(82.5, data_ref[0]['lat'], index=True),
                                find_nearest(360-62.3, data_ref[0]['lon'], index=True)]*1e9-data_alert['ozone'][:-24],
              bins=40, range=(-40,40), histtype='step', color='black', lw=2)
    ax41.set_ylabel("counts", y=-0.25)
    ax42 = plt.subplot(222)
    ax42.set_title("Barrow", ha='left', y=0.85, x=0.02)
    ax42.hist(data_mysic[0]['O3'][:,find_nearest(71.32, data_ref[0]['lat'], index=True),
                                  find_nearest(360-156.61, data_ref[0]['lon'], index=True)]*1e9-data_barrow['ozone'][:-24],
              bins=40, range=(-40,40))
    ax42.hist(data[0]['O3'][:,find_nearest(71.32, data_ref[0]['lat'], index=True),
                            find_nearest(360-156.61, data_ref[0]['lon'], index=True)]*1e9-data_barrow['ozone'][:-24],
              bins=40, range=(-40,40), histtype='step', color='cornflowerblue', ls='--', lw=2)
    ax42.hist(data_ref[0]['O3'][:,find_nearest(71.32, data_ref[0]['lat'], index=True),
                                find_nearest(360-156.61, data_ref[0]['lon'], index=True)]*1e9-data_barrow['ozone'][:-24],
              bins=40, range=(-40,40), histtype='step', color='black', lw=2)
    ax43 = plt.subplot(223)
    ax43.set_title("Zeppelin Mountain", ha='left', y=0.85, x=0.02)
    ax43.hist(data_mysic[0]['O3'][:,find_nearest(78.9, data_ref[0]['lat'], index=True),
                                  find_nearest(11.88, data_ref[0]['lon'], index=True)]*1e9-data_zeppelin['ozone'][:-24],
              bins=40, range=(-40,40))
    ax43.hist(data[0]['O3'][:,find_nearest(78.9, data_ref[0]['lat'], index=True),
                            find_nearest(11.88, data_ref[0]['lon'], index=True)]*1e9-data_zeppelin['ozone'][:-24],
              bins=40, range=(-40,40), histtype='step', color='cornflowerblue', ls='--', lw=2)
    ax43.hist(data_ref[0]['O3'][:,find_nearest(78.9, data_ref[0]['lat'], index=True),
                                find_nearest(11.88, data_ref[0]['lon'], index=True)]*1e9-data_zeppelin['ozone'][:-24],
              bins=40, range=(-40,40), histtype='step', color='black', lw=2)
    ax44 = plt.subplot(224)
    ax44.set_title("Summit", ha='left', y=0.85, x=0.02)
    ax44.hist(data_mysic[0]['O3'][24*30*7+24*16:,find_nearest(72.57, data_ref[0]['lat'], index=True),
                                  find_nearest(360-38.48, data_ref[0]['lon'], index=True)]*1e9-data_summit['ozone'][1776:-24],
              bins=40, range=(-40,40))
    ax44.hist(data[0]['O3'][24*30*7+24*16:,find_nearest(72.57, data_ref[0]['lat'], index=True),
                            find_nearest(360-38.48, data_ref[0]['lon'], index=True)]*1e9-data_summit['ozone'][1776:-24],
              bins=40, range=(-40,40), histtype='step', color='cornflowerblue', ls='--', lw=2)
    ax44.hist(data_ref[0]['O3'][24*30*7+24*16:,find_nearest(72.57, data_ref[0]['lat'], index=True),
                                find_nearest(360-38.48, data_ref[0]['lon'], index=True)]*1e9-data_summit['ozone'][1776:-24],
              bins=40, range=(-40,40), histtype='step', color='black', lw=2)
    ax44.set_xlabel("$O_{3}^{mysic}-O_{3}^{obs}$ (ppb)", x=-0.25)

# Show it
plt.show(block=False)
    

