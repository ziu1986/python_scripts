import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from scipy import fftpack
from read_sunspots import *
from mytools.met_tools import *
from mytools.netcdf_tools import *

# Clean up
plt.close('all')

# Load modules and EBAS reading routine
execfile('read_ebas.py')

# The data
nc_src = os.environ['DATA']+'/astra_data/observations/sunspots/SN_d_tot_V2.0.txt'
src_jergul = os.environ['DATA']+'/astra_data/observations/ozone/Jergul/NO0030R.*ozone*.nas'
src_karasjok = os.environ['DATA']+'/astra_data/observations/ozone/Karasjok/NO0055R.*ozone*.nas'

station_location = {"Jergul":(69.45,24.6),"Karasjok":(69.467,25.217),"Svanvik":(69.45,30.03)}

# Loop through EBAS data and transform them to pandas timeseries
try:
    data_jergul
except NameError:
    data_jergul = []
    data_karasjok = []
    for file in sorted(glob.glob(src_jergul)):
        print("Reading file %s" % (file))
        tmp = read_station_data(file)
        data_jergul.append(pd.Series(tmp['O3'],index=tmp['time']))   
    for file in sorted(glob.glob(src_karasjok)):
        print("Reading file %s" % (file))
        tmp = read_station_data(file)
        data_karasjok.append(pd.Series(tmp['O3'],index=tmp['time']))
    # Concatenate the lists
    data_jergul = pd.concat(data_jergul)
    data_karasjok = pd.concat(data_karasjok)
    data_jerg_kara = pd.concat((data_jergul, data_karasjok))
    # Sunspot data
    data_ss = read_sunspots(nc_src)

data_ss_jerg_kara = data_ss['1988-04-18':'2010-02-28']['Ntot']
data_ss_cumsum = (data_ss_jerg_kara.cumsum()/data_ss_jerg_kara.cumsum().max())

# Selecting two regions of low solar activity
solarlow_1 = data_ss_jerg_kara.where((data_ss_cumsum.resample('1M').mean().diff()<1.2e-3)&(data_ss_cumsum.resample('1M').mean().diff().diff().abs()<1e-4)).dropna()[:'1998'].index
solarlow_2 = data_ss_jerg_kara.where((data_ss_cumsum.resample('1M').mean().diff()<1.2e-3)&(data_ss_cumsum.resample('1M').mean().diff().diff().abs()<1e-4)).dropna()['1998':].index

solarhigh_1 = data_ss_jerg_kara.where((data_ss_cumsum.resample('1M').mean().diff()>7e-3)).dropna()[:'1993'].index
solarhigh_2 = data_ss_jerg_kara.where((data_ss_cumsum.resample('1M').mean().diff()>7e-3)).dropna()['1993':].index

sl_ozone_mean_1 = data_jerg_kara[solarlow_1[0]:solarlow_1[-1]].mean()
sl_ozone_mean_2 = data_jerg_kara[solarlow_2[0]:solarlow_2[-1]].mean()
sl_ozone_std_1 = data_jerg_kara[solarlow_1[0]:solarlow_1[-1]].std()
sl_ozone_std_2 = data_jerg_kara[solarlow_2[0]:solarlow_2[-1]].std()

delta_sl_ozone = sl_ozone_mean_2-sl_ozone_mean_1
std_sl_ozone = np.sqrt((sl_ozone_std_1/np.sqrt(len(data_jerg_kara[solarlow_1[0]:solarlow_1[-1]])))**2+(sl_ozone_std_2/np.sqrt(len(data_jerg_kara[solarlow_2[0]:solarlow_2[-1]])))**2)

sh_ozone_mean_1 = data_jerg_kara[solarhigh_1[0]:solarhigh_1[-1]].mean()
sh_ozone_mean_2 = data_jerg_kara[solarhigh_2[0]:solarhigh_2[-1]].mean()
sh_ozone_std_1 = data_jerg_kara[solarhigh_1[0]:solarhigh_1[-1]].std()
sh_ozone_std_2 = data_jerg_kara[solarhigh_2[0]:solarhigh_2[-1]].std()

delta_sh_ozone = sh_ozone_mean_2-sh_ozone_mean_1
std_sh_ozone = np.sqrt((sh_ozone_std_1/np.sqrt(len(data_jerg_kara[solarhigh_1[0]:solarhigh_1[-1]])))**2+(sh_ozone_std_2/np.sqrt(len(data_jerg_kara[solarhigh_2[0]:solarhigh_2[-1]])))**2)

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("timeseries_ozone_and_sunspots")
ax11 = plt.subplot()
ax11t = ax11.twinx()

data_jerg_kara.dropna().plot(ax=ax11, color='orange', ls='None', marker='+', label='Jergul/Karasjok')
#(data_ss_jerg_kara/data_ss['1988-04-18':'2010-02-28'].max()['Ntot']).plot(ax=ax11, marker='x', alpha=0.25, color='black', label='sunspots')

ax11t.plot(data_ss_jerg_kara.index, (data_ss_jerg_kara/data_ss_jerg_kara.max()), marker='x', alpha=0.25, color='black', label='sunspots')

#ax11t.axvspan(date2num(solarlow_1[0]),date2num(solarlow_1[-1]), alpha=0.10, color='blue')
#ax11t.axvspan(date2num(solarlow_2[0]),date2num(solarlow_2[-1]), alpha=0.10, color='blue')
#ax11t.axvspan(date2num(solarhigh_1[0]),date2num(solarhigh_1[-1]), alpha=0.10, color='red')
#ax11t.axvspan(date2num(solarhigh_2[0]),date2num(solarhigh_2[-1]), alpha=0.10, color='red')

ax11.set_ylabel("[$O_3$] (ppb)", color='orange')
ax11.tick_params(axis='y', colors='orange')
ax11.set_xlabel("Time (year)")
ax11.legend(loc='upper left')

ax11t.set_ylabel("$N^{ss}_{tot}$ (1/$N^{ss}_{max}$)")
ax11t.grid(False)
ax11t.legend(loc=(0.005,0.885))

ax11.axvspan(date2num(dt.datetime.strptime('1986-09','%Y-%m')),date2num(dt.datetime.strptime('1996-08','%Y-%m')),color='linen')
ax11.axvspan(date2num(dt.datetime.strptime('2008-12','%Y-%m')),date2num(dt.datetime.strptime('2011-01','%Y-%m')),color='linen')
ax11.text(date2num(dt.datetime.strptime('1989-11','%Y-%m')),(83),"solar cycle 22")
ax11.plot(date2num(dt.datetime.strptime('1989-11','%Y-%m')),(-5),ls='', marker='^', color='red')
ax11.plot(date2num(dt.datetime.strptime('1996-08','%Y-%m')),(83),ls='', marker='v', color='blue')
ax11.text(date2num(dt.datetime.strptime('2001-11','%Y-%m')),(83),"solar cycle 23")
ax11.plot(date2num(dt.datetime.strptime('2001-11','%Y-%m')),(-5),ls='', marker='^', color='red')
ax11.plot(date2num(dt.datetime.strptime('2008-12','%Y-%m')),(83),ls='', marker='v', color='blue')

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("hist-Nsunspots")
ax21 = plt.subplot(211)
ax21.hist(data_ss_jerg_kara.ravel(), bins=np.arange(-0.5,400,1))
#ax21.hist(data_ss_jerg_kara.ravel(), bins=range(400), histtype='step')
ax21.set_ylabel("counts")
ax21.set_xlabel("$N_{tot}^{ss}$")
ax21.set_xlim(0,40)
ax21.set_ylim(0,150)

ax22 = plt.subplot(212)
ax22t = ax22.twinx()
ax22tt = ax22.twinx()

data_ss_cumsum.plot(ax=ax22)
ax22.set_xlabel("Time (year)")
ax22.set_ylabel("$\Sigma N_{tot}^{ss}$", color='blue')
ax22.tick_params(axis='y', colors='blue')

data_ss_cumsum.resample('1M').mean().diff().plot(ax=ax22t, color='red', alpha=0.25, label='$\\nabla \Sigma N^{ss}$')

#ax22t.set_ylabel("$\\nabla$", color='red')
ax22t.tick_params(axis='y', colors='red')
ax22t.grid(False)

data_ss_cumsum.resample('1M').mean().diff().diff().plot(ax=ax22tt, color='black', alpha=0.25, label='$\\nabla^2 \Sigma N^{ss}$')
#ax22tt.set_ylabel("$\\nabla$", color='black')
ax22tt.tick_params(axis='y', colors='black')
ax22tt.grid(False)

ax22t.legend()
ax22tt.legend(loc=(0.91,0.825))

fig2.tight_layout()

fig3 = plt.figure(3,figsize=(16,9))
fig3.canvas.set_window_title("hist_ozone_solar_low_high")
ax31 = plt.subplot(221)
ax31.set_title("Solar low - Cycle 22/23")
ax31.hist(data_jerg_kara[solarlow_1[0]:solarlow_1[-1]], bins=range(0,80))
ax31.text(58,1200, "$\left<[O_3]\\right>$ = %1.2f ppb" % (sl_ozone_mean_1), size='large')
ax31.text(59.6,1100, "$\sigma_{[O_3]}$ = %1.2f ppb" % (sl_ozone_std_1), size='large')
ax31.text(59.65,1300, "$NDF$ = %d" % (len(data_jerg_kara[solarlow_1[0]:solarlow_1[-1]])), size='large')
ax32 = plt.subplot(222)
ax32.set_title("Solar low - Cycle 23/24")
ax32.hist(data_jerg_kara[solarlow_2[0]:solarlow_2[-1]], bins=range(0,80))
ax32.text(58,1200, "$\left<[O_3]\\right>$ = %1.2f ppb" % (sl_ozone_mean_2), size='large')
ax32.text(59.6,1100, "$\sigma_{[O_3]}$ = %1.2f ppb" % (sl_ozone_std_2), size='large')
ax32.text(59.65,1300, "$NDF$ = %d" % (len(data_jerg_kara[solarlow_2[0]:solarlow_2[-1]])), size='large')
ax33 = plt.subplot(223)
ax33.set_title("Solar high - Cycle 22")
ax33.hist(data_jerg_kara[solarhigh_1[0]:solarhigh_1[-1]], bins=range(0,80))
ax33.text(58,1200, "$\left<[O_3]\\right>$ = %1.2f ppb" % (sh_ozone_mean_1), size='large')
ax33.text(59.6,1100, "$\sigma_{[O_3]}$ = %1.2f ppb" % (sh_ozone_std_1), size='large')
ax33.text(59.65,1300, "$NDF$ = %d" % (len(data_jerg_kara[solarhigh_1[0]:solarhigh_1[-1]])), size='large')
ax34 = plt.subplot(224)
ax34.set_title("Solar high - Cycle 23")
ax34.hist(data_jerg_kara[solarhigh_2[0]:solarhigh_2[-1]], bins=range(0,80))
ax34.text(58,1200, "$\left<[O_3]\\right>$ = %1.2f ppb" % (sh_ozone_mean_2), size='large')
ax34.text(59.6,1100, "$\sigma_{[O_3]}$ = %1.2f ppb" % (sh_ozone_std_2), size='large')
ax34.text(59.65,1300, "$NDF$ = %d" % (len(data_jerg_kara[solarhigh_2[0]:solarhigh_2[-1]])), size='large')

ax33.set_xlabel("[$O_3$] (ppb)", x=1.1)
ax33.set_ylabel("count", y=1)
for ax in fig3.axes:
    ax.set_ylim(0,1400)
# Show it
plt.show(block=False)
