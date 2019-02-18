import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *
from scipy.optimize import curve_fit

func = lambda x, p: p[0]*x+p[1]

# Clean up
plt.close('all')

# Data sources
nc_src = os.environ['DATA']+'/nird_data/results/OsloCTM3/ozone25/ozone25_vmr/vmr_ozone*.nc'
nc_src_old = os.environ['DATA']+'/astra_data/ctm_results/CTM3_oivind/osloctm_ozone*.nc'
src_jergul = os.environ['DATA']+'/astra_data/observations/ozone/Jergul/NO0030R.*ozone*.nas'
src_karasjok = os.environ['DATA']+'/astra_data/observations/ozone/Karasjok/NO0055R.*ozone*.nas'
src_svanvik = os.environ['DATA']+'/astra_data/observations/ozone/Svanvik/NO0047R.*ozone*.nas'
src_svanvik_2018 = os.environ['DATA']+'/astra_data/observations/ozone/Svanvik/NO0047R.*ozone*.xls'

station_location = {"Jergul":(69.45,24.6),"Karasjok":(69.467,25.217),"Svanvik":(69.45,30.03)}

# Loop through EBAS data and transform them to pandas timeseries
try:
    data_jergul
except NameError:
    data_jergul = []
    data_karasjok = []
    data_svanvik = []
    for file in sorted(glob.glob(src_jergul)):
        print("Reading file %s" % (file))
        tmp = read_station_data_ebas(file)
        data_jergul.append(tmp)   
    for file in sorted(glob.glob(src_karasjok)):
        print("Reading file %s" % (file))
        tmp = read_station_data_ebas(file)
        data_karasjok.append(tmp)
    for file in sorted(glob.glob(src_svanvik)):
        print("Reading file %s" % (file))
        tmp = read_station_data_ebas(file)
        data_svanvik.append(tmp)
    data_svanvik_2018 = pd.read_excel(glob.glob(src_svanvik_2018)[0], index_col=0, header=0)
    data_svanvik_2018 = data_svanvik_2018['O3_mugm-3'].where(data_svanvik_2018['O3_mugm-3']>=0).dropna()/2.
    
    # Concatenate the lists
    data_jergul = pd.concat(data_jergul)
    data_karasjok = pd.concat(data_karasjok)
    data_svanvik = pd.concat(data_svanvik)
    # Round time index to full hour
    data_jergul.index = data_jergul.index.round("h")
    data_karasjok.index = data_karasjok.index.round("h")
    data_svanvik.index = data_svanvik.index.round("h")

# Read the model data
try:
    data
except NameError:
    data_list_finnmark = []
    data_list_svan = []
    for ifile in sorted(glob.glob(nc_src)):
        print("Reading %s" % (ifile))
        data = xr.open_dataset(ifile).isel(lev=0)
        data_list_finnmark.append(data['O3'].sel(lat=station_location['Karasjok'][0],method='nearest').sel(lon=station_location['Karasjok'][1],method='nearest')*1e9)
        data_list_svan.append(data['O3'].sel(lat=station_location['Svanvik'][0],method='nearest').sel(lon=station_location['Svanvik'][1],method='nearest')*1e9)
        
    
    data_old = read_data(nc_src_old,var='O3',lev=(1,1))
    # Shifting old data time stemp by 3 hours
    # since there seems to be something off with the daily cycle of ozone
    data_old['time'] = data_old.time.data+np.timedelta64(-3,'h')
# Save selections for multiple use
try:
    data_finnmark
except NameError:
    print("Concatenating...")
    data_finnmark = xr.concat(data_list_finnmark,dim='time')
    data_svan = xr.concat(data_list_svan,dim='time')    
    data_old_finnmark = (data_old.sel(lat=station_location['Karasjok'][0],method='nearest').sel(lon=station_location['Karasjok'][1],method='nearest')*1e9)
   
    data_old_svan = (data_old.sel(lat=station_location['Svanvik'][0],method='nearest').sel(lon=station_location['Svanvik'][1],method='nearest')*1e9)
# Selection of winter and summer data
data_finnmark_winter = data_finnmark.where((data_finnmark.time.dt.month<4) | (data_finnmark.time.dt.month>=10))
data_finnmark_summer = data_finnmark.where((data_finnmark.time.dt.month>=4) & (data_finnmark.time.dt.month<10))
data_old_finnmark_winter = data_old_finnmark.where((data_finnmark.time.dt.month<4) | (data_finnmark.time.dt.month>=10))
data_old_finnmark_summer = data_old_finnmark.where((data_finnmark.time.dt.month>=4) & (data_finnmark.time.dt.month<10))
data_svan_winter = data_svan.where((data_svan.time.dt.month<4) | (data_svan.time.dt.month>=10))
data_svan_summer = data_svan.where((data_svan.time.dt.month>=4) & (data_svan.time.dt.month<10))
data_old_svan_winter = data_old_svan.where((data_old_svan.time.dt.month<4) | (data_old_svan.time.dt.month>=10))
data_old_svan_summer = data_old_svan.where((data_old_svan.time.dt.month>=4) & (data_old_svan.time.dt.month<10))

data_jergul_winter = data_jergul.where((data_jergul.index.month<4) | (data_jergul.index.month>=10))
data_jergul_summer = data_jergul.where((data_jergul.index.month>=4) & (data_jergul.index.month<10))
data_karasjok_winter = data_karasjok.where((data_karasjok.index.month<4) | (data_karasjok.index.month>=10))
data_karasjok_summer = data_karasjok.where((data_karasjok.index.month>=4) & (data_karasjok.index.month<10))
data_svanvik_winter = data_svanvik.where((data_svanvik.index.month<4) | (data_svanvik.index.month>=10))
data_svanvik_summer = data_svanvik.where((data_svanvik.index.month>=4) & (data_svanvik.index.month<10))

# Is the year 2018 exceptional?
delta_svanvik = ((data_svanvik_2018).groupby(data_svanvik_2018.index.dayofyear).apply(np.nanmean)-data_svanvik.groupby(data_svanvik.index.dayofyear).apply(np.nanmean))

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_time_series")
data_finnmark.plot(ls="None", marker='o', fillstyle='none', color='black',label='OsloCTM3 v1.0')
data_old_finnmark.plot(ls="None", marker='d', fillstyle='none', color='blue',label="OsloCTM3 v0.1")
data_jergul.plot(marker='x', ls='none', color='red', alpha=0.15, label='Jergul')
data_karasjok.plot(marker='+', ls='none', color='orange', alpha=0.15, label='Karasjok')
data_svanvik.plot(marker='v', fillstyle='none', ls='none', color='blueviolet', alpha=0.15, label='Svanvik')

ax11 = plt.gca()
ax11.set_ylabel("[$O_3$] (ppb)")
ax11.set_xlabel("Time (year)")
ax11.legend(ncol=3)

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("hist-surface_ozone_concentration_distr")
ax21 = plt.subplot(221)
ax21.set_title("OsloCTM3: Karasjok")
ax22 = plt.subplot(222)
ax22.set_title("OsloCTM3: Svanvik")
ax23 = plt.subplot(223)
ax23.set_title("Observations")
ax24 = plt.subplot(224)
ax24.set_title("Observations")

bins = np.arange(0,80)
ax21.hist(data_finnmark, density=True, bins=bins, histtype='step', color='black', label='OsloCTM3 v1.0')
#ax21.hist(data_finnmark_winter, density=True, bins=bins, histtype='step', color='black', hatch='\\\\', label='OsloCTM3 v1.0 - winter')
#ax21.hist(data_finnmark_summer, density=True, bins=bins, histtype='step', color='black', hatch='//', label='OsloCTM3 v1.0 - summer')
ax21.hist(data_old_finnmark.data.flatten(), density=True, bins=bins, histtype='step', color='blue',label='OsloCTM3 v0.1')


ax22.hist(data_svan, density=True, bins=bins, histtype='step', color='black', label='OsloCTM3 v1.0')
ax22.hist(data_old_svan.data.flatten(), density=True, bins=bins, histtype='step', color='blue', label='OsloCTM3 v0.1')

ax23.hist(data_jergul.dropna(), bins=bins, density=True, histtype="step", color='red', label='Jergul')
ax23.hist(data_karasjok.dropna(), bins=bins, density=True, histtype="step", color='orange', label='Karasjok')

ax24.hist(data_svanvik.dropna(), bins=bins, density=True, histtype="step", color='blueviolet', label='Svanvik')
ax24.hist(data_svanvik_2018.dropna(), bins=bins, density=True, histtype="step", color='fuchsia', label='Svanvik 2018')

for ax in fig2.axes:
    ax.set_ylabel("$\\rho$")
    ax.set_ylim(0,0.06)
    ax.legend()
ax23.set_xlabel("[$O_3$] (ppb)", x=1.1)

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title('ozone_daily_cycle')
ax31 = plt.subplot(211)
# OsloCTM3 v1.0
data_finnmark.groupby('time.hour').mean().plot(ax=ax31, ls="None", marker='o', fillstyle='none',color='black', label='OsloCTM3 v1.0')
# Split into summer and winter
data_finnmark_winter.groupby('time.hour').mean().plot(ax=ax31, ls="--", color='black', alpha=0.5, label='OsloCTM3 v1.0 - winter')
data_finnmark_summer.groupby('time.hour').mean().plot(ax=ax31, ls="-.", color='black', alpha=0.5, label='OsloCTM3 v1.0 - summer')
# OsloCTM3 v0.1
data_old_finnmark.groupby('time.hour').mean().plot(ax=ax31, ls="None", marker='d', fillstyle='none',color='blue', label='OsloCTM3 v0.1')
# Split into summer and winter
data_old_finnmark_winter.groupby('time.hour').mean().plot(ax=ax31, ls="--", color='blue', alpha=0.5, label='OsloCTM3 v0.1 - winter')
data_old_finnmark_summer.groupby('time.hour').mean().plot(ax=ax31, ls="-.", color='blue', alpha=0.5, label='OsloCTM3 v0.1 - summer')
# Observations Jergul
data_jergul.groupby(data_jergul.index.hour).apply(np.nanmean).plot(marker='x', ls='none', color='red', label='Jergul')
# Split into summer and winter
data_jergul_winter.groupby(data_jergul.index.hour).apply(np.nanmean).plot(ax=ax31, ls="--", color='red', alpha=0.5, label='Jergul - winter')
data_jergul_summer.groupby(data_jergul.index.hour).apply(np.nanmean).plot(ax=ax31, ls="-.", color='red', alpha=0.5, label='Jergul - summer')
# Observations Karasjok
data_karasjok.groupby(data_karasjok.index.hour).apply(np.nanmean).plot(marker='+', ls='none', color='orange', label='Karasjok')
# Split into summer and winter
data_karasjok_winter.groupby(data_karasjok.index.hour).apply(np.nanmean).plot(ax=ax31, ls="--", color='orange', alpha=0.5, label='Karasjok - winter')
data_karasjok_summer.groupby(data_karasjok.index.hour).apply(np.nanmean).plot(ax=ax31, ls="-.", color='orange', alpha=0.5, label='Karasjok - summer')

#ax31.set_xlabel("Time (UTC)")
#ax31.set_ylabel("[$O_3$] (ppb)")
ax31.set_ylim(20,50)
ax31.legend(ncol=4)

ax32 = plt.subplot(212)
# OsloCTM3 v1.0
data_svan.groupby('time.hour').mean().plot(ax=ax32, ls="None", marker='o', fillstyle='none',color='black', label='OsloCTM3 v1.0')
# Split into summer and winter
data_svan_winter.groupby('time.hour').mean().plot(ax=ax32, ls="--", color='black', alpha=0.5, label='OsloCTM3 v1.0 - winter')
data_svan_summer.groupby('time.hour').mean().plot(ax=ax32, ls="-.", color='black', alpha=0.5, label='OsloCTM3 v1.0 - summer')
# OsloCTM3 v0.1
data_old_svan.groupby('time.hour').mean().plot(ax=ax32, ls="None", marker='d', fillstyle='none',color='blue', label='OsloCTM3 v0.1')
# Split into summer and winter
data_old_svan_winter.groupby('time.hour').mean().plot(ax=ax32, ls="--", color='blue', alpha=0.5, label='OsloCTM3 v0.1 - winter')
data_old_svan_summer.groupby('time.hour').mean().plot(ax=ax32, ls="-.", color='blue', alpha=0.5, label='OsloCTM3 v0.1 - summer')
# Observations Svanvik
data_svanvik.groupby(data_svanvik.index.hour).apply(np.nanmean).plot(marker='v', fillstyle='none', ls='none', color='blueviolet', label='Svanvik')
# Split into summer and winter
data_svanvik_winter.groupby(data_svanvik.index.hour).apply(np.nanmean).plot(ax=ax32, ls="--", color='blueviolet', alpha=0.5, label='Svanvik - winter')
data_svanvik_summer.groupby(data_svanvik.index.hour).apply(np.nanmean).plot(ax=ax32, ls="-.", color='blueviolet', alpha=0.5, label='Svanvik - summer')

ax32.set_xlabel("Time (UTC)")
ax32.set_ylabel("[$O_3$] (ppb)", y=1)
ax32.set_ylim(20,50)
ax32.legend(ncol=3)

fig4 = plt.figure(4,figsize=(16,9))
fig4.canvas.set_window_title("svanvik-ozone_climatologyVS2018")

ax41 = plt.subplot(121)
(data_svanvik_2018).groupby(data_svanvik_2018.index.dayofyear).apply(np.nanmean).plot(ls='None', marker='^', fillstyle='none', color='fuchsia', label="Svanvik 2018")
data_svanvik.groupby(data_svanvik.index.dayofyear).apply(np.nanmean).plot(ls='None',marker='v', fillstyle='none', color='blueviolet', label='Svanvik 1986-1996')
(data_svanvik.groupby(data_svanvik.index.dayofyear).apply(np.nanmean)+data_svanvik.groupby(data_svanvik.index.dayofyear).apply(np.std)).plot(color='blueviolet', alpha=0.25, label='_')
(data_svanvik.groupby(data_svanvik.index.dayofyear).apply(np.nanmean)-data_svanvik.groupby(data_svanvik.index.dayofyear).apply(np.std)).plot(color='blueviolet', alpha=0.25, label='_')
(data_svan.groupby('time.dayofyear').mean()).plot(marker='o', fillstyle='none', ls='None', color='black', label='OsloCTM3 v1.0')
(data_old_svan.groupby('time.dayofyear').mean()).plot(marker='d', fillstyle='none', ls='None', color='blue', label='OsloCTM3 v0.1')
#(data_svanvik).groupby(data_svanvik.index.dayofyear).apply(np.nanmean).rolling(15).mean().plot()
ax41.set_ylabel("[$O_3$] (ppb)")
ax41.set_xlabel("Time (day of year)")
ax41.legend()

ax42 = plt.subplot(122)
delta_svanvik.hist(bins=np.arange(-20,21,2), orientation='horizontal', color='fuchsia')
ax42.axhline(delta_svanvik.mean(), color='grey', ls='--')
ax42.text(30,0, "<%1.2f> $\pm$ %1.2f" % (delta_svanvik.mean(),delta_svanvik.std()),size='large')
ax42.set_ylabel("$\Delta [O_3]$ (ppb)")
ax42.set_xlabel("count")

#ax43 = plt.subplot(133)
#ax43.hist(data_finnmark.where(data_finnmark>50).groupby('time.dayofyear').count(),bins=np.arange(1,366))
#ax43.set_xlabel("Time (day of year)")
#ax43.set_ylabel("count")

fig5 = plt.figure(5, figsize=(16,9))
fig5.canvas.set_window_title("surface_ozone_obsvsmodel_svanvik1993")
ax52 = plt.subplot(121)
hist_svanvik = ax52.hist2d(data_svan.sel(time=data_svanvik['1993'][::3].dropna().index).data, data_svanvik['1993'][::3].dropna().data.flatten(), bins=range(71), cmap=plt.cm.hot_r)
ax52.plot(np.arange(0,71),np.arange(0,71), color='grey', ls=':')
#ax52.plot(np.arange(0,71),np.arange(0,71)-simple_bias_corr.data, color='blue', ls=':')
cb = fig5.colorbar(hist_svanvik[3], ax=ax52)
cb.set_label("counts")
ax52.set_ylabel("$[O_3]_{obs}$ (ppb)")
ax52.set_xlabel("$[O_3]_{model} (ppb)$")

# Fitting the histrogram
#x = np.arange(0.5,70)
#y = np.arange(0.5,70)
#X, Y = np.meshgrid(x,y)
#weights = np.sqrt(hist_svanvik[0].T)
#popt2, pcov2 = np.polyfit(X.ravel(), Y.ravel(), 1, w=weights.ravel(), cov=True)
#ax52.plot(np.arange(70), func(np.arange(70),popt2), color='black', ls='--')

ax51 = plt.subplot(122)
#(func(data_svan.where(data_svan.time.dt.year==1993,drop=True),popt2)).plot(ls='None', color='black', marker='o', fillstyle='none', label='OsloCTM3 v1.0 - corr')
data_svanvik['1993'].plot(ls='None',marker='v', fillstyle='none', color='blueviolet', label='Svanvik 1986-1996')
ax51.set_xlabel("Time (month)")
ax51.set_ylabel("[$O_3$] (ppb)")
ax51.legend()

print("Correlation: %s" % np.corrcoef(data_svanvik['1993'][::3].dropna().data.flatten(), data_svan.sel(time=data_svanvik['1993'][::3].dropna().index).data))

from scipy import fftpack
fig11 = plt.figure(11,figsize=(16,9))
fig11.canvas.set_window_title("ozone_freq_spectrum")
data_jerg_kara = pd.concat((data_jergul, data_karasjok))

ax111 = plt.subplot(211)
fft_karasjok = fftpack.fft(data_jerg_kara.resample('1M').mean().fillna(method='ffill'))
freqs_karasjok = fftpack.fftfreq(len(fft_karasjok))
ax111.stem(1/freqs_karasjok/12, np.abs(fft_karasjok))
ax111.set_title("Jergul $\mathrm{\cap}$ Karasjok",x=0.9,y=0.85)
ax111.set_xlim(0,25)

ax112 = plt.subplot(212)
fft_finnmark = fftpack.fft(data_finnmark.resample(time='1M').mean())
freqs_finnmark = fftpack.fftfreq(len(fft_finnmark))
ax112.stem(1/freqs_finnmark/12, np.abs(fft_finnmark))
ax112.set_title("OsloCTM3 v1.0",x=0.9,y=0.85)
ax112.set_ylabel("Amplitude")
ax112.set_xlabel("Frequency (years)")
ax112.set_xlim(0,25)

fig12 = plt.figure(12, figsize=(16,9))
fig12.canvas.set_window_title("surface_ozone_jerg_kara_filtered")

rm_data_jerg_kara = data_jerg_kara.fillna(method='ffill').rolling(window=1*365*24,center=True).mean()
rm_seasonal_data_jerg_kara = data_jerg_kara.fillna(method='ffill').rolling(window=1*365*12,center=True).mean()
epoch_index = (rm_data_jerg_kara.dropna().index-pd.Timestamp('1988-01-01')) // pd.Timedelta('1h')

popt_obs_rm, pcov_obs_rm = np.polyfit(epoch_index, rm_data_jerg_kara.dropna(), 1, cov=True)
print("GradientO_3 = %1.3f ppb/a (%1.3f) " % (popt_obs_rm[0]*(epoch_index[-1]-epoch_index[0])/((rm_data_jerg_kara.index[-1]-rm_data_jerg_kara.index[0]) / np.timedelta64(1, 'Y')), popt_obs_rm[0]*(epoch_index[-1]-epoch_index[0])/((rm_data_jerg_kara.index[-1]-rm_data_jerg_kara.index[0]) / np.timedelta64(1, 'Y'))/popt_obs_rm[1]*100))

ax121 = plt.subplot(221)
data_jerg_kara.plot(ax=ax121, color='orange', ls='None', marker='+', label='data')
#data_jerg_kara.resample('1M').mean().fillna(method='ffill').plot(color='red', label="monthly ave, nan filled")
rm_seasonal_data_jerg_kara.plot(color='red', label='seasonal filter')
rm_data_jerg_kara.plot(color='black', label='low-pass filter')
#data_jerg_kara.interpolate(method='spline',order=3, s=3).plot(color='black', ls='--')
ax121.set_ylabel("$[O_3]$ (ppb)")
ax121.set_xlabel("Time (year)")

ax121.legend(ncol=3)

ax123 = plt.subplot(223, sharex=ax121)
rm_seasonal_data_jerg_kara = (data_jerg_kara-rm_data_jerg_kara).fillna(method='ffill').rolling(window=1*365*12,center=True).mean()

ax123.set_title("Data series Jergul/Karasjok short term variations")
(data_jerg_kara-rm_data_jerg_kara-rm_seasonal_data_jerg_kara)['1989':'2010'].plot(ax=ax123)
ax123.set_xlabel("Time (year)")
ax123.set_ylabel("$\Delta[O_3]$ (ppb)")

ax122 = plt.subplot(222, sharex=ax121)
rm_seasonal_data_finnmark = data_finnmark.rolling(time=1*365*4, center=True).mean()
epoch_model = np.arange(0,len(rm_seasonal_data_finnmark.dropna(dim='time'))*3, 3)
popt_model_rm, pcov_model_rm = np.polyfit(epoch_model, rm_seasonal_data_finnmark.dropna(dim='time'), 1, cov=True)
print("GradientO_3 = %1.3f ppb/a" % (popt_model_rm[0]*np.timedelta64((rm_seasonal_data_finnmark.dropna(dim='time').time[-1]-rm_seasonal_data_finnmark.dropna(dim='time').time[0]).data, 'h').astype(float)/23.5))

data_finnmark.plot(color='black', marker='o', fillstyle='none', ls='None', label='OsloCTM3 v1.0')
rm_seasonal_data_finnmark.plot(ax=ax122, color='red')
ax122.set_xticks("")

ax124 = plt.subplot(224, sharex=ax121, sharey=ax123)
(data_finnmark-rm_seasonal_data_finnmark).plot(ax=ax124)
ax124.set_xlabel("Time (year)")

fig13 = plt.figure(13, figsize=(16,9))
fig13.canvas.set_window_title("surface_ozone_correlation_obs_model")
ax131 = plt.subplot(211)
data_sel_jerg_kara = data_jerg_kara[data_finnmark.sel(time=slice('1991-01-01T00','2010-02-28T23')).time.data].dropna()
data_sel_finnmark = data_finnmark.sel(time=data_sel_jerg_kara.index)
hist_sel = ax131.hist2d(data_sel_finnmark, data_sel_jerg_kara, bins=range(0,70), cmap=plt.cm.hot_r)

ax131.plot(np.arange(0,70), color='grey', ls=':')
cb = fig13.colorbar(hist_sel[3], ax=ax131)
cb.set_label("counts")

print("Correlation: %s" % np.corrcoef(data_sel_finnmark, data_sel_jerg_kara))

ax132 = plt.subplot(212)
ax132.set_title("Corrected", x=0.9, y=0.9)
rm_data_sel_jerg_kara = (data_jerg_kara-rm_data_jerg_kara)[(data_finnmark).sel(time=slice('1991-04-02T06','2009-08-31T12')).time.data].dropna()
rm_data_sel_finnmark = (data_finnmark).sel(time=rm_data_sel_jerg_kara.index)
try:
    hist_jerg_kara
except NameError:
    hist_jerg_kara = ax132.hist2d(rm_data_sel_finnmark, rm_data_sel_jerg_kara, bins=range(0,71), cmap=plt.cm.hot_r)
ax132.plot(np.arange(0,71),np.arange(0,71), color='grey', ls=':')
cb = fig13.colorbar(hist_jerg_kara[3], ax=ax132)
cb.set_label("counts")
ax132.set_ylabel("$[O_3]_{obs}$ (ppb)", y=1)
ax132.set_xlabel("$[O_3]_{model} (ppb)$")

print("Correlation: %s" % np.corrcoef(rm_data_sel_finnmark, rm_data_sel_jerg_kara))

# Fitting the histrogram
#x = np.arange(-19.5, 20)
#y = np.arange(-19.5, 20)
#X, Y = np.meshgrid(x,y)
#weights = np.sqrt(hist_jerg_kara[0].T).ravel()
#popt_jerg_kara, pcov_jerg_kara = np.polyfit(X.ravel(), Y.ravel(), 1, w=weights, cov=True)
#ax132.plot(np.arange(-20,21), func(np.arange(-20,21),popt_jerg_kara), ls='--', color='black')
# Show it
plt.show(block=False)
