import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from scipy import fftpack
from read_sunspots import *
from mytools.met_tools import *
from mytools.netcdf_tools import *
#from mytools.cartopy_tools import scale_bar
from station_info import station_location


def load_data(src):
    data = []
    for file in sorted(glob.glob(src)):
        print("Reading file %s" % (file))
        tmp = read_station_data_ebas(file)
        data.append(tmp['O3'])   
    # Concatenate the lists
    print('Concatenating data...')
    data = pd.concat(data)
    # Round to full hours
    data.index = data.index.round("h")
    return(data)

# Clean up
plt.close('all')

# Load functions to read NOAA data
#execfile("../BrXplo/read_station_data.py")
#from read_station_data import read_station_data_noaa
src = os.environ['DATA']+'/astra_data/observations/ozone/'
src_svanvik_2018 = os.environ['DATA']+'/astra_data/observations/ozone/Svanvik/NO0047R.*ozone*.xls'
src_stations = ('Barrow', 'Esrange', 'Janiskoski', 'Jergul', 'Karasjok', 'Pallas', 'Prestebakke', 'Svanvik')

try:
    data_barrow
except NameError:
    data_barrow = []
    for file in sorted(glob.glob(src+src_stations[0]+'/*')):
        if int(file[-4:]) < 2003:
            print(file)
            tmp = (read_station_data_noaa(file, utc=-9, start_data=28))
        elif int(file[-4:]) < 2012:
            print(file)
            tmp = (read_station_data_noaa(file, utc=-9, station='other', column=5))
        data_barrow.append(tmp)
    # Concatenate the lists
    data_barrow = pd.concat(data_barrow)
    # Round time index to full hour
    data_barrow.index = data_barrow.index.round("h")

try:
    data
except NameError:
    data = {}
    for station in src_stations[1:]:
            data.update({station:load_data(src+station+'/*.nas')})
        
data_prestebakke = data['Prestebakke']
data_jergkara = pd.concat((data['Jergul'], data['Karasjok']))
data_svanvik_2018 = pd.read_excel(glob.glob(src_svanvik_2018)[0], index_col=0, header=0)
data_svanvik_2018 = data_svanvik_2018['O3_mugm-3'].where(data_svanvik_2018['O3_mugm-3']>=0).dropna()/2.    
# Spectral analysis
from scipy import fftpack
fft_barrow = fftpack.fft(data_barrow.resample('1M').mean().fillna(method='ffill'))
freqs_barrow = fftpack.fftfreq(len(fft_barrow))

fft_prestebakke = fftpack.fft(data_prestebakke.resample('1M').mean().fillna(method='ffill'))
freqs_prestebakke = fftpack.fftfreq(len(fft_prestebakke))

fft_jergkara = fftpack.fft(data_jergkara.resample('1M').mean().fillna(method='ffill'))
freqs_jergkara = fftpack.fftfreq(len(fft_jergkara))

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_timeseries_ltobs")
ax11 = plt.subplot(311)
ax12 = plt.subplot(312, sharex=ax11)
ax13 = plt.subplot(313, sharex=ax11)
data_barrow.plot(ax=ax11, ls='None', marker='x', label='Utqiagvik (USA)')
data_prestebakke.plot(ax=ax12, ls='None', marker='.', label='Prestebakke (NOR)', color='red')
data_jergkara.plot(ax=ax13, ls='None', marker='+', label='Jergul/Karasjok (NOR)', color='orange')

ax13.set_xlabel("Time (year)")
ax12.set_ylabel("[$O_3$] (ppb)")

for ax in fig1.axes:
    ax.set_ylim(0,100)
    ax.legend()

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("ozone_timeseries_fenoscandicobs")
ax21 = plt.subplot(511)
ax22 = plt.subplot(512, sharex=ax21)
ax23 = plt.subplot(513, sharex=ax21)
ax24 = plt.subplot(514, sharex=ax21)
ax25 = plt.subplot(515, sharex=ax21)
ax21.plot(data['Esrange'].index, data['Esrange'], ls='None', marker='+', label='Esrange (SWE)', color='blue')
ax22.plot(data['Pallas'].index, data['Pallas'], ls='None', marker='+', label='Pallas (FIN)', color='black')
data_jergkara.plot(ax=ax23, ls='None', marker='+', label='Jergul/Karasjok (NOR)', color='orange')
ax24.plot(data['Svanvik'].index, data['Svanvik'], ls='None', marker='+', label='Svanvik (NOR)', color='blueviolet')
#ax24.plot(data_svanvik_2018.index, data_svanvik_2018, ls='None', marker='x', label='Svanvik, 2018 (NOR)', color='blueviolet')
ax25.plot(data['Janiskoski'].index, data['Janiskoski'], ls='None', marker='+', label='Janiskoski (RUS)', color='grey')

ax25.set_xlabel("Time (year)")
ax23.set_ylabel("[$O_3$] (ppb)", y=1)
for ax in fig2.axes:
    ax.set_ylim(0,100)
    ax.legend()

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("ozone_climatology_fenoscandicobs")
ax31 = plt.subplot()

data['Esrange'].groupby(data['Esrange'].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Esrange (SWE)', color='blue')
data['Pallas'].groupby(data['Pallas'].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Pallas (FIN)', color='black')
data_jergkara.groupby(data_jergkara.index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Jergul/Karasjok (NOR)', color='orange')
#data_prestebakke.groupby(data_prestebakke.index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Prestebakke (NOR)', color='red')
data['Svanvik'].groupby(data['Svanvik'].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Svanvik (NOR)', color='blueviolet')
data['Janiskoski'].groupby(data['Janiskoski'].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Janiskoski (RUS)', color='grey')

ax31.set_xlabel("Time (day of year)")
ax31.set_ylabel("[$O_3$] (ppb)")
for ax in fig3.axes:
    ax.set_ylim(0,100)
    ax.legend()

if(False):
    fig2.canvas.set_window_title("fequency_spectrum")
    ax21 = plt.subplot()
    ax21.stem(1/freqs_barrow/12, np.abs(fft_barrow)/np.abs(fft_barrow).max(), label='Utqiagvik (USA)')
    markerline, stemlines, baseline = ax21.stem(1/freqs_prestebakke/12, np.abs(fft_prestebakke)/np.abs(fft_prestebakke).max(), label='Prestebakke (NOR)')
    plt.setp(stemlines, color='red')
    plt.setp(markerline, color='red')
    markerline, stemlines, baseline = ax21.stem(1/freqs_jergkara/12, np.abs(fft_jergkara)/np.abs(fft_jergkara).max(), label='Jergul/Karasjok (NOR)')
    plt.setp(stemlines, color='orange')
    plt.setp(markerline, color='orange')
    ax21.set_xlim(0,40)
    ax21.set_ylabel("Normalized Amplitude")
    ax21.set_xlabel("Frequency (years)")

    ax21.legend()

# Show the stations
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
#import cartopy.geodesic as cgeo
from matplotlib.transforms import offset_copy

fig4 = plt.figure(4)
fig4.canvas.set_window_title("station_map_fennoscandia")
stamen_terrain = cimgt.Stamen('terrain-background')
ax41 = fig4.add_subplot(1,1,1, projection=stamen_terrain.crs) #ccrs.PlateCarree()
ax41.set_extent([19.4,31.4,67.6,71.4], crs=ccrs.PlateCarree())
ax41.add_image(stamen_terrain, 8)

ax41.plot(station_location['Jergul'].lon, station_location['Jergul'].lat, fillstyle='left', marker='o', color='orange', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
ax41.plot(station_location['Karasjok'].lon, station_location['Karasjok'].lat, fillstyle='right', marker='o', color='orange', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
ax41.plot(station_location['Svanvik'].lon, station_location['Svanvik'].lat, marker='o', color='blueviolet', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
ax41.plot(station_location['Esrange'].lon, station_location['Esrange'].lat, marker='o', color='blue', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
ax41.plot(station_location['Pallas'].lon, station_location['Pallas'].lat, marker='o', color='black', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
ax41.plot(station_location['Janiskoski'].lon, station_location['Janiskoski'].lat, marker='o', color='grey', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax41)
text_transform = offset_copy(geodetic_transform, units='dots', x=-25)
text_transform_2 = offset_copy(geodetic_transform, units='dots', y=25)

# Add text.
ax41.text(station_location['Jergul'].lon, station_location['Jergul'].lat, u'Jergul',
          verticalalignment='center', horizontalalignment='right',
          transform=text_transform,
          bbox=dict(facecolor='orange', alpha=0.5, boxstyle='round'))
ax41.text(station_location['Karasjok'].lon, station_location['Karasjok'].lat, u'Karasjok',
          verticalalignment='center', horizontalalignment='right',
          transform=text_transform_2,
          bbox=dict(facecolor='orange', alpha=0.5, boxstyle='round'))
ax41.text(station_location['Svanvik'].lon, station_location['Svanvik'].lat, u'Svanvik',
          verticalalignment='center', horizontalalignment='right',
          transform=text_transform,
          bbox=dict(facecolor='blueviolet', alpha=0.5, boxstyle='round'))
ax41.text(station_location['Esrange'].lon, station_location['Esrange'].lat, u'Esrange',
          verticalalignment='center', horizontalalignment='right',
          transform=text_transform,
          bbox=dict(facecolor='blue', alpha=0.5, boxstyle='round'))
ax41.text(station_location['Pallas'].lon, station_location['Pallas'].lat, u'Pallas',
          verticalalignment='center', horizontalalignment='right',
          transform=text_transform,
          bbox=dict(facecolor='black', alpha=0.5, boxstyle='round'))
ax41.text(station_location['Janiskoski'].lon, station_location['Janiskoski'].lat, u'Janiskoski',
          verticalalignment='center', horizontalalignment='right',
          transform=text_transform,
          bbox=dict(facecolor='grey', alpha=0.5, boxstyle='round'))


fig5 = plt.figure(5, figsize=(16,9))
fig5.canvas.set_window_title("density_distribution")
ax51 = plt.subplot(221)
ax52 = plt.subplot(222)
ax53 = plt.subplot(223)
ax54 = plt.subplot(224)
bins = range(81)
hist2d_1 = ax51.hist2d(data['Esrange'].dropna(), data_jergkara[data['Esrange'].dropna().index], bins=(bins), cmap=plt.cm.hot_r)
hist2d_2 = ax52.hist2d(data['Pallas'].dropna(), data_jergkara[data['Pallas'].dropna().index], bins=(bins), cmap=plt.cm.hot_r)
hist2d_3 = ax53.hist2d(data['Esrange'].dropna(), data['Svanvik'][data['Esrange'].dropna().index], bins=(bins), cmap=plt.cm.hot_r)
hist2d_4 = ax54.hist2d(data['Pallas'].dropna(), data['Svanvik'][data['Pallas'].dropna().index], bins=(bins), cmap=plt.cm.hot_r)

for ax in fig5.axes:
    ax.plot(bins, bins, color='grey',ls='--')
    
# Placeing the colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
for ax,hist in zip(fig5.axes, (hist2d_1,hist2d_2,hist2d_3,hist2d_4)):
    axins = inset_axes(ax,
                       width="2.5%",  # width = 50% of parent_bbox width
                       height="70%",  # height : 5%
                       loc='lower right',
                       #bbox_to_anchor=(1, 0., 1, 1),
                       #bbox_transform=ax.transAxes,
                       #borderpad=0,
    )
    fig5.colorbar(hist[-1], cax=axins, orientation='vertical')
    axins.set_ylabel("Counts")

ax51.set_ylabel("$[O_3]_{Jergul/Karasjok}$ (ppb)")
ax53.set_ylabel("$[O_3]_{Svanvik}$ (ppb)")
ax53.set_xlabel("$[O_3]_{Esrange}$ (ppb)")
ax54.set_xlabel("$[O_3]_{Pallas}$ (ppb)")

ax51.text(60,75, "$r^2 = %1.2f$" % (data['Esrange'].corr(data_jergkara)), size='large')
ax52.text(60,75, "$r^2 = %1.2f$" % (data['Pallas'].corr(data_jergkara)), size='large')
ax53.text(60,75, "$r^2 = %1.2f$" % (data['Esrange'].corr(data['Svanvik'])), size='large')
ax54.text(60,75, "$r^2 = %1.2f$" % (data['Pallas'].corr(data['Svanvik'])), size='large')

fig6 = plt.figure(6, figsize=(16,9))
fig6.canvas.set_window_title("time_lag_correlation")
time_lag = range(-17,18)
lag_1 = []
lag_2 = []
lag_3 = []
lag_4 = []
lag_5 = []
lag_6 = []
lag_7 = []
for i in time_lag:
    #print("%d %1.2f" % (i, time_lagged_corr(data_jergkara, data['Esrange'], lag=i, pandas=True)))
    lag_1.append(time_lagged_corr(data_jergkara, data['Esrange'], lag=i, pandas=True))
    lag_2.append(time_lagged_corr(data_jergkara, data['Pallas'], lag=i, pandas=True))
    lag_3.append(time_lagged_corr(data['Svanvik'], data['Esrange'], lag=i, pandas=True))
    lag_4.append(time_lagged_corr(data['Svanvik'], data['Pallas'], lag=i, pandas=True))
    lag_5.append(time_lagged_corr(data['Svanvik'], data_jergkara, lag=i, pandas=True))
    lag_6.append(time_lagged_corr(data['Svanvik'], data['Janiskoski'], lag=i, pandas=True))
    lag_7.append(time_lagged_corr(data_jergkara, data['Janiskoski'], lag=i, pandas=True))
# Print maximum in lag
for i,lag in zip(range(1,8),(lag_1, lag_2, lag_3, lag_4, lag_5, lag_6, lag_7)):
    print("lag_%d max at %d h" % (i, np.array(time_lag)[np.where(np.array(lag)==np.array(lag).max())[0]]))
    
ax61 = plt.subplot(121)
ax62 = plt.subplot(122)
#ax63 = plt.subplot(133)
ax61.set_title("Jergul/Karasjok")
ax62.set_title("Svanvik")
#ax63.set_title("Svanvik")
ax61.plot(time_lag, lag_1, color='blue', label='Esrange')
ax61.plot(time_lag, lag_2, color='black', label='Pallas')
ax61.plot(np.array(time_lag)*(-1), lag_5, ls='--', color='blueviolet', label='Svanvik')
ax61.plot(time_lag, lag_7, color='grey', ls='--', label='Janiskoski')

ax62.plot(time_lag, lag_3, color='blue', ls='--', label='Esrange')
ax62.plot(time_lag, lag_4, color='black', ls='--', label='Pallas')
ax62.plot(time_lag, lag_5, color='orange', ls='--', label='Jergul/Karasjok')
ax62.plot(time_lag, lag_6, color='grey', label='Janiskoski')

for ax in fig6.axes:
    ax.set_xlabel('Lag (hours)')
    #ax.set_ylim(0,1)
    ax.legend()
ax61.set_ylabel('Correlation Coefficient')
# Show it
plt.show(block=False)



