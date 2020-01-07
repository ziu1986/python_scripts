import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from matplotlib.dates import date2num # Convert dates to matplotlib axis coords
from matplotlib import dates
from scipy import fftpack
from scipy import stats
#from read_sunspots import *
from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.ozone_tools import *
#from mytools.cartopy_tools import scale_bar
from mytools.station_info import station_location

#---------------------------------------------------------------------------------------------------------------------------------
# Read data
src = os.environ['DATA']+'/astra_data/observations/ozone/'
src_svanvik_OzoNorClim = os.environ['DATA']+'/astra_data/observations/ozone/Svanvik/NO0047R.*ozone*.xls'
src_stations = ('Barrow', 'Esrange', 'Janiskoski', 'Jergul', 'Karasjok', 'Pallas', 'Prestebakke', 'Svanvik')

try:
    data
except NameError:
    data = {}
    for station in src_stations:
        if station=='Barrow':
            data.update({station:load_data(src+station+'/*', type="Barrow")})
        else:
            data.update({station:load_data(src+station+'/*.nas')})

    # Concate Jergul and Karasjok data
    data_jergkara = pd.concat((data['Jergul'], data['Karasjok']))

    # Read and convert xls file data
    data_svanvik_OzoNorClim = []
    for file in sorted(glob.glob(src_svanvik_OzoNorClim)):
        tmp_data_svanvik = pd.read_excel(file, index_col=0, header=0)
        data_svanvik_OzoNorClim.append(tmp_data_svanvik['O3_mugm-3'].where(tmp_data_svanvik['O3_mugm-3']>=0).dropna()/2.)
    # Concat data Svanvik data
    data_svanvik_OzoNorClim = pd.concat(data_svanvik_OzoNorClim)

#---------------------------------------------------------------------------------------------------------------------------------
#Control
plot_timeseries = False
plot_timelag = False
plot_correlation = False
plot_splines = True
plot_climatology = False
plot_spectrum = False
plot_map = False
plot_aot = False
plot_rollingsum = True

#---------------------------------------------------------------------------------------------------------------------------------


if plot_timelag:
    execfile("ozone_observation_timelag.py")

execfile("ozone_observation_histograms.py")

execfile("ozone_observation_climatology.py")

# Resample
svanvik_daily = data['Svanvik'].resample('1d').apply(np.nanmean)
svanvik_daily_2018 = data_svanvik_OzoNorClim['2018'].resample('1d').apply(np.nanmean)


# Draw sample from climatology of Jergul/Karsjok, Esrange, Pallas -> fig9
sample = fitSpl_dmean(svanvik_daily.dropna().index.dayofyear)
sample_2018 = fitSpl_dmean(svanvik_daily_2018.dropna().index.dayofyear)
# Draw sample from Svanvik climatology
sample_2018_svanvik = fitSpl_dmean_svanvik(svanvik_daily_2018.dropna().index.dayofyear)


# Draw samples for Svanvik
x_sample, pdf, fit, stat = fit_skew_normal((svanvik_daily.dropna()-sample).values)
x_sample_2018, pdf_2018, fit_2018, stat_2018 = fit_skew_normal((svanvik_daily_2018.dropna()-sample_2018).values)
x_sample_svanvik, pdf_svanvik, fit_svanvik, stat_svanvik = fit_skew_normal((svanvik_daily_2018.dropna()-sample_2018_svanvik).values)

# Select 2018 data -> fig10
esrange_daily_2018 = data['Esrange']['2018'].resample('1d').apply(np.nanmean)
pallas_daily_2018 = data['Pallas']['2018'].resample('1d').apply(np.nanmean)
prestebakke_daily_2018 = data['Prestebakke']['2018'].resample('1d').apply(np.nanmean)

# Sample accordingly from climatology
sample_2018_esrange = fitSpl_dmean(esrange_daily_2018.dropna().index.dayofyear)
sample_2018_pallas = fitSpl_dmean(pallas_daily_2018.dropna().index.dayofyear)
sample_2018_prestebakke = fitSpl_dmean_prestebakke(prestebakke_daily_2018.dropna().index.dayofyear)

# Sample only from June-September
esrange_jja = esrange_daily_2018.where((esrange_daily_2018.index.month>=6) & (esrange_daily_2018.index.month<9)).dropna()
pallas_jja = pallas_daily_2018.where((pallas_daily_2018.index.month>=6) & (pallas_daily_2018.index.month<9)).dropna()
prestebakke_jja = prestebakke_daily_2018.where((prestebakke_daily_2018.index.month>=6) & (prestebakke_daily_2018.index.month<9)).dropna()
svanvik_jja = svanvik_daily_2018.where((svanvik_daily_2018.index.month>=6) & (svanvik_daily_2018.index.month<9)).dropna()

sample_jja_esrange = fitSpl_dmean(esrange_jja.index.dayofyear)
sample_jja_pallas = fitSpl_dmean(pallas_jja.index.dayofyear)
sample_jja_prestebakke = fitSpl_dmean_prestebakke(prestebakke_jja.index.dayofyear)
sample_jja_svanvik = fitSpl_dmean_svanvik(svanvik_jja.index.dayofyear)

# Fit the distributions
x_sample_esrange, pdf_esrange, fit_esrange, stat_esrange = fit_skew_normal((esrange_daily_2018.dropna()-sample_2018_esrange).values)
x_sample_pallas, pdf_pallas, fit_pallas, stat_pallas = fit_skew_normal((pallas_daily_2018.dropna()-sample_2018_pallas).values)
x_sample_prestebakke, pdf_prestebakke, fit_prestebakke, stat_prestebakke = fit_skew_normal((prestebakke_daily_2018.dropna()-sample_2018_prestebakke).values)

# Spectral analysis
from scipy import fftpack
fft_barrow = fftpack.fft(data["Barrow"].resample('1M').mean().fillna(method='ffill'))
freqs_barrow = fftpack.fftfreq(len(fft_barrow))

fft_prestebakke = fftpack.fft(data['Prestebakke'].resample('1M').mean().fillna(method='ffill'))
freqs_prestebakke = fftpack.fftfreq(len(fft_prestebakke))

fft_jergkara = fftpack.fft(data_jergkara.resample('1M').mean().fillna(method='ffill'))
freqs_jergkara = fftpack.fftfreq(len(fft_jergkara))


execfile("plot_ozone_observations.py")
