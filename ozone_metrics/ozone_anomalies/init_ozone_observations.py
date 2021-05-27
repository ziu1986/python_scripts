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
src_rra = os.environ['DATA']+'/nird_data/reanalysis/Copernicus/ensemble_ozone/SCA_ENSa.2018.O3.yearlyrea.nc'

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
        data_svanvik_OzoNorClim.append(tmp_data_svanvik['O3_mugm-3'].where(tmp_data_svanvik['O3_mugm-3']>=0.5).dropna()/2.)
    # Concat data Svanvik data
    data_svanvik_OzoNorClim = pd.concat(data_svanvik_OzoNorClim)
    # Load regional model reanalysis 2018 and set time axis
    data_rra = xr.open_dataset(src_rra)
    data_rra['time'] = pd.date_range("2018-01-01", periods=365*24, freq='H')
    

#---------------------------------------------------------------------------------------------------------------------------------
#Control
plot_timeseries = False
plot_timelag = False
plot_correlation = False
plot_splines = False
plot_climatology = True
plot_spectrum = False
plot_map = False
plot_aot = False
plot_cuo = False
plot_rollingsum = False
plot_residuals = False
plot_ttest = False
plot_svanvik = False

#---------------------------------------------------------------------------------------------------------------------------------


if plot_timelag:
    execfile("ozone_observation_timelag.py")
if plot_spectrum:
    execfile("ozone_obsercation_spectrum.py")

execfile("ozone_observation_histograms.py")

execfile("ozone_observation_climatology.py")

execfile("ozone_observation_sampling.py")

execfile("ozone_observation_ttest.py")

execfile("ozone_observation_fits.py")

if plot_cuo:
    execfile('ozone_observation_gsto.py')
    execfile("ozone_observation_cuo.py")

execfile("plot_ozone_observations.py")
