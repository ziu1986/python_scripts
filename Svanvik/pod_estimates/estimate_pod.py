import os, sys, glob
import matplotlib.pyplot as plt
import datetime as dt
import pandas as pd
#from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.station_info import station_location
from mytools.plot_tools import *

#Load ozone data
def read_date(input_files):
    data_list = []
    print("Reading")
    for file_name in sorted(glob.glob(input_files)):
        print(file_name)
        data = xr.open_dataset(file_name)
        data_list.append(data)
        

    # Concat list
    data = xr.concat(data_list, dim='time')
    return(data)

try:
    data_macc
except NameError:
    #nc_src = os.environ['DATA']+"/astra/ECMWF/MACC_reanalysis/netcdf/VMR/*climatology.nc"
    nc_src = os.environ['DATA']+"/vmr_macc_r_o3_ml60_3h_climatology.nc" 
    data_macc = read_date(nc_src)
    data_macc = (data_macc['go3']).sel(latitude=station_location["Svanvik"].lat, longitude=station_location["Svanvik"].lon, method='nearest')
    #nc_src = os.environ['DATA']+"/astra/ECMWF/CAMS_reanalysis/netcdf/VMR/*climatology.nc"
    nc_src = os.environ['DATA']+"/vmr_cams_r_o3_ml60_climatology.nc" 
    data_cams = read_date(nc_src)
    data_cams = (data_cams['go3']).sel(latitude=station_location["Svanvik"].lat, longitude=station_location["Svanvik"].lon, method='nearest')

    #nc_src = os.environ['DATA']+"/nird/reanalysis/Copernicus/ensemble_ozone/*climatology.nc"
    nc_src = os.environ['DATA']+"/SCA_ENSa.yearlyrea.climatology.nc"
    data_rra = read_date(nc_src)
    data_rra['time'] = pd.date_range("2016-01-01", periods=366*24, freq='H')
    data_rra = (data_rra['O3']*0.5*1e-9).sel(lat=station_location["Svanvik"].lat, lon=station_location["Svanvik"].lon, method='nearest')

    #src = os.environ['DATA']+'/astra/input_data/DO3SE_input/svanvik_ozone_2018.csv'
    src = os.environ['DATA']+'/astra/input_data/DO3SE_input/svanvik_ozone_climatology.csv'
    data_svanvik = pd.read_csv(src)
    data_svanvik.index = pd.to_datetime(data_svanvik['Year_Mnth_Date_Time(NMT)'].values, format='%Y-%m-%d %H:%M:%S')
    data_svanvik = data_svanvik.drop(columns=["Year_Mnth_Date_Time(NMT)"])



from func_def_plot_javis import *

# Data sources
src_svanvik = os.environ['DATA']+'/astra/observations/metdata_svanvik/Svanvik_temp_relhum_wind_*.csv'
src_svanvik_rad = os.environ['DATA']+'/astra/observations/metdata_svanvik/svanvik_glob_rad_*.csv'   
    
# Set up the different species
# Evergreen
evergreen = JavisModel('coniferous', Tmin=0, Tmax=200, Topt=20, fmin=0.1, Dmax=0.8, Dmin=2.8, alpha=0.006, gmax=125)
# Birch
birch = JavisModel('deciduous', Tmin=5, Tmax=200, Topt=20, fmin=0.1, Dmax=0.5, Dmin=2.7, alpha=0.0042, gmax=240)
# Grassland
grassland = JavisModel('grassland', Tmin=10, Tmax=36, Topt=24, fmin=0.1, Dmax=1.75, Dmin=4.5, alpha=0.011, gmax=190)

# idata.resample('1H').mean().interpolate('linear').
shift_macc = data_macc[1:].to_pandas().shift(periods=(6*365+1)*24,freq='H')
shift_macc = shift_macc.resample('1H').mean().interpolate('linear')
shift_cams = data_cams[1:].to_pandas().shift(periods=(6*365+1)*24,freq='H')
shift_cams = shift_cams.resample('1H').mean().interpolate('linear')
shift_camsraq = data_rra.to_pandas().shift(periods=(2*365+1)*24,freq='H')
# Load data
try:
    data_temp
except NameError:
    data_temp = import_data(src_svanvik)
    data_rad = import_data(src_svanvik_rad)
# Comput vpd
vpd = VPD(data_temp.iloc[:,1], data_temp.iloc[:,0])/kilo

# P/(RT) at std P and 20deg C to convert ozone ppb -> mmol m-3
conversion =  1.013e5/((data_temp.iloc[:,0]+273.15)*8.31447)['2018-05':'2018-09']*1e-6
results_pod = []
results_pod1 = []
for species in (birch, evergreen, grassland):
    for idata, iname in zip((data_svanvik['O3 (ppb)'], shift_macc*1e9, shift_cams*1e9, shift_camsraq*1e9),("obs", "macc", "cams", "camsraq")):
        tmp = get_f_function(species, data_temp, vpd, data_rad)[-1]['2018-05':'2018-09']*species.gmax*(idata*conversion)/60**2
        tmp1 = tmp.apply(lambda x: x-1e-6)
        print("POD0", species.name, iname, tmp.sum())
        results_pod.append(tmp.sum())
        print("POD1", species.name, iname, tmp1[tmp1 > 0].sum())
        results_pod1.append(tmp1[tmp1 > 0].sum())

for i in np.arange(1,4):
    print((np.array(results_pod[i::4])-np.array(results_pod[0::4]))/np.array(results_pod[0::4])*100)
    print((np.array(results_pod1[i::4])-np.array(results_pod1[0::4]))/np.array(results_pod1[0::4])*100)


    

