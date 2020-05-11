import os, sys, glob
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.station_info import *

def buck_vapor_press(temperature):
    coeff_a = 0.61121
    coeff_b = 18.678
    coeff_c = 234.5
    coeff_d = 257.14

    press = coeff_a * np.exp((coeff_b-temperature/coeff_c)*(temperature/(coeff_d+temperature)))
    # Return saturation water vapor pressure in Pa
    return(press*1e3)

def VPD(relhum, temperature):
    Ps = buck_vapor_press(temperature)
    vpd = Ps * (1-relhum/100)
    # Return VPD in Pa
    return(vpd)

# Clean up
plt.close('all')

# Source
src_svanvik_clim = os.environ['DATA']+'/astra_data/observations/Svanvik_press_relhum_wind_1987-2019.csv'

# Load data
data_list = []
try:
    data_svanvik_clim
except NameError:
    data_svanvik_clim = pd.read_csv(src_svanvik_clim, na_values=(-9999,-9999.0,'x'), parse_dates=[[1,2,3,4]], date_parser=lambda y,m,d,h : pd.datetime(int(y), int(m), int(d), int(h)-1))
    data_svanvik_clim.index = data_svanvik_clim.iloc[:,0]
    data_svanvik_clim = data_svanvik_clim.drop(columns=[u'Year_Mnth_Date_Time(NMT)', u'St.no'])

    data_svanvik_2018_2019 = data_svanvik_clim['2018':'2019'].copy()
    data_svanvik_2018_2019.loc[:,'day'] = data_svanvik_2018_2019.index.day.values
    data_svanvik_2018_2019.loc[:,'month'] = data_svanvik_2018_2019.index.month.values
    data_svanvik_2018_2019.loc[:,'hour'] = data_svanvik_2018_2019.index.hour.values

    data_svanvik_clim.loc[:,'day'] = data_svanvik_clim.index.day.values
    data_svanvik_clim.loc[:,'month'] = data_svanvik_clim.index.month.values
    data_svanvik_clim.loc[:,'hour'] = data_svanvik_clim.index.hour.values

svanvik_clim = data_svanvik_clim[:'2017'].groupby(['month','day','hour']).mean()

vpd_clim = VPD(svanvik_clim['UU'], svanvik_clim['TA'])
vpd_2018 = VPD(data_svanvik_2018_2019['2018']['UU'], data_svanvik_2018_2019['2018']['TA'])
vpd_2019 = VPD(data_svanvik_2018_2019['2019']['UU'], data_svanvik_2018_2019['2019']['TA'])

wind_clim = svanvik_clim[['FF','FF2']]
wind_2018 = data_svanvik_2018_2019['2018'][['FF', 'FF2']]
wind_2019 = data_svanvik_2018_2019['2019'][['FF', 'FF2']]



# Save data to file
vpd_2018.reindex(data_svanvik_clim['2018'].index).to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_vpd_2018.csv")
vpd_2019.reindex(data_svanvik_clim['2019'].index).to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_vpd_2019.csv")
# Drop Feb 29 from climatology and reindex it
save_data = pd.DataFrame({'VPD (Pa)':vpd_clim.drop(vpd_clim.loc[2,29,:].index).values}, index=data_svanvik_clim['2018'].index)
save_data.to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_vpd-climatology.csv")

wind_2018.reindex(data_svanvik_clim['2018'].index).to_csv("svanvik_wind_2018.csv")
wind_2019.reindex(data_svanvik_clim['2019'].index).to_csv("svanvik_wind_2019.csv")
# Drop Feb 29 from climatology and reindex it
save_data = pd.DataFrame({'U10m (ms^-1)':wind_clim.drop(wind_clim.loc[2,29,:].index)['FF'].values, 'U2m (ms^-1)':wind_clim.drop(wind_clim.loc[2,29,:].index)['FF2'].values}, index=data_svanvik_clim['2018'].index)
save_data.to_csv("svanvik_wind-climatology.csv")

