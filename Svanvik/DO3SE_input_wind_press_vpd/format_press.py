import os, sys, glob
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.station_info import *

# Clean up
plt.close('all')

b_save = True

# Source
src_svanvik_clim = os.environ['DATA'] + '/astra_data/observations/Kirkenes_airport_press_*.csv'
src_svanvik = os.environ['DATA'] + '/astra_data/observations/metdata_svanvik/Svanvik_press_2013-2020.xls'

# Load data
data_list = []
try:
    data_svanvik_clim
except NameError:
    for file in sorted(glob.glob(src_svanvik_clim)):
        data_svanvik_clim = pd.read_csv(file, na_values=(-9999,-9999.0,'x'), parse_dates=[[1,2,3,4]], date_parser=lambda y,m,d,h : pd.datetime(int(y), int(m), int(d), int(h)))
        data_svanvik_clim.index = data_svanvik_clim.iloc[:,0]
        data_svanvik_clim = data_svanvik_clim.drop(columns=[u'Year_Mnth_Date_Time(UTC)', u'St.no'])

        data_list.append(data_svanvik_clim)
    
    data_svanvik_clim = pd.concat(data_list)

    data_svanvik = pd.read_excel(src_svanvik, sheet_name=[2,3])
    data_svanvik_2018 = data_svanvik[2].copy()
    data_svanvik_2018.index = data_svanvik_2018.iloc[:,0]
    data_svanvik_2018 = data_svanvik_2018.drop(columns=['Fra-tid', 'Til-tid'])['2018']
    data_svanvik_2019 = data_svanvik[3].copy()
    data_svanvik_2019.index = data_svanvik_2019.iloc[:,0]
    data_svanvik_2019 = data_svanvik_2019.drop(columns=['Fra-tid', 'Til-tid'])['2019']
    data_svanvik_2018_2019 = pd.concat((data_svanvik_2018, data_svanvik_2019))
    #data_svanvik_2018_2019 = data_svanvik_clim['2018':'2019'].copy()
    #data_svanvik_2018_2019.loc[:,'day'] = data_svanvik_2018_2019.index.day.values
    #data_svanvik_2018_2019.loc[:,'month'] = data_svanvik_2018_2019.index.month.values
    #data_svanvik_2018_2019.loc[:,'hour'] = data_svanvik_2018_2019.index.hour.values

    data_svanvik_clim.loc[:,'day'] = data_svanvik_clim.index.day.values
    data_svanvik_clim.loc[:,'month'] = data_svanvik_clim.index.month.values
    data_svanvik_clim.loc[:,'hour'] = data_svanvik_clim.index.hour.values

svanvik_clim = data_svanvik_clim[:'2017'].groupby(['month','day','hour']).mean()


# Save data to file
if b_save:
    data_svanvik_2018_2019['2018'].iloc[:,0].to_csv("svanvik_press_2018.csv")
    data_svanvik_2018_2019['2019'].iloc[:,0].to_csv("svanvik_press_2019.csv")
    # Drop Feb 29 from climatology and reindex it
    save_data = pd.DataFrame({'Pressure (hPa)':svanvik_clim['PO'].drop(svanvik_clim.loc[2,29,:].index).values}, index=data_svanvik_clim['2019']['PO'].index)
    save_data.to_csv("svanvik_press-climatology.csv")


