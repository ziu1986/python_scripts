import os, sys, glob
import numpy as np
import xarray as xr
from mytools.met_tools import *

src = os.environ['DATA']+"/nird_data/reanalysis/Copernicus/ensemble_ozone/SCA_ENSa.2013*"

try:
    data
except NameError:
    data = xr.open_dataset(glob.glob(src)[0])
    # Reset time
    year = (data.time/10000).astype(int)
    month = (data.time/100-year*100).astype(int) 
    day = (data.time-year*10000-month*100).astype(int)
    hour = (data.time-year*10000-month*100-day)*24
    new_time = [pd.datetime(year[i],month[i],day[i], hour[i]) for i in range(year.size)]
    data["time"] = new_time

    
