import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *
from collections import OrderedDict
# Clean up
plt.close('all')
species = 'O_3'
names = {'O_3':'ozone', 'HNO_3':'nitric_acid', 'SO_2':'sulfur_dioxid', 'NH_3':'ammonia', 'SO_4':'sulfate'}
title_name = names[species]
# Data source
b_ozone = True
scav_dir = "scavenging_daily/"
experiment = (#'C3RUN_default/',
              #'C3RUN_mOSaic/',
              #'C3RUN_mOSaic_offLight/',
              #'C3RUN_mOSaic_offPhen/',
              #'C3RUN_mOSaic_SWVL1/',
              #'C3RUN_mOSaic_ice/',
              #'C3RUN_mOSaic_desert/',
              'C3RUN_mOSaic_emis2014/',
    'C3RUN_mOSaic_hough/'
              #'C3RUN_oDD/',
              #'C3RUN_emep_full/',
              #'C3RUN_emep_offLight/',
              #'C3RUN_emep_offPhen/',
              #'C3RUN_emep_SWVL4/',
              #'C3RUN_emep_ppgs/',
              #'C3RUN_emep_ppgssh/',
              #'C3RUN_emep_ppgssh_ice/',
              #'C3RUN_emep_ppgs_2005/'
)

data_dir = os.environ['DATA']+'/astra_data/ctm_results/' 
    
labels = (#'Wesely_type',
          #'mOSaic',
          #'mOSaic_offLight','mOSaic_offPhen','mOSaic_SWVL1',
          #'mOSaic_ice',
          #'mOSaic_desert',
          'mOSaic_emis2014',
    'mOSaic_hough'
    #'OsloCTM3: Wesely type',
          #'OsloCTM3: EMEP_swgd',
          #'OsloCTM3: EMEP_full',
          #'OsloCTM3: EMEP_offLight','OsloCTM3: EMEP_offPhen','OsloCTM3: EMEP_SWVL4',
          #'OsloCTM3: EMEP_ppgs','OsloCTM3: EMEP_ppgssh',
          #'OsloCTM3: EMEP_ppgssh_ice',
          #'OsloCTM3: EMEP_ppgs_2005'
)
colors = np.concatenate((('black',),
                         ('blue',),
                         np.repeat('cornflowerblue',3),
                         np.repeat('purple',2),
                         ('grey',),
                         ('orange',)) )


for iexp in experiment:
    subdir = data_dir+iexp+scav_dir+'scavenging_daily_2d_*.nc'
    print(data_dir+iexp+scav_dir)
    for file in sorted(glob.glob(subdir)):
        newfile_name = "vo3"+os.path.basename(file)[10:]
        print(newfile_name)
        data = xr.open_dataset(file)
        v_avg = 1/((1/data['VRaO3_avg'].where(data['VRaO3_avg']>0))+(1/data['VRbO3_avg'].where(data['VRbO3_avg']>0))+(1/data['VRcO3_avg'].where(data['VRcO3_avg']>0)))
        output = xr.Dataset({'VO3':v_avg})
        output.to_netcdf(data_dir+iexp+scav_dir+newfile_name)
