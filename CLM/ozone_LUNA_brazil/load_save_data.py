import os, sys, glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all, plot_error_bands
from mytools.clm_tools import *

def load_date():
    # Load reference simulation
    brazil_src = run_archive + case[0] + land_hist
    brazil_test.update({threshold[0]:load_data(brazil_src)})
    brazil_test_40.update({threshold[0]:load_data(brazil_src.replace('100', '40'))})
    brazil_test_20.update({threshold[0]:load_data(brazil_src.replace('100', '20'))})
    brazil_test_60.update({threshold[0]:load_data(brazil_src.replace('100', '60'))})
    brazil_test_80.update({threshold[0]:load_data(brazil_src.replace('100', '80'))})
    brazil_test_ozone.update({ozone[0]:load_data(brazil_src)})
    brazil_test_ozone.update({ozone[1]:load_data(brazil_src.replace('_ozone_luna_100', ''))})
    brazil_ref_ozone.update({ozone[0]:load_data(brazil_src.replace('_luna_100', ''))})
    brazil_ref_ozone.update({ozone[1]:load_data(brazil_src.replace('_ozone_luna_100', ''))})

    # Load sensitivity tests
    for ithresh in threshold[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test.update({ithresh:load_data(brazil_src)})
    for ithresh in threshold_2[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test_40.update({ithresh:load_data(brazil_src.replace('100', '40'))})
    for ithresh in threshold_2[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test_20.update({ithresh:load_data(brazil_src.replace('100', '20'))})
    for ithresh in threshold_2[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test_60.update({ithresh:load_data(brazil_src.replace('100', '60'))})
    for ithresh in threshold_2[1:]:
        brazil_src = run_archive + case[1] + "%s" % ithresh + land_hist
        brazil_test_80.update({ithresh:load_data(brazil_src.replace('100', '80'))})


    for iozone in ozone[2:]:
        brazil_src = run_archive + case[2] + "%s" % iozone + land_hist
        brazil_test_ozone.update({iozone:load_data(brazil_src)})
        brazil_ref_ozone.update({iozone:load_data(brazil_src.replace('luna_',''))})
    for iozone in ozone_deep:
        brazil_src = run_archive + case[2] + "%s" % iozone + land_hist
        brazil_test_ozone_deep.update({iozone:load_data(brazil_src)})

    # Merge fine resolution scan ozone deep
    brazil_test_ozone.update(brazil_test_ozone_deep)

def save_data(**karg):
    directory = karg.pop('dir', os.environ['DATA'] + "/preprocessed_data/CLM50_ozone_luna_brazil")
    for ithr in threshold:
        brazil_test[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_100_%s.nc" % (directory, ithr))
        brazil_test[ithr][['NPP','GPP']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_100_%s.nc" % (directory, ithr))
    for ithr in threshold_2:
        brazil_test_40[ithr][['NPP','GPP']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_40_%s.nc" % (directory, ithr))
        brazil_test_20[ithr][['NPP','GPP']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_20_%s.nc" % (directory, ithr))
        brazil_test_60[ithr][['NPP','GPP']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_60_%s.nc" % (directory, ithr))
        brazil_test_80[ithr][['NPP','GPP']].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_80_%s.nc" % (directory, ithr))
        brazil_test_40[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_40_%s.nc" % (directory, ithr))
        brazil_test_20[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_20_%s.nc" % (directory, ithr))
        brazil_test_60[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_60_%s.nc" % (directory, ithr))
        brazil_test_80[ithr].groupby('time.month').mean().to_netcdf("%s/brazil_80_%s.nc" % (directory, ithr))
    for ioz in ozone:
        brazil_ref_ozone[ioz].apply(lambda x: x.groupby('time.month').sum()/np.unique(x.time.dt.year).size).to_netcdf("%s/brazil_npp_%s_ref.nc" % (directory,ioz))
        brazil_ref_ozone[ioz].groupby('time.month').mean().to_netcdf("%s/brazil_%s_ref.nc" % (directory,ioz))

# Source
ref_data_src = os.environ['PY_SCRIPTS'] + '/plant_model/test.cvs'
try:
    run_archive = os.environ['CESM_RUN'] + '/archive/'
except KeyError:
    run_archive = os.environ['DATA'] + '/astra_data/clm_results/'
    
land_hist = '/lnd/hist/*.clm2.h0.*.nc'
case = ('brazil_2000_ozone_luna_100', 'brazil_2000_ozone_luna_100_thresh_', 'brazil_2000_ozone_luna_')
threshold = (0.8, 0.0, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.7, 0.85, 0.9, 1, 2, 3, 4, 5)
threshold_2 = (0.8, 0, 0.2, 0.5, 1, 2, 3, 4, 5)
ozone = (100, 0, 40, 60, 80)
ozone_deep = np.arange(42, 60, 2)

# Load data

ref_data = pd.read_csv(ref_data_src)

brazil_test = {}
brazil_test_40 = {}
brazil_test_20 = {}
brazil_test_60 = {}
brazil_test_80 = {}
brazil_test_ozone = {}
brazil_test_ozone_deep = {}
brazil_ref_ozone = {}

load_date()
