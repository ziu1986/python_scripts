import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *

def read_infile(filename):
    print("Reading: %s" % (os.path.basename(filename)))
    data = xr.open_dataset(filename)
    # Define new time coordinates and drop the old one
    data['time'].reset_coords(drop=True)
    year = int(os.path.basename(filename)[3:7])
    doy = int(os.path.basename(filename)[-6:-3])
    data.coords['time'] = ([dt.datetime(year, 1, 1, 0) + dt.timedelta(doy - 1, i*60**2) for i in np.arange(0,24,24/len(data['time']))])
   
    return(data)

# Data source
project = '/abel/C3RUN_emep_SWVL4.231018.32284/'
nc_src = os.environ['DATA']+project
sub_dir = ("trop_tracer/trp*.nc","air_density/air*.nc")
molweights_src = os.environ['DATA']+project+'monthly_means/avgsav_20050101_20050201.nc'
Mair = 28.949 # g/mol
# Read the data
try:
    molweights
except NameError:
    molweights = get_molarweight(xr.open_dataset(glob.glob(molweights_src)[0]))

for sub_i, sub_j in zip(sorted(glob.glob(nc_src+sub_dir[0])),sorted(glob.glob(nc_src+sub_dir[1]))):
    data = (read_infile(sub_i), read_infile(sub_j))      
    mair = data[1]['air_densit'].isel(lev=0) # AIRdnst
    Mtracer = molweights.sel(name="O3")
    mtracer = data[0]['O3'].isel(lev=0)
          
    vmr = mtracer/mair*Mair/Mtracer
    vmr.attrs['units'] = 'mol/mol'

    vmr_over40 = vmr/1e-9-40
    vmr_over40.attrs['units'] = 'ppb'
    vmr_over40 = vmr_over40.where(vmr_over40>=0)

    dataset = xr.Dataset({'AOT40':vmr_over40, 'O3':vmr})
    outfile = "aot40_%s" % (os.path.basename(sub_i)[3:])
    print("Writing: %s" % (outfile))
    dataset.to_netcdf(outfile)

