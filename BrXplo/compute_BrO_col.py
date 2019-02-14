import os, glob
import numpy as np
import pandas as pd
import xarray as xr
from mytools.med_tools import compute_column

dataset = 'BrXplo_mysic_rsnow'
year = 2000
month = np.arange(9,13)
# Access environment variable for directory
nc_data = os.environ['DATA']
# Data directory
#nc_subd = '/BrXplo'
nc_subd = '/ic2'

for imonth in month:
    nc_src = '/%s/%s*%s%s01_????_tr_bromine.nc' % (dataset, dataset[:-4], year, str(imonth).zfill(2))
    print "Reading file(s): "
    print nc_data+nc_subd+nc_src
    data = xr.open_dataset(sorted(glob.glob(nc_data+nc_subd+nc_src))[0])
    press = data['hyam']+data['hybm']*data['aps']
    # Compute BrO column density
    BrO_col = compute_column(press, data['BrO'])
    # Set some meta-data
    BrO_col.name = 'BrO_column'
    BrO_col.attrs['units'] = 'molecules/cm2'
    # Write
    BrO_col.to_netcdf("BrO_col_%s%s_%s.nc" % (year, str(imonth).zfill(2), dataset))
    # Close data file
    data.close()
