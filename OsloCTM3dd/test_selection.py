import os, glob, sys
import numpy as np
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt

nc_src = os.environ['DATA']+'/CTM3_input_data/EMIS/CEDS0517/NOx-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_200001-201412.nc'

data = xr.open_dataset(nc_src)

test_data = data.sel(lon=slice(-10,55),lat=slice(30,85),drop=True)
#test_data.to_netcdf("test_selection.nc")

