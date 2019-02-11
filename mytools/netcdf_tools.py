import os # Access environment variables
import glob
#import datetime as dt  # Python standard library datetime  module
import numpy as np
import xarray as xr
from scipy.constants import * # Pysical constants, units conversion
import datetime as dt

# -----------------------------------------------------------------------------

def average_ds(self, dim=None, weights=None):
    """
    weighted average for Datasets

    Parameters
    ----------
    dim : str or sequence of str, optional
        Dimension(s) over which to apply average.
    weights : DataArray
        weights to apply. Shape must be broadcastable to shape of data.

    Returns
    -------
    reduced : Dataset
        New Dataset with average applied to its data and the indicated
        dimension(s) removed.

    """

    if weights is None:
        return self.mean(dim)
    else:
        return self.apply(average_da, dim=dim, weights=weights)

def average(data, dim=None, weights=None):
    """
    weighted average for xr objects

    Parameters
    ----------
    data : Dataset or DataArray
        the xr object to average over
    dim : str or sequence of str, optional
        Dimension(s) over which to apply average.
    weights : DataArray
        weights to apply. Shape must be broadcastable to shape of data.

    Returns
    -------
    reduced : Dataset or DataArray
        New xr object with average applied to its data and the indicated
        dimension(s) removed.

    """

    if isinstance(data, xr.Dataset):
        return average_ds(data, dim, weights)
    elif isinstance(data, xr.DataArray):
        return average_da(data, dim, weights)
    else:
        raise ValueError("date must be an xr Dataset or DataArray")


def read_data(src,**karg):
    '''
    Extract data from netCDF files.
    '''
    datatype = karg.pop('datatype','')
    variable = karg.pop('var','')
    level = karg.pop('lev', 'all')
    data_list = []
    if len(sorted(glob.glob(src))) == 0:
        print("Error: No file provided")
        return
    for file in sorted(glob.glob(src)):
        print("Reading %s" % (os.path.basename(file)))
        data = xr.open_dataset(file)
        if datatype == 'osloctm':
            # Defining new time coordinates
            data['time'].reset_coords(drop=True)
            data.coords['time'] = ([dt.datetime(data['YEAR'], data['MONTH'], 15),])
        elif datatype=='openifs':
            data.coords['time'] = dt.datetime.strptime(str(data['data_date'].data),"%Y%m%d%H")
        elif datatype=='osloctm_ozone':
            # Defining new time coordinates
            year = (int)(file[file.rfind('/')+4:-7])
            days_since = (int)(file[file.rfind('/')+9:-3])-1
            hour = (0,6,12,18)
            new_time = [dt.datetime(year,1,1,ihour)+dt.timedelta(days_since) for ihour in hour]
            data['time'].reset_coords(drop=True)
            data.coords['time'] = (new_time)
        if variable == 'NOx':
            try:
                data = data['NOx']
                
            except:
                unit = data['NO3'].attrs['units']
                data = data['NO3'] + data['N2O5'] + data['NO'] + data['NO2']
                data.attrs['units'] = unit
            
        elif variable != '':
            data = data[variable]
         
        if level=='all':
            data_list.append(data)
        else:
            data_list.append(data.sel(lev=slice(level[0],level[-1])))
    # Return the concatenated data
    if len(data_list) > 1:
        data_list = xr.concat(data_list, dim='time')
    else:
        data_list = data_list[0]
    return(data_list)
    
