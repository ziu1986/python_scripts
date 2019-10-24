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

def covariance(x, y, dims=None):
    return xr.dot(x - x.mean(dims), y - y.mean(dims), dims=dims) / x.count(dims)

def corrrelation(x, y, dims=None):
    return covariance(x, y, dims) / (x.std(dims) * y.std(dims))

from scipy import stats
def fit_skew_normal(data, **karg):
    '''
    Fit a Fit a skew normal distribution. Only works with plain numpy arrays.
    '''
    b_stats = karg.pop('stats', True)
    b_fit = karg.pop('fit', True)
    ae, loce, scalee = stats.skewnorm.fit(data)
    fit = {"shape":ae, "location":loce, "scale":scalee}
    max_x = data.max().round()
    min_x = data.min().round()
    sample_x = np.linspace(min_x, max_x)
    p = stats.skewnorm.pdf(sample_x, ae, loce, scalee)
    mean, variance, skewness, kurtosis = stats.skewnorm.stats(ae, loce, scalee, moments='mvsk')
    median = stats.skewnorm.median(ae, loce, scalee)
    # Calculate mode of the skew norm
    delta = ae/np.sqrt(1+ae**2)
    muz = np.sqrt(2/pi)*delta
    sigmaz = np.sqrt(1-muz**2)
    mo = muz - skewness*sigmaz*0.5-np.sign(ae)*0.5*np.exp(-2*pi/np.abs(ae))
    mode = loce + scalee*mo
    stat = {"mean":mean, "variance":variance, "skew":skewness, "kurtosis":kurtosis, "median":median, "mode":mode}
    print("Fit skew normal:\n shape = %3.1f\n location = %3.1f\n scale = %3.1f" % (ae, loce, scalee))
    print("Stats:\n mean = %3.1f\n variance = %3.1f\n skewness = %3.1f\n mode = %3.1f" % (mean, variance, skewness, mode))
    if (b_stats and b_fit):
        return(sample_x, p, fit, stat)
    elif b_stats:
        return(sample_x, p, stat)
    elif b_fit:
        return(sample_x, p, fit)
    else:
        return(sample_x, p)

    
def stats_text(ax, stat, fit, **karg):
    '''
    Set a stats box for the fit.
    '''
    xpos = karg.pop("xpos", 0.01)
    ypos = karg.pop("ypos", 0.9)
    name = karg.pop('name',"statistics")
    name = name.replace(" ", "\,")
    quantile_interval = stats.skewnorm.interval(0.5, fit["shape"], fit["location"], fit["scale"])
    std = stats.skewnorm.std(fit["shape"], fit["location"], fit["scale"])
    
    ax.text(xpos, ypos, "$%s$\n mean = %3.1f\n std = %3.1f\n median = %3.1f\n (q1, q3) = (%3.1f, %3.1f)\n mode =  %3.1f" % (name, stat["mean"], std, stat["median"], quantile_interval[0], quantile_interval[1],  stat["mode"]), transform=ax.transAxes)


