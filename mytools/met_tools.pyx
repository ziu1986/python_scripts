import os
import numpy as np
import pandas as pd
import calendar
from matplotlib.colors import Normalize
from scipy.constants import * # Pysical constants, units conversion
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d, interp2d
import datetime as dt  # Python standard library datetime  module
import matplotlib.pyplot as plt
import cartopy as cp        # Globe projections

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

# USstd atmosphere 1976: Pressure from height
def USstdP(height):
    '''
    USstd atmosphere 1976: Pressure from height
    '''
    if height<=11*kilo:
        return(atm*(1-height/44329.)**5.255876);
    elif height<=20*kilo:
        return(atm*(0.223361)*np.exp((10999.-height)/6341.4));
    elif height<=32*kilo:
        return(atm*(0.988626+height/198903.)**-34.16319);
    elif height<=47*kilo:
        return(atm*(0.898309+height/55280.)**-12.20114);
    elif height<=51*kilo:
        return(atm*(0.00109456)*np.exp((46998.-height)/7922.));
    elif height<=86*kilo:
        return(atm*(0.838263-height/176142.)**12.20114);

def draw_altitude_axis(ax):
    '''
    Add additional y-axis. Altitude computed as reverse function of height->pressure relationship in USstd atmosphere.
    The pressure-axis ax should be in hPa!
    '''
    ax_x_interval = ax.get_xbound()
    ax_y_interval = ax.get_ybound()
    ax_alt = ax.twinx()
    ax_alt.set_ylabel("Altitude (km)")
    # Get metric meassure of the height
    altitude = np.arange(0, 86*kilo, 5) # top of USstd atmosphere parametrization
    pressure_from_altitude = np.array([USstdP(each) for each in altitude])
    fPAlt = InterpolatedUnivariateSpline(-pressure_from_altitude, altitude) # x has to be increasing! y in meters
    # Define range for the axis (x_min, x_max, y2_min, y2_max)
    ax_alt.axis([ax_x_interval[0], ax_x_interval[1], 
              fPAlt(-ax_y_interval[1]*hecto)/kilo, fPAlt(-ax_y_interval[0]*hecto)/kilo]) 
    plt.grid(0)
    return ax_alt

def set_pressure_axis(ax, **kargs):
    '''
    Adjust the y-axis of a plot to represent pressure levels.
    kargs: limits, ticks, label
    '''
    limits = sorted(kargs.pop("limits", (1000., USstdP(30*kilo)/hecto)))[::-1] # reverse sorted
    ticks = kargs.pop("ticks", [20, 50, 100, 200, 500, 1000])
    label = kargs.pop("label", "Pressure (hPa)")
    ax.set_ylabel(label)
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_yticks(ticks)
    ax.invert_yaxis()
    ax.set_ylim(limits)

def datetime_from_str(htime):
    '''
    Returns np.array of datetime objects from list of EMAC date string (yyyymmdd,Minute/(60*24)).
    Doesn't work for monthly mean.
    '''
    # Getting hour from htime
    cdef int delta_time = int(round((htime[2]-htime[1])*24))
    hour = (htime-htime.astype(int))*24 # np.array of floats
    rHound = [] # round to full hour
    cdef double iHour
    for iHour in hour:
        rHound.append(int(round(iHour)))
    cdef int minute = int(round((rHound[0]-hour[0])*60))
    d = [] # new datetime string
    cdef int i
    for i in np.arange(len(htime)):
        if rHound[i] >= 10:
            if rHound[i] == 24:
                d.append("%s 00:%s" % (htime.astype(int)[i+1], minute))
            else:
                d.append("%s %s:%s" % (htime.astype(int)[i], rHound[i], minute))
        else:
            d.append("%s 0%s:%s" % (htime.astype(int)[i], rHound[i], minute))

    x_time = np.array([dt.datetime.strptime(dTime, '%Y%m%d %H:%M') for dTime in d])
    x_time = x_time + dt.timedelta(minutes=delta_time)
    return x_time, delta_time

def datetime_from_time(start_date, time):
    '''
    Returns np.array of datetime objects from list of EMAC time (days since yyyy-mm-dd 00:00).
    Doesn't work for monthly mean.
    '''
    # Getting hour from time
    cdef int delta_time = int(round((time[2]-time[1])*24))
    hour = (time-time.astype(int))*24 # np.array of floats
    day = time.astype(int) # np.array
    rHound = [] # round to full hour
    cdef double iHour
    for iHour in hour:
        rHound.append(int(round(iHour)))
    cdef int minute = int(round((rHound[0]-hour[0])*60))
     # Generate the date vector
    dtime = np.array([ dt.timedelta(days=day[i], hours=hour[i], minutes=minute, seconds=0) for i in range(len(time)) ])
    x_time = start_date + dtime
    return x_time, delta_time

class FieldConversion:
    '''
    Convert input fields to standardized pressur levels.\n
    Expects:  \n
    atmospheric field variable, pressure ([P] = Pa); \n
    optional: userdefined evaluation "levels" ([lev] = Pa), "model"\n
    If left empty "model" falls back to ecmwf levels,
    If "levels" are provided "model" will be ignored. 
    If values of levels are outside the interpolation range nan will be written to the resulting field.\n
    Use np.nanmean()!\n
    Available models: ecmwf, gdas, ncep_r2\n
    '''
    dic_models = {"ecmwf":np.array([1,2,3,5,7]+[10,20,30,50,70]+range(100,250,25)+range(250,750,50)+range(750,1001,25)), \
                  "gdas":np.array([20,] + range(50,900,50) + range(900,1001,25)), \
                  "ncep_r2":np.array([10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000]), \
                  "escimo":np.array([  9.90000000e-03,   4.28000000e-02,   1.11200000e-01, \
                                        2.31500000e-01,   4.25900000e-01,   7.35500000e-01, \
                                        1.21530000e+00,   1.92920000e+00,   2.95710000e+00, \
                                        4.40010000e+00,   6.37480000e+00,   9.01090000e+00, \
                                        1.24507000e+01,   1.68459000e+01,   2.23544000e+01, \
                                        2.91362000e+01,   3.74306000e+01,   4.75176000e+01, \
                                        5.96433000e+01,   7.40594000e+01,   9.10194000e+01, \
                                        1.10655000e+02,   1.32390000e+02,   1.55785000e+02, \
                                        1.81158000e+02,   2.08668000e+02,   2.38335000e+02, \
                                        2.70154000e+02,   3.04024000e+02,   3.39867000e+02, \
                                        3.77551000e+02,   4.16867000e+02,   4.57703000e+02, \
                                        4.99976000e+02,   5.43502000e+02,   5.88081000e+02, \
                                        6.33529000e+02,   6.79605000e+02,   7.25949000e+02, \
                                        7.72111000e+02,   8.17539000e+02,   8.61415000e+02, \
                                        9.02643000e+02,   9.39937000e+02,   9.71661000e+02, \
                                        9.95671000e+02,   1.00935000e+03])}# hPa

    def __init__(self, **kargs):
        self.Levels = kargs.pop('levels', None)
        self.Model = kargs.pop('model', 'ecmwf')
            
    def get_model_level(self):
        '''
        Return the used model levels
        '''
        if self.Levels is None:
            return self.dic_models[self.Model]*hecto
        else:
            return self.Levels

    def __convert_vector(self, atm_field, pressure, interp_levels):
        inter_x = np.log(pressure/atm)
          
        f1 = interp1d(inter_x, atm_field, bounds_error=False)
        tmp_field = f1(np.log(interp_levels/atm))

        return (tmp_field)
    
    def convert_field(self, atm_field, pressure):
        '''
        Convert input nD-fields to standardized pressure levels.\n
        Calls __convert_vector which is a wrapper for interp1d
        
        '''
        # Interpolation for the whole grid
        # move level axis to the back
        atm_field2d = np.transpose(atm_field, (0,2,3,1)) # todo check dimensions!
        press2d = np.transpose(pressure, (0,2,3,1))
        # reshape the array 4d->2d (date x lat x lon, lev)
        shape_old = atm_field2d.shape
        shape_old_red = shape_old[:-1]
        shape_new = (np.multiply.reduce(shape_old_red), shape_old[-1])
    
        atm_field2d = np.reshape(atm_field2d, shape_new) 
        press2d = np.reshape(press2d, shape_new) 
     
        interp_levels = self.get_model_level()

        # reserve space for the converted field
        cdef int n_levels = len(interp_levels)
        cdef int n_locations = shape_new[0]
        shape_tmp = (n_locations, n_levels)
        atm_field_new = np.zeros(shape_tmp)
    
        # iterate the flattend globe
        cdef int iloc
        for iloc in np.arange(n_locations):
            atm_field_new[iloc,:] = self.__convert_vector(atm_field2d[iloc,:], press2d[iloc,:], interp_levels)
   
        # reshape and transpose it to the original shape
        shape_new2 = shape_old_red + (n_levels,)
        atm_field_new = np.reshape(atm_field_new, shape_new2)
        atm_field_new = np.transpose(atm_field_new, (0,3,1,2))
    
        return (atm_field_new)

def convert_to_std_lev(atm_field, pressure, levels=None, model="ecmwf"):
    '''
    Convert input 1D-fields to standardized pressure levels.\n
    Type convert_to_std_lev(field, pressure, levels=[], model="name")\n
    Expects:  [pressure] = Pa, field/pressure 1D 
    "levels" and "model" are optional; if left empty "model" falls back to ecmwf levels,
    if "levels" are provided "model" will be ignored. 
    If values of levels are outside the interpolation range nan will be written to the resulting field.\n
    Available models: ecmwf, gdas, ncep_r2\n
   
    '''
    Models = FieldConversion # todo
    inter_x = np.log(pressure/atm)
    inter_y = atm_field

    if levels is None:
        interp_levels = Models.dic_models[model]*hecto
    else:
        interp_levels = levels
    
    f1 = interp1d(inter_x, inter_y, bounds_error=False)
    log_tmp = f1(np.log(interp_levels/atm))
    return (log_tmp, interp_levels)
   
def convert_field(atm_field, pressure, levels=None, model='ecmwf'):
    '''
    Convert input nD-fields to standardized pressure levels.\n
    Calls convert_to_std_lev
    '''
    Models = FieldConversion # todo

    # Interpolation for the whole grid
    # move level axis to the back
    atm_field2d = np.transpose(atm_field, (0,2,3,1))
    press2d = np.transpose(pressure, (0,2,3,1))
    
    # reshape the array 4d->2d (date x lat x lon, lev)
    shape_old = atm_field2d.shape
    shape_new = (np.multiply.reduce(shape_old[:-1]), shape_old[-1])
    
    atm_field2d = np.reshape(atm_field2d, shape_new) 
    press2d = np.reshape(press2d, shape_new) 

    if levels is None:
        interp_levels = Models.dic_models[model]*hecto
    else:
        interp_levels = levels

    # reserve space for the converted field
    shape_tmp = (shape_new[0], len(interp_levels))
    atm_field_new = np.zeros(shape_tmp)
    
    # iterate the flattend globe
    for ilev in np.arange(0, shape_new[0]):
        atm_field_new[ilev,:] = convert_to_std_lev(atm_field2d[ilev,:], press2d[ilev,:])[0]
   
    # reshape and transpose it to the original shape
    shape_new2 = shape_old[:-1] + (len(interp_levels),)
    atm_field_new = np.reshape(atm_field_new, shape_new2)
    atm_field_new = np.transpose(atm_field_new, (0,3,1,2))
    return (atm_field_new, interp_levels)

def weighted_mean_std(variances, weights, axis=0):
    '''
    Calculate the standard deviation of a weighted mean. Kown as np.average.
    '''
    var = np.add.reduce(variances**2*weights**2, axis=axis)
    return np.sqrt(var)

def make_box_plot(data):
    '''
    Generate array from 4D-tracegase field to be used as boxplot.
    '''
    level = data.shape[1]
    boxData = []
    for iLevel in np.arange(level):
        tmp_box = np.reshape(np.transpose(data, axes=(1,0,2,3)), (level, -1))[iLevel, :]
        if isinstance(data, np.ma.core.MaskedArray):
            boxData.append(tmp_box[~tmp_box.mask])
        else:
            boxData.append(tmp_box)
    return np.array(boxData)

def find_nearest(a, a0, index=False):
    '''
    Element in nd-array 'a' closest to the scalar value 'a0'
    '''
    idx = np.abs(a - a0).argmin()
    if index==False:
        return a.flat[idx]
    else:
        return idx

def print_all(**kargs):
    '''
    Plots all figures. Make sure all windows have a meaningfull name. Use fig.canvas.set_window_title()!
    options (standard values): 
    target='./plots'
    file_type=('pdf', 'svg')
    '''
    fig_path = kargs.pop('target', './plots')
    file_type = kargs.pop('type', ('pdf','svg','png'))
    # Check if directory exists
    if not os.path.isdir(fig_path):
        os.mkdir(fig_path)
    for i in plt.get_fignums():
        fig = plt.figure(i)
        w_title = fig.canvas.get_window_title()
        print "Print " + w_title
        for itype in file_type:
            if not os.path.isdir(fig_path+'/'+itype):
                os.mkdir(fig_path+'/'+itype)
            fig.savefig(fig_path+'/'+itype+'/'+w_title+'.'+itype)

def zonal_band_fraction(lat_ref, lat_grid):
    '''
    Computes the right from reference latitude fraction of the zonal band.
    The fraction can be used directly for the upper boundary.
    1-fraction for lower boundary.
    '''
    cdef double lat_i = find_nearest(lat_grid, lat_ref)
    cdef int lat_idx_i = find_nearest(lat_grid, lat_ref, index=True)
    cdef double lat_left = (lat_i-lat_grid[lat_idx_i-1])/2+lat_i
    cdef double lat_right = (lat_i-lat_grid[lat_idx_i+1])/2+lat_i
    cdef double A_i = np.sin(lat_right*pi/180)-np.sin(lat_left*pi/180)
    cdef double A_ref = np.sin(lat_ref*pi/180)-np.sin(lat_left*pi/180)
    #print lat_i, lat_1, lat_2
    return A_ref/A_i

def seconds_in_month(month, year):
    '''
    Compute the corresponding seconds for a month.
    seconds_in_month(month, year)
    '''
    nomatter, daysinmonth = calendar.monthrange(year, month)
    return daysinmonth * 24 * 60 * 60
        
def plot_error_bands(ax, x, y, error, **kargs):
    '''
    Function to create a errorband like plot.
    '''
    color = kargs.pop('color', 'blue')
    ls = kargs.pop('ls', '-')
    ax.plot(x, y, color=color, ls=ls)
    ax.fill_between(x, y-error, y+error,
                    color=color, alpha=0.5)

#def compute_column(press, atm_var):
#    cdef double mdry = 0.0289645    # molec. wt. of dry air, kg/mol
#    cdef double constant = N_A/(g*mdry) # molecules*s2/(m*kg)
#    press_new = [0,]
#    atm_new = np.zeros_like(atm_var)
#    for i in range(len(press)-1):
#        press_new.append(press[i]+(press[i+1]-press[i])*0.5)
#    press_new.append(press[i+1]+(press[i+1]-press[i])*0.5)
#    for i in range(1,len(press_new)):
#        atm_new[:,i-1,:,:] = abs(press_new[i]-press_new[i-1])*atm_var[:,i-1,:,:]
#    return np.ma.add.reduce(atm_new, axis=1)*constant*1e-2/2.69e16 # 1 DU := 2.69 molecules/cm2; factor 1e-2 from units conversation
def compute_column_density(press, atm_var, **kargs):
    '''
    Compute the column density. 
    Warning: Make sure pressure is in units of hPa!
    '''
    # Keyword arguments
    unit = kargs.pop('unit', 'molecules/cm2')
    # Conversion factors
    cdef double Mdry = 0.0289645    # molec. wt. of dry air, kg/mol
    cdef double constant = N_A/(g*Mdry) # molecules*s2/(m*kg)
    cdef double cm2 = 1e-4
    cdef double DU = 1e-2/2.69e16
    # Computation
    # If press.data would not be used, labels may interfer in regrouping.
    # Shift the pressure field and multiply it with the unshifted atmospheric variable,
    # but start with second and end one before last.
    p1 = press[2:]
    p2 = press[:-2]
    if np.all(np.equal(data['O3'].shape, data['lev'].shape)):
        atm_new = 0.5*(p1-p2)*atm_var[1:-1]
    else:
        atm_new = xr.DataArray.copy(atm_var)
        for i in range(len(p1)):
            atm_new[:,i,:,:] = np.fabs(0.5*(p1[i]-p2[i]))*atm_var[:,i+1,:,:]
    # Sum it up == integration
    if unit == 'DU':
        atm_new = atm_new.sum(dim='lev')*constant*DU
    else:
        atm_new = atm_new.sum(dim='lev')*constant*cm2
    # Set the new units
    atm_new.attrs['units'] = unit
    return atm_new

def draw_parallels(ax, parallels, **kargs):
    '''
    Draw parallels in cartopy plot.
    '''
    cdef double pc = 66.57 # deg
    b_pc = kargs.pop('polarcircle', False)
    cdef int resolution = kargs.pop('resolution', 100)
    lon = np.linspace(-180, 180, resolution)
    for ilat in parallels:
        lat = np.linspace(ilat, ilat, resolution)
        ax.plot(lon, lat, color='lightgrey', linewidth=.8, transform=cp.crs.PlateCarree())
    if b_pc:
        ax.plot(lon, np.linspace(pc,pc,resolution), color='lightgrey', linewidth=.8, ls='--', transform=cp.crs.PlateCarree())
        ax.plot(lon, np.linspace(-pc,-pc,resolution), color='lightgrey', linewidth=.8, ls='--', transform=cp.crs.PlateCarree())

def draw_meridians(ax, meridians, **kargs):
    '''
    Draw meridians in cartopy plot.
    '''
    polarcap = kargs.pop('polarcap', 80)
    resolution = kargs.pop('resolution', 100)
    lat = np.linspace(-polarcap, polarcap, resolution)
    for ilon in meridians:
        lon = np.linspace(ilon,ilon, resolution)
        ax.plot(lon, lat, color='lightgrey', linewidth=.8, transform=cp.crs.PlateCarree())
    ax.plot(np.linspace(-180, 180, resolution),
            np.linspace(polarcap, polarcap, resolution),
            color='lightgrey', linewidth=.8, transform=cp.crs.PlateCarree())
    ax.plot(np.linspace(-180, 180, resolution),
            np.linspace(-polarcap, -polarcap, resolution),
            color='lightgrey', linewidth=.8, transform=cp.crs.PlateCarree())
    # Adding text
    #for ilon in np.unique(np.fabs(meridians)):
    #    ax.text(-ilon, 50, '%d$^\circ$W' % (ilon), horizontalalignment='center', transform=cp.crs.Geodetic())
import cartopy.util as ccrs_util  # Add cyclic
import xarray as xr
def addcyclicpoint(data, lon):
    '''
    Wrapps cartopy.util.add_cyclic_point() to xarray.
    Works only on 1-level data.
    '''
    cyclic_data, cyclic_lon = ccrs_util.add_cyclic_point(data.values, lon)
    try:
        cyclic_da = xr.DataArray(cyclic_data, coords=[data['time'], data['lat'], cyclic_lon], dims=['time','lat', 'lon'])
    except KeyError:
        cyclic_da = xr.DataArray(cyclic_data, coords=[data['time'], data['latitude'], cyclic_lon], dims=['time','lat', 'lon'])
    
    return cyclic_da

# Global maps
# Set the map and axis attributes
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
def plot_ozone_drydep(ax, data, **kargs):
    '''
    Plotting routine for the ozone dry deposition.
    '''
    verbose = kargs.pop('v',False)
    mode = kargs.pop('mode', 'sum')
    glob = kargs.pop('glob', True)
    region = kargs.pop('region', (-180,180,-90,90))
    title = kargs.pop('title',"")
    #cbar_kwargs = kargs.pop('cbar_kwargs', {'fraction':0.046, 'pad':0.04,'aspect':30})
    if glob:
        if verbose:
            print("Setting globe.")
        ax.set_global()   # Expands the map to fit the tick labels!
        ax.set_xticks(np.arange(-180, 181, 45), crs=cp.crs.PlateCarree())
        ax.set_yticks(np.arange(-80, 81, 20), crs=cp.crs.PlateCarree())
        lon_formatter = LongitudeFormatter()
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.coastlines()
    else:
        if verbose:
            print("Setting region: Lon (%d,%d) Lat (%d,%d)" % (region[0],region[1],region[2],region[3]) )
        ax.set_extent(region, cp.crs.PlateCarree())
        # Adding lines
        draw_parallels(ax, np.arange(60,91,10))
        draw_meridians(ax, np.arange(-180,181,45))    
        ax.coastlines(resolution='50m')
    
    # Generate cyclic data
    cyclic = addcyclicpoint(data, data['lon'])
    try:
        cyclic.attrs['units'] = data.attrs['units']
    except KeyError:
        try:
            cyclic.attrs['units'] = data.attrs['unit']
        except KeyError:
            cyclic.attrs['units'] = ''
    if mode=='None':
        data.plot(ax=ax,**kargs)
    elif mode=='sum':
        data.sum(dim='time').plot(ax=ax,**kargs)
    elif mode=='mean':
        data.mean(dim='time').plot(ax=ax,**kargs)
    ax.set_title(title)

def get_molarweight(data):
    '''
    Extract dictonary of tracers* molarweights from OsloCTM3 file.
    '''
    tracer_names = np.squeeze(data['tracer_name'].data)          # Extract the names of the tracers
    tracer_names = [name.strip() for name in tracer_names]       # Trim whitespace
    tracer_molweight = np.squeeze(data['tracer_molweight'].data) # Extract the molarweight
    return(xr.DataArray(tracer_molweight, coords=[('name', tracer_names)])) # Return DataArray

def get_vmr(data, **karg):
    '''
    Compute volume mixing ratio (VMR) from mass mixing ratio (MMR).
    '''
    tracer = karg.pop('tracer', 'none')
    unit = karg.pop('unit', 'mol/mol')
    Mtracer = karg.pop('molw',  get_molarweight(data).sel(name=tracer)) # g/mol
    mair =  karg.pop('mair', data['AIR']) # g/cm3
    Mair = 28.949 # g/mol
    mtracer = data[tracer] # g/cm3
    
    vmr = mtracer/mair*Mair/Mtracer
    if unit=='ppm':
        vmr = vmr*10**6
    elif unit=='ppb':
        vmr = vmr*10**9
    elif unit=='ppt':
        vmr = vmr*10**12
    else:
        unit = 'mol/mol'    # Fall back to default
    vmr.attrs['units'] = unit
    return(vmr)

    
def get_month_name(mid,**karg):
    '''
    Generate month names. Options: length (int)
    '''
    length = karg.pop('length', 0)
    month_name = {1:'January', 2:'February', 3:'March', 4:'April', 5:'May', 6:'June', 7:'July', 8:'August', 9:'September', 10:'October', 11:'November', 12:'December'}
    if length==0:
        return(month_name[mid])
    else:
        return(month_name[mid][:length])
    

def time_lagged_corr(test_data, truth, **kargs):
    '''
    Compute time lagged correlation and returns correlation cooeficient.
    Takes test_data, truth as arguments. test_data is the one that gets shifted.
    kwargs:
       lag (0) - the time-lag in hours (can be positive or negative)
       v (False) - be a bit more verbose
       pandas (False) - use pandas dataframe object as input instead of numpy arrays
    '''
    lag = kargs.pop('lag',0)
    verbose = kargs.pop('v', False)
    pandas = kargs.pop('pandas', False)
    if pandas:
        corr_coef = truth.corr(test_data.shift(lag))
        if verbose:
            corr_sign = corr_coef*np.sqrt((len(test_data)-np.fabs(lag)-2)/(1-corr_coef**2))
            print "%d %d %1.2f -> %2.2f" % (lag, len(test_data), corr_coef, corr_sign)
        return(corr_coef)
    else:
        if lag >= 0:
            corr_coef = np.ma.corrcoef(np.roll(test_data, lag)[lag:], truth[lag:])
        else:
            corr_coef = np.ma.corrcoef(np.roll(test_data, lag)[:lag], truth[:lag])
        corr_sign = corr_coef[0,1]*np.sqrt((len(test_data)-np.fabs(lag)-2)/(1-corr_coef[0,1]**2))
        if verbose:
            print "%d %d %1.2f -> %2.2f" % (lag, len(test_data), corr_coef[0,1], corr_sign)
        return corr_coef

def read_station_data_noaa(infile, **kargs):
    '''
    Import of NOAA station data and return as pandas Series.
    **kargs utc - timeshift wrt UTC
            start_data - line number in which data block starts
            column - column index in which data is found
            station - standard, southpole, other format
    '''
    import codecs
    # Key words
    cdef int start_data = kargs.pop('start_data', 0)
    cdef double utc = kargs.pop('utc', 0)
    cdef int column = kargs.pop('column', 2)
    station = kargs.pop('station', 'standard')
    # Reading data
    lines = codecs.open(infile, "r", 'utf-8').read().splitlines()
    # Compute the timeshift
    timeshift = dt.timedelta(hours=utc)
    # Data structure
    time = []
    ozone_data_raw = []
    for each in lines[start_data:]:
        # Split the lines
        col = each.split()
        # Seperate ozone data
        ozone_data_raw.append(float(col[column]))
        # Generate timestamp
        if (station=='southpole') | (station=='other'):
            if  int(col[4]) == 24:
                temp_time = dt.datetime.strptime("%s-%s-%s %s" % (col[1], col[2], col[3], '00'), '%Y-%m-%d %H')+dt.timedelta(days=1)
            else:
                temp_time = dt.datetime.strptime("%s-%s-%s %s" % (col[1], col[2], col[3], col[4]), '%Y-%m-%d %H')
        else:
            if int(col[1][:2]) == 24:
                temp_time = dt.datetime.strptime("%s" % (col[0]+" 00:00"), '%Y-%m-%d %H:%M')+dt.timedelta(days=1)
            else:
                temp_time = dt.datetime.strptime("%s" % (col[0]+" "+col[1]), '%Y-%m-%d %H:%M')
        time.append(temp_time-timeshift)
    # Mask the data (negative ozone values)
    ozone_data = np.ma.masked_where(np.array(ozone_data_raw)<0, np.array(ozone_data_raw))
    x_time = np.ma.masked_where(np.array(ozone_data_raw)<0, np.array(time))
    # Output pandas series
    pd_series = pd.Series(ozone_data, index=x_time)
    return pd_series


# Get the read-in routine for EBAS data
def read_station_data_ebas(infile, **kargs):
    '''
    Conversion of NASA aimes files from nilo.no data sets.
    '''
    import nappy as nap # Read and write NASA data files
    # Open file (read header only)
    nas_file = nap.openNAFile(infile)
    # Read actual data
    nas_file.readData()
    # Access the ozone data
    data_raw = np.array(nas_file.getNADict()['V'])
    # Close the nas-file
    nas_file.close()
    # Keywords
    conversion = kargs.pop('conversion', True)
    tracer = kargs.pop('tracer', ('O3',))
    verbose = kargs.pop('v', False)
    if (data_raw.shape[0]-1)/2 < len(tracer):
        for each in tracer:
            if nas_file.getNADict()['NCOM'][-1].find(each) > 0:
                return(read_station_data(infile, tracer=(each,)))
    # Conversion = 1/air_dens(25degC)/mass_fraction
    # 1/2 for O3 (0.5)
    # ca. 0.38 for SO2 
    # ca. 0.81 for NO
    cdef double M_air = 28.949                        # [g/mol] 
    cdef double air_dens_25 = 1.1839                  # [kg/m3]
    mass_fraction = {'O3':3*15.9994/M_air,
                     'SO2':32.065/M_air,# ug(S)/m3
                     'SO4':32.065/M_air,
                     'NO':14.0067/M_air,# ug(N)/m3
                     'NO2':14.0067/M_air}
       
    nas_date = nas_file.getNADict()['DATE']
    start_date = dt.datetime.strptime("%s-0%s-0%s 00:00:00" % (nas_date[0], nas_date[1], nas_date[2]), '%Y-%m-%d %H:%M:%S')
    # Filter the data (FLAG: 0 - valid, >0 - invalid)
    # Add day fraction at stop time to start date 
    x_time_station = np.ma.masked_where(data_raw[2]>0, 
                                        datetime_from_time(start_date, data_raw[0])[0])
    data_dic = {}
    if verbose:
        for each in tracer:
            print(each, data_raw.shape,
                  (nas_file.getNADict()['NCOM'][-1]).split().index(each),
                  (nas_file.getNADict()['NCOM'][-1]).split()[(nas_file.getNADict()['NCOM'][-1]).split().index(each)])
    #i = 1
    for each in tracer:
        # Find column in which tracer appears
        # Split table header and look for it
        # Table header has start/end time while
        # nas_file.getNADict()['V'] only gives one time field
        i = (nas_file.getNADict()['NCOM'][-1]).split().index(each)-1
        if conversion:
            data_dic[each] = np.ma.masked_where(data_raw[i+1]>0, data_raw[i]*1/air_dens_25/mass_fraction[each])
        else:
            data_dic[each] = np.ma.masked_where(data_raw[i+1]>0, data_raw[i])
        #i += 1
    pd_frame = pd.DataFrame(data_dic, index=x_time_station)
    return (pd_frame)
