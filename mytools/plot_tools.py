# Tools for easier plotting

def get_month_name(mid,**karg):
    '''
    Get the names of the months. 
    Parameters
    ----------
    mid : int
       The cadinal number of the month starting with 1.
    Keyword arguments
    -----------------
    length : int
       Number of letters which will be returned for each month
    Returns
    -------
    month_name : string
    '''
    length = karg.pop('length', 0)
    month_name = {1:'January', 2:'February', 3:'March', 4:'April', 5:'May', 6:'June', 7:'July', 8:'August', 9:'September', 10:'October', 11:'November', 12:'December'}
    if length==0:
        return(month_name[mid])
    else:
        return(month_name[mid][:length])
    
def plot_month_span(ax, **kargs):
    '''
    Indicate the span of a month in a day-of-year type plot.
    Parameters
    ----------
    ax : axis object
       Axis object to be altered.
    Keyword arguments
    -----------------
    year : int
       The year for which the plot is made.
       Disciminates autmatically between leap and non-leap years.
    fill : boolean
       If True axspans are drawn else axvlines are used.
    
    '''
    import numpy as np
    
    year = kargs.pop('year', 2006)
    fill = kargs.pop('fill', True)
    daysinmonth = np.array([calendar.monthrange(year, imonth)[1] for imonth in range(1,13)])
    if fill:
        for i in np.arange(0,12,2):
            ax.axvspan(daysinmonth[:i].sum(), daysinmonth[:i+1].sum(), color='linen')
    else:
        for i in np.arange(0,12):
            ax.axvline(daysinmonth[:i+1].sum(), color='grey', ls='--')
            
def plot_month_name(ax, **kargs):
    '''
    Mark a month by name in a day-of-year type plot.
    Parameters
    ----------
    ax : axis object
       Axis object to be altered.
    Keyword arguments
    -----------------
    year : int
       The year for which the plot is made.
       Disciminates autmatically between leap and non-leap years.
    ypos : float
       Determines the y position at which the text will be drawn.
       Uses first guess from y-axis tick locations but this may fail.
    mlength : int
       Number of letters that shall be shown for each months.
    size : string
       Size of the text fonts. Standard "medium"
    '''
    import numpy as np
    
    year = kargs.pop('year', 2006)
    ypos = kargs.pop('ypos', ax.yaxis.get_ticklocs()[-1]+ax.yaxis.get_tick_space()*0.1)
    mlength = kargs.pop('mlength', 3)
    size = kargs.pop('size', 'medium')
    daysinmonth = np.array([calendar.monthrange(year, imonth)[1] for imonth in range(1,13)])
    xpos = [daysinmonth[:i].sum() for i in np.arange(0,12)]
    for i in range(1,13):
        ax.text(xpos[i-1]+1, ypos, get_month_name(i, length=mlength), size=size)

def set_pressure_axis(ax, **kargs):
    '''
    Adjust the y-axis of a plot to represent pressure levels.
    Parameters
    ----------
    ax : axis object
    Keyword arguments
    -----------------
    limits : (float, float) 
        A tuple of float to indicate the lower and upper limit.
    ticks : list of int 
        A list of int to replace the standard ticks. (20, 50, 100, 200, 500, 1000)
    label : string
        Label text that shall be shown.
    '''
    from met_tools import USstdP
    import matplotlib.pyplot as plt
    
    limits = sorted(kargs.pop("limits", (1000., USstdP(30*kilo)/hecto)))[::-1] # reverse sorted
    ticks = kargs.pop("ticks", [20, 50, 100, 200, 500, 1000])
    label = kargs.pop("label", "Pressure (hPa)")
    ax.set_ylabel(label)
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_yticks(ticks)
    ax.invert_yaxis()
    ax.set_ylim(limits)

def draw_altitude_axis(ax):
    '''
    Add additional y-axis. Altitude computed as reverse function of 
    height->pressure relationship in USstd atmosphere.
    The pressure-axis ax should be in hPa!
    Parameters
    ----------
    ax : axis object
    Returns
    -------
    New axis object
    '''
    from met_tools import USstdP
    import numpy as np
    
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

def print_all(**kargs):
    '''
    Plot all active figures. 
    Make sure all windows have a meaningfull name:
    Use fig.canvas.set_window_title()!
    Keyword arguments
    -----------------
    target : string
        Directory to save plots to. Standard: './plots'
    type : tuple of strings
        The types the plot hsall be saved as.
        Standard: ('pdf', 'svg', 'png')
    '''
    import matplotlib.pyplot as plt
    import os
    
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

def plot_error_bands(ax, x, y, error, **kargs):
    '''
    Function to create a errorband like plot.
    Parameters
    ----------
    ax : Axis object
       The axis where the band shall be drawn
    x : numpy array
       A numpy array of x values.
    y : numpy array
       A numpy array of y values.
    error : numpy array
       A numpy array of corresponding uncertainties in y.
    Keyword arguments
    -----------------
    color : string
       Define the color of the band. Standard: 'blue'
    ls : string
       Line style. Standard: '-'
    '''
    color = kargs.pop('color', 'blue')
    ls = kargs.pop('ls', '-')
    ax.plot(x, y, color=color, ls=ls)
    ax.fill_between(x, y-error, y+error,
                    color=color, alpha=0.5)


def draw_meridians(ax, meridians, **kargs):
    '''
    Draw meridians in cartopy plot.
    Parameters
    ----------
    ax : axis object
    meridians : list of float/int
    Keyword arguments
    -----------------
    polarcap : int
       Position of ploarcap. Standard: 80 deg
    resolution : int
       Number of meridians.
    '''
    import cartopy as cp        # Globe projections
    import numpy as np
    
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


def addcyclicpoint(data, lon):
    '''
    Wrapps cartopy.util.add_cyclic_point() to xarray.
    Works only on 1-level data.
    Parameters
    ----------
    data : xarray DataArray
    lon : list of longitudes
    Returns
    -------
    cyclic_da : xarray DataArray
    '''
    import cartopy.util as ccrs_util  # Add cyclic
    import xarray as xr
    
    cyclic_data, cyclic_lon = ccrs_util.add_cyclic_point(data.values, lon)
    try:
        cyclic_da = xr.DataArray(cyclic_data, coords=[data['time'], data['lat'], cyclic_lon], dims=['time','lat', 'lon'])
    except KeyError:
        cyclic_da = xr.DataArray(cyclic_data, coords=[data['time'], data['latitude'], cyclic_lon], dims=['time','lat', 'lon'])
    
    return cyclic_da


def plot_ozone_drydep(ax, data, **kargs):
    '''
    Plotting routine for global maps of ozone dry deposition.
    Set the map and axis attributes
    Parameters
    ----------
    ax : axis object
    data : xarray DataArray
    Keyword arguments
    -----------------
    v : boolean
       Verbose output.
    mode : string
       Processing of data. Valid choices: "None", "sum", "mean"
    glob : boolean
       Make a global plot if True else a region needs to be defined.
    region : tuple
       Defines a region (lon west, lon east, lat south, lat north)
    title : string
       Set the title of the plot.
    '''

    import cartopy as cp        # Globe projections
    import cartopy.util as ccrs_util  # Add cyclic
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    
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
