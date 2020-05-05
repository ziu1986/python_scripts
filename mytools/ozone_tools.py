# Tools used for ozone and dose analysis
# F.e. in plot_ozone_observations
def load_data(src, **karg):
    '''
    Load ozone station data and round to full hours.
    Parameters
    ----------
    src : string
        Path to the files that shall be concatenated.
    Keyword arguments
    -----------------
    species : string
        Species that shall be extracted from file.
        Standard: "O3"
    type : string
        Type of source. Choices: ebas/barrow
    Returns
    -------
    pandas Timeseries.
    '''
    import os, glob, sys
    from mytools.met_tools import read_station_data_ebas, read_station_data_noaa
    import pandas as pd
    
    species = karg.pop("species", "O3")
    src_type = karg.pop("type", "ebas")
    data = []
    for file in sorted(glob.glob(src)):
        print("Reading file %s" % (file))
        if src_type=="ebas":
            tmp = read_station_data_ebas(file)
            data.append(tmp[species]) 
        elif ((src_type=="Barrow") | (src_type=="barrow")) :
            if int(file[-4:]) < 2003:
                tmp = (read_station_data_noaa(file, utc=-9, start_data=28))
            elif int(file[-4:]) < 2012:
                tmp = (read_station_data_noaa(file, utc=-9, station='other', column=5))
            data.append(tmp)   
    # Concatenate the lists
    print('Concatenating data...')
    data = pd.concat(data)
    # Round to full hours
    data.index = data.index.round("h")
    return(data)


def compute_aot(data, **karg):
    '''
    Compute AOTx and SUMx.
    Parameters
    ----------
    data : pandas Timeseries
        Keyword arguments
    -----------------
    level : int
        Critical threshold. Standard: 40 ppb
    month_start : int
        Number of the month to start the integration.
        Standard: 5 (May)
    month_end : int
        Number of the month to stop the integration (inclusive).
        Standard: 8 (August)
    time_start : int
        Time of the day to start the integration (0-24).
        Standard: 8
    time_end : int
        Time of the day to stop the integration (0-24).
        Standard: 20
    rolling : boolean
        If True, compute SUMx
    Returns
    -------
    aotX : float
        Daily integrated ozone exposion.
    sumX : float
        Seasonal integrated ozone exposion.
    '''
    import numpy as np
    
    threshold = karg.pop('level', 40)
    month_start = karg.pop('month_start', 5)
    month_end = karg.pop('month_end', 8)
    time_start = karg.pop('time_start', 8)
    time_end = karg.pop('time_end', 20)
    rolling = karg.pop('rolling', False)

    selection = data.where(
        (data.index.hour>=time_start)&
        (data.index.hour<=time_end)&
        (data.index.month>=month_start)&
        (data.index.month<=month_end)).dropna()
    
    delta = selection-threshold
    aot = delta.where(delta>0).dropna().resample('1D').sum()
    
    #print(aot)
    if rolling:
        sumX = {}
        print("Applying rolling window 90 days")
        for iyear in aot.index.year.unique():
            sumX[str(iyear)] = (aot.where(aot.index.year==iyear).dropna()).rolling(3*30, center=True).apply(np.sum).shift(-42).dropna()
            #print(sumX[str(iyear)])
    else:
        sumX = aot.groupby(aot.index.year).sum()
    return(sumX)


def compute_climatology(data, **karg):
    '''
    Compute daily ozone climatology from observation.
    Parameters
    ----------
    data : pandas Timeseries
    Keyword arguments
    -----------------
    mode : string
        mean: Compute mean climatology (default).
        min/max: Climatology of daily min/max.
        hourly: Compute hourly climatology.
    Returns
    -------
    clim_ozone : pandas Timeseries
        Daily ozone climatology.
    clim_ozone_std : pandas Timeseries
        Uncertainty on daily ozone climatology.
    clim_ozone_stderr : pandas Timeseries
        Standard error on daily ozone climatology.
    '''
    import numpy as np
    import pandas as pd
    mode = karg.pop("mode", "mean")
    if mode=="mean":
        clim_ozone = data.groupby(data.index.dayofyear).apply(np.nanmean)
        clim_ozone_std = data.groupby(data.index.dayofyear).apply(np.nanstd)
        clim_ozone_stderr = data.groupby(data.index.dayofyear).apply(lambda x: x.mean()/np.sqrt(x.count()))
    elif mode=='max':
        data_res = data.resample("1D").apply(np.nanmax)
        clim_ozone = data_res.groupby(data_res.index.dayofyear).apply(np.nanmean)
        clim_ozone_std = data_res.groupby(data_res.index.dayofyear).apply(np.nanstd)
        clim_ozone_stderr = data_res.groupby(data_res.index.dayofyear).apply(lambda x: x.mean()/np.sqrt(x.count()))
    elif mode=='min':
        data_res = data.resample("1D").apply(np.nanmin)
        clim_ozone = data_res.groupby(data_res.index.dayofyear).apply(np.nanmean)
        clim_ozone_std = data_res.groupby(data_res.index.dayofyear).apply(np.nanstd)
        clim_ozone_stderr = data_res.groupby(data_res.index.dayofyear).apply(lambda x: x.std()/np.sqrt(x.count()))
    elif mode=='hourly':
        tmp = pd.DataFrame({'O3':data.values}, index=data.index)
        #creating the hour, day, month, & day columns
        tmp.loc[:,'hour'] = tmp.index.hour.values
        tmp.loc[:,'day'] = tmp.index.day.values
        tmp.loc[:,'month'] = tmp.index.month.values

        #create groups and calculate the mean of each group
        clim_ozone = tmp.groupby(['month','day','hour']).mean()
        clim_ozone_std = tmp.groupby(['month','day','hour']).std()
        clim_ozone_stderr = tmp.groupby(['month','day','hour']).apply(lambda x: x.std()/np.sqrt(x.count()))

    return(clim_ozone, clim_ozone_std, clim_ozone_stderr)


def flunder(x, **kwarg):
    '''
    Flatten any kind of list of lists, numpy.arrays, and numbers.
    Parameters
    ----------
    x : list of lists etc.
    Keyword arguments
    -----------------
    verbose : bool
    Returns
    -------
    Flat numpy array.
    '''
    import numpy as np
    verbose = kwarg.pop('verbose', False)
    result = []
    for elem in x:
        try:
            for num in elem:
                if verbose:
                    print(num)
                result.append(num)
        except TypeError:
            if verbose:
                print(elem)
            result.append(elem)
    return(np.array(result))

def get_molarweight(data):
    '''
    Extract dictonary of all tracers molarweights from OsloCTM3 file.
    Parameters
    ----------
    data : xarray
        Source OsloCTM3 avgsav file.
    Returns
    -------
    tracer_molweight : a dictionary with corresponding molmasses.
    '''
    tracer_names = np.squeeze(data['tracer_name'].data)          # Extract the names of the tracers
    tracer_names = [name.strip() for name in tracer_names]       # Trim whitespace
    tracer_molweight = np.squeeze(data['tracer_molweight'].data) # Extract the molarweight
    return(xr.DataArray(tracer_molweight, coords=[('name', tracer_names)])) # Return DataArray

def get_vmr(data, **karg):
    '''
    Compute volume mixing ratio (VMR) from mass mixing ratio (MMR).
    Parameters
    ----------
    data : xarray
        Source OsloCTM3 avgsav file.
    Keyword arguments
    -----------------
    tracer : string
        Valid name of a valid tracer.
    unit : string
        Traget unit for the conversion.
        Standard: mol/mol
        Option: ppm, ppb, ppt
    Mtracer : float
        Molecular mass of the tracer.
        Standard: infer from data.
        Option: user defined number.
    mair : float
        Molecular density of air
        Standard: use field 'AIR' in data
        Option: user defined input.
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



def compute_column_density(press, atm_var, **kargs):
    '''
    Compute the column density.
    Parameters
    ----------
    press : xarray
        The atmospheric pressure field in hPa!
    atm_var : xarray
        The atmospheric tracer field in mol/mol.
    Keyword arguments
    -----------------
    unit : sting
        Define the unit for the output.
        Standard: molecules/cm3
        Option: DU
    accumulate : boolean
        Accumulate all levels.
        Standard: True
    Returns
    -------
    atm_new : xarray
        Global field of the tracer.    
    '''
    # Keyword arguments
    unit = kargs.pop('unit', 'molecules/cm2')
    b_accu = kargs.pop('accumulate', True)
    # Conversion factors
    Mdry = 0.0289645    # molec. wt. of dry air, kg/mol
    constant = N_A/(g*Mdry) # molecules*s2/(m*kg)
    cm2 = 1e-4
    DU = 1e-2/2.69e16
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
    if b_accu:
        atm_new = atm_new.sum(dim='lev')
    if unit == 'DU':
        atm_new = atm_new*constant*DU
    else:
        atm_new = atm_new*constant*cm2
    # Set the new units
    atm_new.attrs['units'] = unit
    return atm_new



