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
        Standard: Compute mean climatology
        Climatology of daily min/max.
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
        clim_ozone_stderr = data_res.groupby(data_res.index.dayofyear).apply(lambda x: x.mean()/np.sqrt(x.count()))

    return(clim_ozone, clim_ozone_std, clim_ozone_stderr)



