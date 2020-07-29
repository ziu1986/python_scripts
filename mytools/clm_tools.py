def load_data(src, **karg):
    import glob, os, sys
    import xarray as xr
    '''
    Load CLM history files and select varables.
    Parameters
    ----------
    src : string
        The directory where the hist files are stored.
    Keywords
    --------
    var : list of strings
        A list of valid variable names found in the files.
    Returns
    -------
    data : xarray DataFrame
        Concatenated data.
    '''
    var_list = karg.pop('var', ['GSSHA', 'GSSUN', 'JMX25T', 'Jmx25Z', 'PSNSHA', 'PSNSUN', 'RSSHA','RSSUN', 'VCMX25T', 'Vcmx25Z', 'TOTVEGC', 'TOTVEGN', 'NPP', 'GPP'])
    data_list = []
    for file in sorted(glob.glob(src)):
        print("Loading... %s" % (file))
        data = xr.open_dataset(file)
        data_list.append(data[var_list])
    data = xr.concat(data_list, dim='time')
    return(data)

