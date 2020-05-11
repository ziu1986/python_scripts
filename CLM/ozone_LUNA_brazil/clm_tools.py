def load_data(src):
    import glob, os, sys
    import xarray as xr

    data_list = []
    for file in sorted(glob.glob(src)):
        print("Loading... %s" % (file))
        data = xr.open_dataset(file)
        data_list.append(data[['GSSHA', 'GSSUN', 'JMX25T', 'Jmx25Z', 'PSNSHA', 'PSNSUN', 'RSSHA','RSSUN', 'VCMX25T', 'Vcmx25Z']])
    data = xr.concat(data_list, dim='time')
    return(data)
