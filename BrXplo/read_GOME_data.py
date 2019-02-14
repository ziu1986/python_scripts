import numpy as np
import datetime as dt  # Python standard library datetime module
import xarray as xr

def read_gome_data(infile, **kargs):
    '''
    Read GOME satellite data from ASCII files.
    Returns xarray DataArray with additional meta data.
    '''
    verbose = kargs.pop('v', False)
    lines = open(infile, "r").readlines()
    l_header = kargs.pop('l_header', 23)
    # Meta data
    start_date = lines[3].find('Dates : ')+1
    date = lines[3].split()[start_date]
    if verbose:
        print lines[3], date, len(date)
    x_time = dt.datetime.strptime("%s20%s" % (date[:-2], date[-2:]), '%d.%m.%Y')
        
    var_name = lines[4].split()[1:]
    var_name = "%s %s %s" % (var_name[0], var_name[1], var_name[2])
    scaling = float(lines[l_header-5].split()[-1])
    errorflag = float(lines[l_header-4].split()[-1])
    start_lat = float(lines[l_header-3].split()[4][:-1])
    step_lat = float(lines[l_header-3].split()[8][:-1])
    n_lat = int(lines[l_header-3].split()[15][:-1])
    start_lon = float(lines[l_header-2].split()[4][:-1])
    step_lon = float(lines[l_header-2].split()[8][:-1])
    n_lon = int(lines[l_header-2].split()[15][:-1])
    lat = np.arange(start_lat, start_lat+n_lat*step_lat, step_lat)
    lon = np.arange(start_lon, start_lon+n_lon*step_lon, step_lon)
    data_raw = []
    if verbose:
        print var_name, scaling, errorflag
        print start_lat, step_lat, n_lat
        print start_lon, step_lon, n_lon, lat, lon
    for line in lines[l_header:]:
        cols = line.split()
        try:
            cols[0]
        except IndexError:
            if verbose:
                print len(cols), "characters in line -> continue"
            continue
        if cols[0] != '*':
            for each in cols:
                data_raw.append(float(each))
    data_raw = np.expand_dims(np.array(data_raw).reshape(n_lat, n_lon), axis=0)
    # Mask the flagged data
    data = np.ma.masked_where(data_raw==errorflag, data_raw)
    # Output dictonary
    data_dic = xr.DataArray(data, name='BrO', coords={'time':[x_time,],'lat':lat, 'lon':lon},
                            dims=['time','lat','lon'], 
                            attrs={'units':'molecules/cm2', 'scaling':scaling, 'long_name':var_name})
    if verbose:
        print data_dic
    return data_dic
