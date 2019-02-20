import os, glob # Access environment variables
import numpy as np
import datetime as dt  # Python standard library datetime module
import xarray as xr
import matplotlib.pyplot as plt

src = os.environ['DATA'] # Access environment variable for directory
subdir = '/Shipping/'

def read_shipping_fc(infile, **kargs):
    '''
    Read the fuel consumption of global shipping from ASCII file.
    '''
    verbose = kargs.pop('v', False)
    file = open(infile, 'r')
    lines = file.readlines()
    l_header = kargs.pop('l_header', 3)
    header = lines[l_header][1:].split()
    year = []
    raw_data = []
    if verbose:
        print header

    for line in lines[l_header:]:
        cols = line.split()
        try:
            cols[0]
        except IndexError:
            if verbose:
                print len(cols), "characters in line -> continue"
            continue
        if cols[0] != '%':
            year.append(int(cols[0]))
            for each in cols[1:]:
                raw_data.append(float(each))
    # Close the file
    file.close()
    data = np.array(raw_data).reshape(len(year),3)
    data = np.ma.masked_values(data, -1000)
    #data_dic = dict(zip(header, [year, data[:,0], data[:,1], data[:,2]))
    return(year, data, header)
    
infile = 'fuel_consumption_2007-2015.dat'

year, data, header = read_shipping_fc(src+subdir+infile)


# Plotting
plt.close("all")
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("Fuel_consumption")

ax11 = plt.subplot()
for each in range(0,3):
    ax11.scatter(year,data[:,each],label=header[each+1])

ax11.set_xlabel("Time (years)")
ax11.set_ylabel("Fuel consumption ($10^9$t)")
ax11.legend(loc='best')

#Show it  
plt.show(block=False)
