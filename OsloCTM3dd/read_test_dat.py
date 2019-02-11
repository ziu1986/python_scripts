import os, glob # Access environment variables
import numpy as np
import datetime as dt  # Python standard library datetime module
import xarray as xr
import matplotlib.pyplot as plt

src = os.environ['DATA'] # Access environment variable for directory
infile = '/VGstO3_test.dat'

file = open(src+infile, 'r')
lines = file.readlines()

raw_data = []
for line in lines:
    cols = line.split()
    for each in cols:
        raw_data.append(float(each))

data = np.array(raw_data).reshape(80,160)
lat = np.linspace(-90,90,80)
lon = np.linspace(0,360,160)
# Plot it
fig1 = plt.figure(1)
ax11 = plt.subplot()

ax11.contourf(data)

# Show it
plt.show(block=False)
