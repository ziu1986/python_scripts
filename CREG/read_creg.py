import os, sys, glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl

#plt.close('all')

# Read data
src = os.environ['DATA']+'/CREG12.L75-REF08_y2000m02d19.5d_gridU.nc'
src_bath = os.environ['DATA']+'/bathym.nc'

data = xr.open_dataset(src)
data_bath = xr.open_dataset(src_bath)

velo = data['vozocrtx'].squeeze()

# Plot data
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(131)
ax12 = plt.subplot(132)
ax13 = plt.subplot(133)

velo.min(dim='depthu').plot(ax=ax11)
data_bath['mbathy'].plot(ax=ax12)
# Select level based on bathometry and collaps to 2 dim
velo.where(data_bath['mbathy']).sum(dim='depthu').plot(ax=ax13, vmax=20, vmin=-20, cmap=plt.cm.coolwarm)

# Show it
plt.show(block=False)


