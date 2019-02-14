import numpy as np
#import pandas as pd
#import xarray as xr
import matplotlib.pyplot as plt

def local_solar_time(time, longitude):
    '''
    Longitude > 0 eastward, < 0 westward
    '''
    return np.array([time + lon/15. for lon in longitude]) #

# Generate test data
time = np.arange(-1,23)
lon = np.arange(-180, 181, 2.5)

lst = local_solar_time(time, lon)

x, y = np.meshgrid(time, lon)
x1, y1 = np.where((lst<11) & (lst>10))
# Plot it
plt.close('all')
fig1 = plt.figure(1)
cm11 = plt.pcolormesh(x,y,lst)
plt.scatter(time[y1], lon[x1], marker='x', color='black')
ax11 = plt.gca()
ax11.set_xlabel('Time (hours)')
#ax11.set_xlim(0,24)
ax11.set_xticks(np.arange(24))
ax11.set_ylabel('Longitude (deg)')
#ax11.set_ylim(-180,180)
ax11.set_yticks(np.arange(-180,190,10))
cbar = plt.colorbar(cm11)
cbar.set_label("Local solar time")

# Show it
#plt.show(block=False)
