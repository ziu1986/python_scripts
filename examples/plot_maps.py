import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
#import cartopy.geodesic as cgeo
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
from mytools.station_info import station_location

plt.close("all")

fig1 = plt.figure(1)
fig1.canvas.set_window_title("map")
stamen_terrain = cimgt.Stamen('terrain-background')
ax11 = fig1.add_subplot(1,1,1, projection=stamen_terrain.crs) #ccrs.PlateCarree()
center = np.array((30.03, 69.45))
delta = np.array((0.01,0.005))
print([center[0]-delta[0], center[0]+delta[0], center[1]-delta[1], center[1]+delta[1]])
ax11.set_extent([center[0]-delta[0], center[0]+delta[0], center[1]-delta[1], center[1]+delta[1]], crs=ccrs.PlateCarree())
ax11.plot(station_location['Svanvik'].lon, station_location['Svanvik'].lat, marker='o', color='blueviolet', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
#ax11.set_extent([19.4,31.4,67.6,71.4], crs=ccrs.PlateCarree())
# Fennoscadia 19.4,31.4,67.6,71.4
# Middle east 36, 55, 25,42
ax11.add_image(stamen_terrain, 8)

plt.show(block=False)
