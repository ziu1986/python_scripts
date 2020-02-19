import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
#import cartopy.geodesic as cgeo
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy

plt.close("all")

fig1 = plt.figure(1)
fig1.canvas.set_window_title("map")
stamen_terrain = cimgt.Stamen('terrain-background')
ax11 = fig1.add_subplot(1,1,1, projection=stamen_terrain.crs) #ccrs.PlateCarree()
ax11.set_extent([36, 55, 25,42], crs=ccrs.PlateCarree()) #19.4,31.4,67.6,71.4
ax11.add_image(stamen_terrain, 8)

plt.show(block=False)
