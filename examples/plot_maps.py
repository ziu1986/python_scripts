import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
#import cartopy.geodesic as cgeo
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
from mytools.station_info import station_location
from mytools.plot_tools import print_all



def plot_stations(ax):
    def plot_frame(ax, start_lon, end_lon, start_lat, end_lat):
        # Plot frame
        linewidth = 15
        ax.plot(np.array((start_lon,end_lon)), np.array((start_lat, start_lat)), color='black', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        ax.plot(np.array((start_lon,end_lon)), np.array((end_lat, end_lat)), color='black', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        for x in np.arange(np.ceil(start_lon),end_lon,2):
            ax.plot((x+linewidth/100.,x+1-linewidth/100.), np.array((start_lat, start_lat)), color='white', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
            ax.plot((x+linewidth/100.,x+1-linewidth/100.), np.array((end_lat, end_lat)), color='white', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        
        ax.plot(np.array((start_lon, start_lon)), np.array((start_lat, end_lat)), color='black', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        ax.plot(np.array((end_lon, end_lon)), np.array((start_lat, end_lat)), color='black', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        for y in np.arange(np.ceil(start_lat),end_lat,2):
            ax.plot((start_lon, start_lon), np.array((y+0.06, y+1-0.04)), color='white', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
            ax.plot((end_lon, end_lon), np.array((y+0.06, y+1-0.04)), color='white', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())

    
    ax.plot(station_location['Jergul'].lon, station_location['Jergul'].lat, fillstyle='left', marker='o', color='orange', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax.plot(station_location['Karasjok'].lon, station_location['Karasjok'].lat, fillstyle='right', marker='o', color='orange', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax.plot(station_location['Svanvik'].lon, station_location['Svanvik'].lat, marker='o', color='blueviolet', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax.plot(station_location['Esrange'].lon, station_location['Esrange'].lat, marker='o', color='blue', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax.plot(station_location['Pallas'].lon, station_location['Pallas'].lat, marker='o', color='black', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax.plot(station_location['Janiskoski'].lon, station_location['Janiskoski'].lat, marker='o', color='grey', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
    text_transform = offset_copy(geodetic_transform, units='dots', x=-25)
    text_transform_2 = offset_copy(geodetic_transform, units='dots', y=40)

    # Add text.
    ax.text(station_location['Jergul'].lon, station_location['Jergul'].lat, u'Jergul',
              verticalalignment='center', horizontalalignment='right', size='xx-large',
              transform=text_transform,
              bbox=dict(facecolor='orange', alpha=0.5, boxstyle='round'))
    ax.text(station_location['Karasjok'].lon, station_location['Karasjok'].lat, u'Karasjok',
              verticalalignment='center', horizontalalignment='right', size='xx-large',
              transform=text_transform_2,
              bbox=dict(facecolor='orange', alpha=0.5, boxstyle='round'))
    ax.text(station_location['Svanvik'].lon, station_location['Svanvik'].lat, u'Svanvik',
              verticalalignment='center', horizontalalignment='right', size='xx-large',
              transform=text_transform,
              bbox=dict(facecolor='blueviolet', alpha=0.5, boxstyle='round'))
    ax.text(station_location['Esrange'].lon, station_location['Esrange'].lat, u'Esrange',
              verticalalignment='center', horizontalalignment='right', size='xx-large', color='white',
              transform=text_transform_2,
              bbox=dict(facecolor='blue', alpha=0.5, boxstyle='round'))
    ax.text(station_location['Pallas'].lon, station_location['Pallas'].lat, u'Pallas',
              verticalalignment='center', horizontalalignment='right', size='xx-large', color='white',
              transform=text_transform,
              bbox=dict(facecolor='black', alpha=0.5, boxstyle='round'))
    ax.text(station_location['Janiskoski'].lon, station_location['Janiskoski'].lat, u'Janiskoski',
              verticalalignment='center', horizontalalignment='right', size='xx-large',
              transform=text_transform,
              bbox=dict(facecolor='grey', alpha=0.5, boxstyle='round'))

    plot_frame(ax, 19.4,31.4,67.6,71.4)
    

       
plt.close("all")

fig1 = plt.figure(1)
fig1.canvas.set_window_title("map")
stamen_terrain = cimgt.Stamen('terrain-background')
ax11 = fig1.add_subplot(1,1,1, projection=stamen_terrain.crs) #ccrs.PlateCarree()
#center = np.array((30.03, 69.45))
#delta = np.array((0.01,0.005))
#print([center[0]-delta[0], center[0]+delta[0], center[1]-delta[1], center[1]+delta[1]])
#ax11.set_extent([center[0]-delta[0], center[0]+delta[0], center[1]-delta[1], center[1]+delta[1]], crs=ccrs.PlateCarree())
#ax11.plot(station_location['Svanvik'].lon, station_location['Svanvik'].lat, marker='o', color='blueviolet', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
ax11.set_extent([19.4,31.4,67.6,71.4], crs=ccrs.Geodetic())#PlateCarree())
# Fennoscadia 19.4,31.4,67.6,71.4
# Middle east 36, 55, 25,42
ax11.add_image(stamen_terrain, 8)

plot_stations(ax11)

plt.show(block=False)
