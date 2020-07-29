import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
#import cartopy.geodesic as cgeo
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
from mytools.station_info import station_location
from mytools.plot_tools import print_all

def read_data(src):
    data = pd.read_csv(src, encoding='utf-8')
    return(data)


def plot_stations(fig, extend, data):
    import codecs
    
    def plot_frame(ax, extend):
        start_lon = extend[0]
        end_lon = extend[1]
        start_lat = extend[2]
        end_lat = extend[3]
        # Plot frame
        linewidth = 15
        # Longitudes
        ax.plot(np.array((start_lon,end_lon)), np.array((start_lat, start_lat)), color='black', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        ax.plot(np.array((start_lon,end_lon)), np.array((end_lat, end_lat)), color='black', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        for x in np.arange(np.ceil(start_lon),end_lon,2):
            ax.plot((x+linewidth/100.,x+1-linewidth/100.), np.array((start_lat, start_lat)), color='white', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
            ax.plot((x+linewidth/100.,x+1-linewidth/100.), np.array((end_lat, end_lat)), color='white', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        # Latitudes
        ax.plot(np.array((start_lon, start_lon)), np.array((start_lat, end_lat)), color='black', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        ax.plot(np.array((end_lon, end_lon)), np.array((start_lat, end_lat)), color='black', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
        for y in np.arange(np.ceil(start_lat),end_lat,2):
            ax.plot((start_lon, start_lon), np.array((y+0.06, y+1-0.04)), color='white', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())
            ax.plot((end_lon, end_lon), np.array((y+0.06, y+1-0.04)), color='white', ls='-', linewidth=linewidth, transform=ccrs.Geodetic())

    # Plot background
    stamen_terrain = cimgt.Stamen('terrain-background')
    ax = fig.add_subplot(1,1,1, projection=stamen_terrain.crs) #ccrs.PlateCarree()
    ax.set_extent(extend, crs=ccrs.Geodetic())#PlateCarree())
    ax.add_image(stamen_terrain, 8)

    geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
    text_transform = offset_copy(geodetic_transform, units='dots', x=-15, y=20)
    text_transform_2 = offset_copy(geodetic_transform, units='dots', x=15, y=20)
    text_transform_3 = offset_copy(geodetic_transform, units='dots', x=-15, y=-20)

    for index, station in data.iterrows():
        # Decode the text into unicode
        name = station.Name #codecs.decode(station.Name)
        # Plot station
        ax.plot(station.lon, station.lat, marker='o', ls='None', color='None', markersize=12, alpha=0.7, transform=ccrs.Geodetic(), label="%2.0f %s" % (index+1, name))
        ax.plot(station.lon, station.lat, marker='o', ls='None', color='black', markersize=12, alpha=0.7, transform=ccrs.Geodetic(), label="_")
        if station.Name=="Vindeln":
            tt = text_transform_2
        elif station.Name=="Rosinedal":
            tt = text_transform_3
        else:
            tt = text_transform
            
        ax.text(station.lon, station.lat, index+1, #unicode(station.Name),
                verticalalignment='center', horizontalalignment='right', size='xx-large', color='white',
                transform=tt,
                bbox=dict(facecolor='black', alpha=0.5, boxstyle='round'))
        
    legend = ax.legend(bbox_to_anchor=(0.95, 0.95))

    #ax.plot(station_location['Jergul'].lon, station_location['Jergul'].lat, fillstyle='left', marker='o', color='orange', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    #ax.plot(station_location['Karasjok'].lon, station_location['Karasjok'].lat, fillstyle='right', marker='o', color='orange', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    #ax.plot(station_location['Svanvik'].lon, station_location['Svanvik'].lat, marker='o', color='blueviolet', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    #ax.plot(station_location['Esrange'].lon, station_location['Esrange'].lat, marker='o', color='blue', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    #ax.plot(station_location['Pallas'].lon, station_location['Pallas'].lat, marker='o', color='black', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    #ax.plot(station_location['Janiskoski'].lon, station_location['Janiskoski'].lat, marker='o', color='grey', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    
    # Add text.
    #ax.text(station_location['Jergul'].lon, station_location['Jergul'].lat, u'Jergul',
    #          verticalalignment='center', horizontalalignment='right', size='xx-large',
    #          transform=text_transform,
    #          bbox=dict(facecolor='orange', alpha=0.5, boxstyle='round'))
    #ax.text(station_location['Karasjok'].lon, station_location['Karasjok'].lat, u'Karasjok',
    #          verticalalignment='center', horizontalalignment='right', size='xx-large',
    #          transform=text_transform_2,
    #          bbox=dict(facecolor='orange', alpha=0.5, boxstyle='round'))
    #ax.text(station_location['Svanvik'].lon, station_location['Svanvik'].lat, u'Svanvik',
    #          verticalalignment='center', horizontalalignment='right', size='xx-large',
    #          transform=text_transform,
    #          bbox=dict(facecolor='blueviolet', alpha=0.5, boxstyle='round'))
    #ax.text(station_location['Esrange'].lon, station_location['Esrange'].lat, u'Esrange',
    #          verticalalignment='center', horizontalalignment='right', size='xx-large', color='white',
    #          transform=text_transform_2,
    #          bbox=dict(facecolor='blue', alpha=0.5, boxstyle='round'))
    #ax.text(station_location['Pallas'].lon, station_location['Pallas'].lat, u'Pallas',
    #          verticalalignment='center', horizontalalignment='right', size='xx-large', color='white',
    #          transform=text_transform,
    #          bbox=dict(facecolor='black', alpha=0.5, boxstyle='round'))
    #ax.text(station_location['Janiskoski'].lon, station_location['Janiskoski'].lat, u'Janiskoski',
    #          verticalalignment='center', horizontalalignment='right', size='xx-large',
    #          transform=text_transform,
    #          bbox=dict(facecolor='grey', alpha=0.5, boxstyle='round'))

    # Plot the frame
    #plot_frame(ax, extend)
    

       
plt.close("all")

data = read_data("stations_report.csv")
fig1 = plt.figure(1)
fig1.canvas.set_window_title("map")

# Fennoscadia 19.4,31.4,67.6,71.4
# Middle east 36, 55, 25,42
plot_stations(fig1, [6.4,31.4,56.6,71.4], data)

plt.show(block=False)
