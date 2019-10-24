import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *
from mytools.station_info import station_location

def read_date(input_files):
    data_list = []
    print("Reading")
    for file_name in sorted(glob.glob(input_files)):
        print(file_name)
        data = xr.open_dataset(file_name)
        data_list.append(data)
        

    # Concat list
    data = xr.concat(data_list, dim='time')
    return(data)

def plot_sumxi(fig, ax, data, color):
    year = np.unique(data.time.dt.year)[0]
    data['AOT30_8-20'].plot(ax=ax, color=color, ls='-', label='AOT30_8-20')
    data['AOT40_8-20'].plot(ax=ax, color=color, ls='--', label='AOT40_8-20')
    data['AOT30_1-23'].plot(ax=ax, color=color, ls='-.', label='AOT30_1-23')
    data['AOT40_1-23'].plot(ax=ax, color=color, ls=':', label='AOT40_1-23')
    ax.set_ylim(0,20000)
    ax.set_ylabel("SUMx$_i$ (ppb h)")

    ax.axhline(3000, color='black', ls=':')

    text_y = ax.get_ylim()[1]
    ax.axvline(dt.date(day=11, month=5, year=year), color='black', ls=':')
    ax.axvspan(dt.date(day=10, month=5, year=year), dt.date(day=12, month=5, year=year), color='black', alpha=0.25)
    ax.text(dt.date(day=10, month=5, year=year), text_y-1000, "Extrapol. SGS 2100")

    ax.axvline(dt.date(day=1, month=7, year=year), color='orange', ls=':')
    ax.text(dt.date(day=1, month=7, year=year), text_y-1000, "SGS present", color='orange')
    ax.legend(ncol=1, bbox_to_anchor=(3.15, -2.05), loc='lower right', borderaxespad=0.)
    ax.set_title(year)

def plot_aotxi(fig, ax, data, color):
    year = np.unique(data.time.dt.year)[0]
    data['AOT30_8-20'].plot(ax=ax, color=color, ls='-', label='AOT30_8-20')
    data['AOT40_8-20'].plot(ax=ax, color=color, ls='--', label='AOT40_8-20')
    data['AOT30_1-23'].plot(ax=ax, color=color, ls='-.', label='AOT30_1-23')
    data['AOT40_1-23'].plot(ax=ax, color=color, ls=':', label='AOT40_1-23')
    #ax.set_ylim(0,15000)
    ax.set_ylabel("AOTx$_i$ (ppb h)")

    text_y = ax11.get_ylim()[1]
    ax.axvline(dt.date(day=11, month=5, year=year), color='black', ls=':')
    ax.axvspan(dt.date(day=10, month=5, year=year), dt.date(day=12, month=5, year=year), color='black', alpha=0.25)
    ax.text(dt.date(day=10, month=5, year=year), text_y-1000, "Extrapol. SGS 2100")

    ax.axvline(dt.date(day=1, month=7, year=year), color='orange', ls=':')
    ax.text(dt.date(day=1, month=7, year=year), text_y-1000, "SGS present", color='orange')
    #ax.legend(ncol=2)
    ax.set_title(year)

    
def plot_maps(fig, data, metric, **kwarg):
    maximum = kwarg.pop('max', 20000)
    minimum = kwarg.pop('min', 3000)
    stepwidth = kwarg.pop('step', 1000)
    year = np.unique(data.time.dt.year)[0]
    ax1 = plt.subplot(131, projection=ccrs.PlateCarree())
    ax2 = plt.subplot(132, projection=ccrs.PlateCarree())
    ax3 = plt.subplot(133, projection=ccrs.PlateCarree())
    levels = np.arange(minimum,maximum+10,stepwidth)
    levels_diff = np.arange(-2.05,2.1,0.1)

    data.sel(time=np.datetime64('%d-07-01' % (year)))[metric].plot.contourf(ax=ax1, transform=ccrs.PlateCarree(), levels=levels, cmap=plt.cm.hot_r)
    data.sel(time=np.datetime64('%d-05-11' % (year)))[metric].plot.contourf(ax=ax2, transform=ccrs.PlateCarree(), levels=levels, cmap=plt.cm.hot_r)
    (data.sel(time=np.datetime64('%d-05-11' % (year)))[metric]/
     data.sel(time=np.datetime64('%d-07-01' % (year)))[metric]-1).plot.contourf(ax=ax3, transform=ccrs.PlateCarree(), levels=levels_diff, cmap=plt.cm.seismic)

    for ax in fig.axes[:3]:
        ax.set_extent([0,40,56,81], crs=ccrs.PlateCarree())
        ax.set_aspect('auto')
        ax.add_feature(crs.feature.OCEAN, zorder=100, edgecolor='k')
        # Adding lines
        #draw_parallels(ax, np.arange(60,91,5))
        #draw_meridians(ax, np.arange(0,40,5))   
        ax.coastlines()

def rolling(data):
    ag_data = {}
    for each in data:
        ag_data[each] = data[each].rolling(time=30*3, center=True).reduce(np.sum).shift(time=-42)
    return(ag_data)

def select(ag_data, lat,lon):
    ag_data_sel = {}
    for each in ag_data:
        ag_data_sel[each] = ag_data[each].sel(latitude=lat,longitude=lon, method='nearest')
    return(ag_data_sel)

#def main():
    #Data source
try:
    data
except NameError:
    data = {}
    for year in range(2003,2013):
        nc_src = os.environ['DATA']+"/astra_data/ECMWF/MACC_reanalysis/netcdf/aot/*%s0[5-9]*.nc" % str(year).zfill(2)
        data['%d' % year] = read_date(nc_src)

       
    ag_data = rolling(data)
    for istation in ('Esrange', 'Pallas', 'Jergul'):
        ag_data_sel[istation] = select(ag_data, lat=station_location[istation].lat, lon=station_location[istation].lon)
        data_sel[istation] = select(data, lat=station_location[istation].lat, lon=station_location[istation].lon)
# Clean up
plt.close('all')
# Plot it
#xShift = (dt.date(2000, 1, 1), dt.date(2001, 1, 1))
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("sum40_int3month_change_greeningseason")
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)

plot_sumxi(fig1, ax11, ag_data_sel['Jergul']['2005'], 'orange')
plot_sumxi(fig1, ax11, ag_data_sel['Pallas']['2005'], 'black')
plot_sumxi(fig1, ax11, ag_data_sel['Esrange']['2005'], 'blue')
plot_sumxi(fig1, ax12, ag_data_sel['Jergul']['2006'], 'orange')
plot_sumxi(fig1, ax12, ag_data_sel['Pallas']['2006'], 'black')
plot_sumxi(fig1, ax12, ag_data_sel['Esrange']['2006'], 'blue')


fig7 = plt.figure(7, figsize=(16,9))
fig7.canvas.set_window_title("aot40_int3month_change_greeningseason")
ax71 = plt.subplot(211)
ax72 = plt.subplot(212)

plot_aotxi(fig7, ax71, data_sel['Esrange']['2005'], 'blue')
plot_aotxi(fig7, ax72, data_sel['Esrange']['2006'], 'red')


fig2 =  plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("sum40_int3month_change_greeningseason_2003-2012")

for i in range(len(data)):
    ax = plt.subplot(3,4,i+1)
    plot_sumxi(fig2, ax, ag_data_sel['%d' % (2003+i)], 'blue')

for ax in fig2.axes[1:]:
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.get_legend().remove()

for ax in fig2.axes[:-2]:
    ax.set_xticklabels("")
    ax.set_xlabel("")

for ax in fig2.axes[1:4]:
    ax.set_yticklabels("")

for ax in fig2.axes[5:8]:
    ax.set_yticklabels("")
fig2.axes[-1].set_yticklabels("")
#ax.legend(('-','--','-.',':'),('AOT30_8-20','AOT40_8-20','AOT30_1-23','AOT40_1-23'),bbox_to_anchor=(1.05, 1), loc='lower right', borderaxespad=0.)


    
#fig2.axes[1].set_ylabel(ax.get_ylabel(), y=-10)
import cartopy as crs
import cartopy.crs as ccrs

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("AOT40_8-20_int3month_change_greeningseason_map_2012")
plot_maps(fig3, ag_data['2012'], 'AOT40_8-20')

fig4 = plt.figure(4, figsize=(16,9))
fig4.canvas.set_window_title("AOT30_8-20_int3month_change_greeningseason_map_2012")
plot_maps(fig4, ag_data['2012'], 'AOT30_8-20', max=30000)

fig5 = plt.figure(5, figsize=(16,9))
fig5.canvas.set_window_title("AOT40_1-23_int3month_change_greeningseason_map_2012")
plot_maps(fig5, ag_data['2012'], 'AOT40_1-23')

fig6 = plt.figure(6, figsize=(16,9))
fig6.canvas.set_window_title("AOT30_1-23_int3month_change_greeningseason_map_2012")
plot_maps(fig6, ag_data['2012'], 'AOT30_1-23', max=30000)


# Show it
plt.show(block=False)

#if __name__ == "__main__":
#    main()
