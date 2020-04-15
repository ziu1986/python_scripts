import os, sys, glob
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.station_info import *

def plot_climadiagram(temp, precip, **kwarg):
    bfrost = kwarg.pop('frost', False)
    bdrought = kwarg.pop('drought', False)
    verbose = kwarg.pop('verbose', False)
    temp_err = kwarg.pop('temp_err', None)
    precip_err = kwarg.pop('precip_err', None)

    fig1 = plt.figure(1)
    fig1.subplots_adjust(right=0.9)
    fig1.canvas.set_window_title("climatediagram")
    ax11 = plt.subplot()
    ax12 = ax11.twinx()
    
    precip.plot.bar(ax=ax12, color='blue', alpha=0.5, label='')
    #if temp_err==0:
    #    ax11.plot(ax12.get_xticks(), temp, color='red', label='')
    #else:
    ax11.plot(ax12.get_xticks(), temp, color='red', label='')
    
    if temp_err is not None:
        plot_error_bands(ax11, ax12.get_xticks(), temp, error=temp_err, color='red', ls='None')
        
    if precip_err is not None:
        ax12.errorbar(ax12.get_xticks(), precip, yerr=precip_err, ls='None', color='blue', label='')

    if precip_err is None:
        precip_max = np.ceil(precip.values.max())
        ax11.text(9, precip_max*0.5-5, "$_{max}: %2.1f\,mm$\n$_{min}: %2.1f\,mm$" % (precip.max(), precip.min()), color='blue', size=14)
    else:
        precip_max = np.ceil(precip.values.max()+precip_err.max())
        ax11.text(8.5, precip_max*0.5-5, "$_{max}: (%2.1f\pm %2.1f)\,mm$\n$_{min}: (%2.1f\pm %2.1f)\,mm$" % (precip.max(), precip_err[np.where(precip==precip.max())[0]], precip.min(), precip_err[np.where(precip==precip.min())[0]]), color='blue', size=14)


    temp_max = np.ceil(temp.max())
    temp_min = np.floor(temp.min())
    
    if temp_err is not None:
        temp_max = np.ceil(temp_max+temp_err.max())
        temp_min = np.floor(temp_min-temp_err.max())
        ax11.text(0, precip_max*0.5-5, "$_{max}: (%2.1f\pm %2.1f)^\circ C$\n$_{min}: (%2.1f\pm %2.1f)^\circ C$" % (temp.max(), temp_err[np.where(temp==temp.max())[0]], temp.min(), temp_err[np.where(temp==temp.min())[0]]), color='red', size=14)
    else:
        ax11.text(0, precip_max*0.5-5, "$_{max}: %2.1f^\circ C$\n$_{min}: %2.1f^\circ C$" % (temp_max, temp_min), color='red', size=14)
    
    if verbose:
        print(temp_max, temp_min, precip_max)
        print(np.min((0,temp_min)), np.max((temp_max, precip_max*0.5)))
        
    ax11.set_ylim(np.min((0,temp_min)), np.max((temp_max, precip_max*0.5)))
    ax11.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)])
    ax11.set_ylabel("Temperature ($^\circ\,C$)", color='red')
    ax12.set_ylabel("Precipitation (mm)", color='blue')
    ax12.set_ylim(np.min((0,temp_min*2)), np.max((temp_max*2, precip_max)))

    if bfrost:
        frost = np.where(temp<=0)[0]
        if verbose:
            print(frost+1)
        ax11.scatter(frost, np.zeros_like(frost), marker='*', s=12**2, label='frost months')

    if bdrought:
        drought = np.where(precip<temp*2)[0]
        if verbose:
            print(drought+1, temp[drought])
            ax11.fill_between(drought, 0, temp[drought], color='black')
    if bfrost or bdrought:
        ax11.legend()

# Clean up
plt.close('all')

src_svanvik_clim = os.environ['DATA']+'/astra_data/observations/temp_prec_rad/temp_svanvik*'
src_svanvik_prec_clim = os.environ['DATA']+'/astra_data/observations/temp_prec_rad/precip_svanvik_1992_2020.csv'

data_list = []
try:
    data_svanvik_clim
except NameError:
    for file in sorted(glob.glob(src_svanvik_clim)):
        print(file)
        data_svanvik_clim = pd.read_csv(file)
        data_svanvik_clim.index = pd.date_range("%s" % data_svanvik_clim['Time measured'][0][:-3], "%s" % data_svanvik_clim['Time measured'].iloc[-1][:-3], freq='H')
        data_svanvik_clim = data_svanvik_clim.drop(columns=['Time measured'])
        data_list.append(data_svanvik_clim)
        
    data_svanvik_clim = pd.concat(data_list)
    data_svanvik_precip_clim = pd.read_csv(src_svanvik_prec_clim)
    data_svanvik_precip_clim.index = pd.date_range("%s" % data_svanvik_precip_clim['Time Measured'][0][:-3], "%s" % data_svanvik_precip_clim['Time Measured'].iloc[-1][:-3], freq='D')
    data_svanvik_precip_clim = data_svanvik_precip_clim.drop(columns=['Time Measured'])
            
    data_svanvik_precip_clim.loc[:,'month'] = data_svanvik_precip_clim.index.month.values
    data_svanvik_precip_clim.loc[:,'year'] = data_svanvik_precip_clim.index.year.values

# Plot it
temperature = data_svanvik_clim.iloc[:,0].groupby(data_svanvik_clim.index.month)
precipitation = (data_svanvik_precip_clim[:'2012'].groupby(['year','month']).sum()).groupby(['month'])
plot_climadiagram(temperature.mean().values,
                  precipitation.mean()['RR'],
                  temp_err=temperature.std().values/np.sqrt(data_svanvik_clim.index.year.unique().size),
                  precip_err=precipitation.std()['RR'].values/np.sqrt(data_svanvik_precip_clim.index.year.unique().size),
                  verbose=True, frost=True)
#
# Show it
plt.show(block=False)
