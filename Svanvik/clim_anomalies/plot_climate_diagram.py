import os, sys, glob
from turtle import title
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.ozone_tools import *
from mytools.station_info import *

def plot_climate_diagram(temp, precip, **kwarg):
    bfrost = kwarg.pop('frost', False)
    bdrought = kwarg.pop('drought', False)
    verbose = kwarg.pop('verbose', False)
    temp_err = kwarg.pop('temp_err', None)
    precip_err = kwarg.pop('precip_err', None)
    precip_plot_margin = kwarg.pop('pp_margin',0)

    fig1 = plt.figure(1, figsize=(12,8))
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
        ax11.text(8.25, precip_max*0.5-6.75, "$Precip_{max}: %2.1f\,mm$\n$Precip_{min}: %2.1f\,mm$" % (precip.max(), precip.min()), color='blue', size=16)
    else:
        precip_max = np.ceil(precip.values.max()+precip_err.max())
        ax11.text(7.75, precip_max*0.5-9.25, "$Precip_{max}: (%2.1f\pm %2.1f)\,mm$\n$Precip_{min}: (%2.1f\pm %2.1f)\,mm$\n$\Sigma: (%2.1f\pm %2.1f)\,mm$" % (precip.max(), precip_err[np.where(precip==precip.max())[0]], precip.min(), precip_err[np.where(precip==precip.min())[0]], precip.sum(), np.sqrt((precip_err**2).sum())), color='blue', size=14)


    temp_max = np.ceil(temp.max())
    temp_min = np.floor(temp.min())
    
    if temp_err is not None:
        temp_max = np.ceil(temp_max+temp_err.max())
        temp_min = np.floor(temp_min-temp_err.max())
        ax11.text(0, precip_max*0.5-7.5, "$T_{max}: (%2.1f\pm %2.1f)^\circ C$\n$T_{min}: (%2.1f\pm %2.1f)^\circ C$" % (temp.max(), temp_err[np.where(temp==temp.max())[0]], temp.min(), temp_err[np.where(temp==temp.min())[0]]), color='red', size=14)
    else:
        ax11.text(0, precip_max*0.5-6.5, "$T_{max}: %2.1f^\circ C$\n$T_{min}: %2.1f^\circ C$" % (temp_max, temp_min), color='red', size=16)
    
    if verbose:
        print(temp_max, temp_min, precip_max)
        print(np.min((0,temp_min)), np.max((temp_max, precip_max*0.5)))
        
    ax11.set_ylim(np.min((0,temp_min)), np.max((temp_max, (precip_max+precip_plot_margin)*0.5)))
    ax11.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)], size=18)
    ax11.set_ylabel("Temperature ($^\circ\,C$)", color='red')
    ax12.set_ylabel("Precipitation (mm)", color='blue')
    ax12.set_ylim(np.min((0,temp_min*2)), np.max((temp_max*2, (precip_max+precip_plot_margin))))

    if bfrost:
        frost = np.where(temp<=0)[0]
        if verbose:
            print(frost+1)
        ax11.scatter(frost, np.zeros_like(frost), marker='*', s=16**2, label='frost months')

    if bdrought:
        drought = np.where(precip<temp*2)[0]
        if verbose:
            print(drought+1, temp[drought])
            ax11.fill_between(drought, 0, temp[drought], color='black')
    if bfrost or bdrought:
        ax11.legend(loc='lower center', prop={'size': 14})

def make_boxplot(var, ax, **kwarg):
    color = kwarg.pop('color', 'black')
    notch = kwarg.pop('notch', False)

    box_var = var.boxplot(ax=ax, by='month', showmeans=True,
                notch=notch, bootstrap=1000, return_type='dict')
 
    [[item.set_color(color) for item in box_var[key]['medians']] for key in box_var.keys()]
    [[item.set_markerfacecolor(color) for item in box_var[key]['means']] for key in box_var.keys()]
    [[item.set_linewidth(4)for item in box_var[key]['medians']] for key in box_var.keys()]

    return(box_var)
    

def plot_climate_diagram_box(temp, precip, rad, **kwarg):

    notch = kwarg.pop('notch', False)

    fig1 = plt.figure(figsize=(20,6))
    fig1.canvas.set_window_title("climatediagram_box")
    
    ax11 = plt.subplot(131)
    ax12 = plt.subplot(132)
    ax13 = plt.subplot(133)

    box_temp = make_boxplot(temp,ax11,notch=notch, color='r')
    box_precip = make_boxplot(precip,ax12,notch=notch, color='b')
    box_rad = make_boxplot(rad,ax13,notch=notch, color='black')

    
    for ax, ititle in zip(fig1.axes, ["(a)", "(b)", "(c)"]):
        ax.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)], size=16)
        ax.set_title(ititle, size=18) #y=-0.18
        ax.set_xlabel("")

    ax11.set_ylabel("Temperature ($^\circ\,C$)")
    ax12.set_ylabel("Precipitation (mm)")
    ax13.set_ylabel("Gobal radiation ($W\,m^{-2}$)")

    plt.suptitle('')
    plt.tight_layout()

def plot_climate_diagram_ozone_box(ozone, **kwarg):
    notch = kwarg.pop('notch', False)

    fig1 = plt.figure(figsize=(10,6))
    fig1.canvas.set_window_title("climatediagram_ozone_box")
    
    ax11 = plt.subplot()
    
    box_ozone = make_boxplot(ozone,ax11,notch=notch, color='blueviolet')

    ax11.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)], size=16)
    ax11.set_xlabel("")
    ax11.set_title('')

    ax11.set_ylabel("Ozone (ppb)")
    plt.suptitle('')
    plt.tight_layout()

# Clean up
plt.close('all')

src_svanvik_clim = os.environ['DATA']+'/astra/observations/metdata_svanvik/temp_svanvik*'
src_svanvik_prec_clim = os.environ['DATA']+'/astra/observations/metdata_svanvik/precip_svanvik_1992_2020.csv'
src_svanvik_rad_clim = os.environ['DATA']+'/astra/observations/metdata_svanvik/svanvik_glob_rad_2000_2017.csv'
src_svanvik_ozone_clim = os.environ['DATA']+'/astra/observations/ozone/Svanvik/*.nas'

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

    data_svanvik_rad_clim = pd.read_csv(src_svanvik_rad_clim)
    data_svanvik_rad_clim.index = pd.date_range("%s" % data_svanvik_rad_clim['Time measured'][0][:10], "%s" % data_svanvik_rad_clim['Time measured'].iloc[-1][:-3], freq='H')
    data_svanvik_rad_clim = data_svanvik_rad_clim.drop(columns=['Time measured'])

    data_svanvik_rad_clim.loc[:,'month'] = data_svanvik_rad_clim.index.month.values
    data_svanvik_rad_clim.loc[:,'year'] = data_svanvik_rad_clim.index.year.values

    data_svanvik_ozone_clim = pd.DataFrame({'ozone':load_data(src_svanvik_ozone_clim)})

    data_svanvik_ozone_clim.loc[:,'month'] = data_svanvik_ozone_clim.index.month.values
    data_svanvik_ozone_clim.loc[:,'year'] = data_svanvik_ozone_clim.index.year.values



# Plot it
temperature = data_svanvik_clim.iloc[:,0].groupby(data_svanvik_clim.index.month)
precipitation = (data_svanvik_precip_clim[:'2012'].groupby(['year','month']).sum())

plot_climate_diagram(temperature.mean().values,
                  precipitation.groupby(['month']).mean()['RR'],
                  temp_err=temperature.std().values, 
                  #precip_err=precipitation.groupby(['month']).std()['RR'].values, 
                  verbose=True, frost=True, pp_margin=5)
##/np.sqrt(data_svanvik_clim.index.year.unique().size), 
##/np.sqrt(data_svanvik_precip_clim.index.year.unique().size),

tmp_temp = pd.DataFrame({'temp': data_svanvik_clim.iloc[:,0]})
tmp_temp.loc[:,'month'] = tmp_temp.index.month.values

tmp_rad = pd.DataFrame({'q0': data_svanvik_rad_clim.iloc[:,0]})
tmp_rad.loc[:,'month'] = tmp_rad.index.month.values

tmp_ozone = pd.DataFrame({'ozone': data_svanvik_ozone_clim.iloc[:,0]})
tmp_ozone.loc[:,'month'] = tmp_ozone.index.month.values



plot_climate_diagram_box(tmp_temp, precipitation, tmp_rad[:'2012'], notch=True)

plot_climate_diagram_ozone_box(tmp_ozone, notch=True)

# Show it
plt.show(block=False)
