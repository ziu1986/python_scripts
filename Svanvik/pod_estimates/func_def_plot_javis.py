import os, glob, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt   # Plotting data
#from cmcrameri import cm
import datetime as dt
import copy                       # for deep copy of object
from scipy.constants import *     # Get physics constants
from javis_model.javis_model import *
from mytools.plot_tools import *
from mytools.ozone_tools import VPD

def import_data(src):
    data_list = []
    # Read data
    for file in sorted(glob.glob(src)):
        print(file)
        data_svanvik = pd.read_csv(file)
        data_svanvik.index = pd.date_range("%s" % data_svanvik['Time measured'][0][:-3], "%s" % data_svanvik['Time measured'].iloc[-1][:-3], freq='H')
        data_svanvik = data_svanvik.drop(columns=["Time measured"])

        data_list.append(data_svanvik)
        
    data_svanvik = pd.concat(data_list)

    return(data_svanvik)

def gs_noon(par):
     return(par.where((par.index.month>=5)&(par.index.month<9)&
               (par.index.hour>=11)&(par.index.hour<=13)).dropna())

def gs_morning(par):
     return(par.where((par.index.month>=5)&(par.index.month<9)&
               (par.index.hour>=6)&(par.index.hour<=8)).dropna())

def get_f_function(javis_model, temp, vpd, rad):
    f_temp = temp.iloc[:,0].apply(lambda x: javis_model.f_temp(x))
    f_vpd = vpd.apply(lambda x: javis_model.f_vpd(x))
    f_light = rad.iloc[:,0].apply(lambda x: javis_model.f_light(x))

    max_f = (np.maximum(javis_model.f_min, f_temp*f_vpd))
    prod_f = (f_light*np.maximum(javis_model.f_min, f_temp*f_vpd))
    
    return(f_temp, f_vpd, f_light, max_f, prod_f)

def get_stats(prod_f, **karg):
    var_type = karg.pop('type','noon')
    stats = karg.pop('stats', 'var')
    
    if var_type == 'noon':
        selection = gs_noon(prod_f)
    else:
        selection = gs_morning(prod_f)

    if stats == 'var':
        stats_result =  selection.var()
    elif stats == 'mean':
        stats_result =  selection.mean()
    elif stats == 'std':
        stats_result =  selection.std()
    elif stats == 'stderr':
        stats_result = selection.mean()/np.sqrt(selection.dropna().size)
        
    print(var_type, stats, stats_result)
    return(stats_result)
   
        


def plot_f_functions(javis_model, fig_i, **karg):
    start_date = karg.pop('start', '2019-05')
    end_date = karg.pop('end', '2019-08')

    # Compute f_functions
    f_temp, f_vpd, f_light, max_f, prod_f = get_f_function(javis_model)
    
    # Plot it
    fig = plt.figure(fig_i, figsize=(10,12))
    fig.canvas.set_window_title("javis_funcs_%s" % javis_model.name)
    plt.subplots_adjust(right=0.9)

    ax1 = plt.subplot(411)
    ax2 = plt.subplot(412)
    ax3 = plt.subplot(413)
    ax4 = plt.subplot(414)

    f_temp[start_date:end_date].plot(ax=ax1)
    f_vpd[start_date:end_date].plot(ax=ax2)
    f_light[start_date:end_date].plot(ax=ax3)
    
    #max_f[start_date:end_date].plot(ax=ax3, color='red')
    (prod_f[start_date:end_date]*javis_model.gmax).plot(ax=ax4, color='darkgrey')

    ax1.set_ylabel("f_temp")
    ax2.set_ylabel("f_vpd")
    ax3.set_ylabel("f_light")
    ax4.set_ylabel("$g_{sto}$ $(mmol\,m^{-2}\,s^{-1})$")
    ax4.set_xlabel("Time (months)")

    for ax in fig.axes[:-1]:
        ax.set_ylim(0,1)

    # Additional axis
    ax11 = ax1.twinx()
    (data_temp.iloc[:,0][start_date:end_date]).plot(ax=ax11, color='red', alpha=0.5)
    ax11.set_ylabel('$T$ $(^\circ C)$', color='red')
    ax11.grid(b=False)
    ax21 = ax2.twinx()
    vpd[start_date:end_date].plot(ax=ax21, color='blue', alpha=0.5)
    ax21.set_ylabel('$VPD$ $(Pa)$', color='blue')
    ax21.grid(b=False)

    
def plot_optimal(results, **karg):
    stats = karg.pop('stats', 'var')
    err = karg.pop('err', None)
    
    fig = plt.figure(1, figsize=(10,12))
    fig.canvas.set_window_title("javis_func_opt_%s" % stats)
    for i, iname, ilabel in zip((0,1,2), ('coniferous', 'deciduous', 'grassland'), ('(a)', '(b)', '(c)')):
        ax = plt.subplot(3,1,i+1)
        ax.set_title(ilabel)
        
        noon_results = results[::2][9*i:9*(i+1)]
        morning_results = results[1::2][9*i:9*(i+1)]

        if err:
            noon_err = err[::2][9*i:9*(i+1)]
            morning_err = err[1::2][9*i:9*(i+1)]
        else:
            noon_err = None
            morning_err = None
        
        ax.errorbar(np.arange(len(noon_results))-0.1, noon_results, yerr=noon_err, ls='None', marker='o', label='noon', zorder=2)
        ax.errorbar(np.arange(len(morning_results))+0.1, morning_results, yerr=morning_err, ls='None', marker='o', fillstyle='none', label='morning', zorder=2)
    for ax in fig.axes:
        if stats == 'var':
            ax.set_ylim(0,0.11)
            ax.set_ylabel('Variance')
        else:
            ax.set_ylim(0,1.1)
            ax.set_ylabel('$\\frac{\\left<g_{sto}\\right>}{g_{max}}$')

        ax.set_xticks(np.arange(9))
        ax.set_xlim(-0.2,8.2)
        ax.axvspan(-0.2, 2.5, color='linen', zorder=1)
        ax.axvspan(5.5, 8.2, color='linen', zorder=1)
        ax.set_xticklabels((r'subarctic', r'PPFD0.8', r'PPFD1.2',
                            r'cold', r'PPFD0.8', r'PPFD1.2',
                            r'MM', r'PPFD0.8', r'PPFD1.2'))
        ax.legend(loc='lower left')

    ax.set_xlabel('Categories')
    

def plot_f_light(javis_model, alpha_variation, **karg):
    fig = plt.figure(figsize=(10,8))
    fig.canvas.set_window_title("javis_funcs_f_light_%s" % javis_model.name)
    ax1 = plt.subplot()
    plt.plot(javis_model.f_light(np.arange(1000)), label="$\alpha=%1.4f$" % javis_model.alpha_light)
    for each in alpha_variation:
        tmp_model = copy.deepcopy(javis_model)
        tmp_model.alpha_light = each
    plt.plot(tmp_model.f_light(np.arange(1000)), label="$\alpha=%1.4f$" % tmp_model.alpha_light)

    ax1.set_xlabel("PPFD $(W\,m^{-2}\,s^{-1})$")
    ax1.set_ylabel("$f_{light}$")

    ax1.axhline(0.5, ls=':', color='red')

def plot_temperature_histogram(temperature, **karg):
    javis_model = karg.pop('javis_model', None)
    
    fig = plt.figure(figsize=(10,6))
    ax1 = plt.subplot()
    fig.canvas.set_window_title("javis_funcs_temp_hist")
    
    temp_summer_all = temperature.where((temperature.index.month>=5)&(temperature.index.month<9))
    temp_summer_90 = temperature.where((temperature.index.month>=5)&(temperature.index.month<9))[:'2000']
    temp_summer_00 = temperature.where((temperature.index.month>=5)&(temperature.index.month<9))['2001':'2010']
    temp_summer_10 = temperature.where((temperature.index.month>=5)&(temperature.index.month<9))['2011':]
    if not javis_model:
        temperature.hist(ax=ax1, bins=np.arange(-40, 40), label='all years')
        print("Mean all years summer: %2.3f\\n Mean 90s summer: %2.3f\\n Mean 00s summer: %2.3f" % (temp_summer_all.mean(), temp_summer_90.mean(), temp_summer_00.mean()))
        
        temp_summer_all.dropna().hist(ax=ax1, bins=np.arange(-40, 40), histtype='stepfilled', hatch='//', color='red', label='May-Sep')
        temp_summer_90.dropna().hist(ax=ax1, bins=np.arange(-40, 40), histtype='stepfilled', hatch='\\', color='blue', label='May-Sep (1992-2000)')
        temp_summer_00.dropna().hist(ax=ax1, bins=np.arange(-40, 40), histtype='step', hatch='\\', color='black', label='May-Sep (2001-2010)')
        ax1.set_ylabel("Counts")
    else:
        fig.canvas.set_window_title("javis_funcs_temp_hist_%s" % javis_model[0].name)
        plt.subplots_adjust(right=0.9)
        temp_summer_90.dropna().hist(ax=ax1, density=True, bins=np.arange(-4, 40), color='darkgrey', label='May-Sep (1992-2000)')
        temp_summer_10.dropna().hist(ax=ax1, density=True, bins=np.arange(-4, 40), histtype='step', hatch='\\', color='black', linewidth=2, label='May-Sep (2011-2019)')
        ax11 = ax1.twinx()
        for each in javis_model:
            ax11.plot(np.arange(-4, 40), each.f_temp(np.arange(-4,40)), color='blue', linewidth=2, label='%s' % each.name.replace('_', ' '))
        
        ax11.grid(b=False)
        ax11.set_ylim(0,1.4)
        ax11.legend(fontsize='x-large')
        ax11.tick_params(axis='y', colors='blue')
        ax11.set_ylabel('$f_{T}$', color='blue')
        ax1.set_ylabel("PDF")
        ax1.set_ylim(0,0.09)
        
        
    ax1.set_xlabel("$T$ $(^\circ C)$")
    ax1.legend(fontsize='x-large', loc='upper left')


def plot_histogram(variable, **karg):
    var_mod = karg.pop('var', 'temperature')
    javis_model = karg.pop('javis_model', None)
    
    fig = plt.figure(figsize=(10,6))
    ax1 = plt.subplot()
    fig.canvas.set_window_title("javis_funcs_%s_hist" % (var_mod))
    
    variable_summer_all = variable.where((variable.index.month>=5)&(variable.index.month<9))
    variable_summer_90 = variable.where((variable.index.month>=5)&(variable.index.month<9))[:'2000']
    variable_summer_00 = variable.where((variable.index.month>=5)&(variable.index.month<9))['2001':'2010']
    variable_summer_10 = variable.where((variable.index.month>=5)&(variable.index.month<9))['2011':]

    if var_mod=='temperature':
        var_range = np.arange(-4, 40)
    elif var_mod=='light':
        var_range = np.arange(0, 1000, 10)
    elif var_mod=='vpd':
        var_range = np.arange(0, 5, 0.1)
        
    fig.canvas.set_window_title("javis_funcs_%s_hist_%s" % (var_mod, javis_model[0].name))
    plt.subplots_adjust(right=0.9)
    variable_summer_90.dropna().hist(ax=ax1, density=True, bins=var_range, color='darkgrey', label='May-Sep (1992-2000)')
    variable_summer_10.dropna().hist(ax=ax1, density=True, bins=var_range, histtype='step', hatch='\\', color='black', linewidth=2, label='May-Sep (2011-2019)')
    ax11 = ax1.twinx()
    for each in javis_model:
        if var_mod=='temperature':
            ax11.plot(var_range, each.f_temp(var_range), color='blue', linewidth=2, label='%s' % each.name.replace('_', ' '))
            
        elif var_mod=='light':
            ax11.plot(var_range, each.f_light(var_range), color='blue', linewidth=2, label='%s' % each.name.replace('_', ' '))
        
        elif var_mod=='vpd':
            ax11.plot(var_range, each.f_vpd(var_range), color='blue', linewidth=2, label='%s' % each.name.replace('_', ' '))

    if var_mod=='temperature':
        ax11.set_ylabel('$f_{T}$', color='blue')
        ax1.set_xlabel("$T$ $(^\circ C)$")
        ax1.set_ylim(0,0.09)
                
    elif var_mod=='light':
        ax11.set_ylabel('$f_{light}$', color='blue')
        ax1.set_xlabel("$PPFD$ $(\mu mol\,m^{-2}\,s^{-1})$")
        ax1.set_ylim(0,0.006)
        ax1.set_xlim(var_range[0], var_range[-1])
        
    elif var_mod=='vpd':
        ax11.set_ylabel('$f_{VPD}$', color='blue')
        ax1.set_xlabel("$VPD$ $(kPa)$")
        ax1.set_xlim(var_range[0], var_range[-1])
        ax1.set_ylim(0,3)
       
    
    ax11.grid(b=False)
    ax11.set_ylim(0,1.4)
    ax11.legend(fontsize='x-large')
    ax11.tick_params(axis='y', colors='blue')
    ax1.set_ylabel("PDF")
        
    ax1.legend(fontsize='x-large', loc='upper left')
