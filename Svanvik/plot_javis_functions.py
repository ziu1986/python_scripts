import os, glob, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt   # Plotting data
import datetime as dt
import copy                       # for deep copy of object
from scipy.constants import *     # Get physics constants
from javis_model.javis_model import *
from mytools.plot_tools import *
from mytools.ozone_tools import VPD

# Data sources
src_svanvik = os.environ['DATA']+'/astra_data/observations/metdata_svanvik/Svanvik_temp_relhum_wind_*.csv'
src_svanvik_rad = os.environ['DATA']+'/astra_data/observations/metdata_svanvik/svanvik_glob_rad_*.csv'


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

def get_f_function(javis_model):
    f_temp = data_temp.iloc[:,0].apply(lambda x: javis_model.f_temp(x))
    f_vpd = vpd.apply(lambda x: javis_model.f_vpd(x))
    f_light = data_rad.iloc[:,0].apply(lambda x: javis_model.f_light(x))

    max_f = (np.maximum(javis_model.f_min, f_temp*f_vpd))
    prod_f = (f_light*np.maximum(javis_model.f_min, f_temp*f_vpd))
    
    return(f_temp, f_vpd, f_light, max_f, prod_f)

def get_variance(prod_f, **karg):
    var_type = karg.pop('type','noon')
    
    if var_type == 'noon':
        noon = gs_noon(prod_f)
        print('noon var', noon.var())
        return(noon.var())
    else:
        morning = gs_morning(prod_f)
        print('morning var', morning.var())
        return(morning.var())


def plot_f_functions(javis_model, fig_i, **karg):
    start_date = karg.pop('start', '2019-05')
    end_date = karg.pop('end', '2019-08')

    # Compute f_functions
    f_temp, f_vpd, f_light, max_f, prod_f = get_f_function(javis_model)
    
    # Plot it
    fig = plt.figure(fig_i, figsize=(10,12))
    fig.canvas.set_window_title("javis_funcs_%s" % javis_model.name)
    ax1 = plt.subplot(411)
    ax2 = plt.subplot(412)
    ax3 = plt.subplot(413)
    ax4 = plt.subplot(414)

    f_temp[start_date:end_date].plot(ax=ax1)
    f_vpd[start_date:end_date].plot(ax=ax2)
    f_light[start_date:end_date].plot(ax=ax3)

    max_f[start_date:end_date].plot(ax=ax3, color='red')
    (prod_f[start_date:end_date]*javis_model.gmax).plot(ax=ax4, color='blue')

    ax1.set_ylabel("f_temp")
    ax2.set_ylabel("f_vpd")
    ax3.set_ylabel("f_light")
    ax4.set_ylabel("$g_{sto}$ $(mmol m^{-2} s^{-1})$")
    ax4.set_xlabel("Time (months)")

    for ax in fig.axes[:-1]:
        ax.set_ylim(0,1)

def plot_optimal(results):
    fig = plt.figure(1, figsize=(10,12))
    fig.canvas.set_window_title("javis_func_opt")
    for i, icolor, iname in zip((0,1,2), ('darkgreen', 'lightgreen', 'green'), ('evergreen', 'birch', 'grassland')):
        ax = plt.subplot(3,1,i+1)
        ax.set_title(iname)
        ax.plot(results[::2][9*i:9*(i+1)], ls='None', marker='o', color=icolor, label='noon')
        ax.plot(results[1::2][9*i:9*(i+1)], ls='None', marker='s', color=icolor, label='morning')
    for ax in fig.axes:
        ax.set_ylim(0,0.11)
        ax.set_ylabel('Variance')
        ax.set_xticklabels((r'', r'mm', r'$\alpha_{+20}$', r'$\alpha_{-20}$',
                            r'boreal', r'$\alpha_{+20}$', r'$\alpha_{-20}$',
                            r'cold', r'$\alpha_{+20}$', r'$\alpha_{-20}$'))
        ax.legend()

    ax.set_xlabel('Categories')
    

def plot_temperature_histogram(temperature, **karg):
    fig = plt.figure()
    ax1 = plt.subplot()
    fig.canvas.set_window_title("javis_funcs_temp_hist")
    temperature.hist(ax=ax1, bins=np.arange(-40, 40), label='all years')
    temp_summer_all = temperature.where((temperature.index.month>=5)&(temperature.index.month<9))
    temp_summer_90 = temperature.where((temperature.index.month>=5)&(temperature.index.month<9))[:'2000']
    temp_summer_00 = temperature.where((temperature.index.month>=5)&(temperature.index.month<9))['2001':'2010']

    print("Mean all years summer: %2.3f\\n Mean 90s summer: %2.3f\\n Mean 00s summer: %2.3f" % (temp_summer_all.mean(), temp_summer_90.mean(), temp_summer_00.mean()))
    temp_summer_all.dropna().hist(ax=ax1, bins=np.arange(-40, 40), histtype='stepfilled', hatch='//', color='red', label='May-Sep')
    temp_summer_90.dropna().hist(ax=ax1, bins=np.arange(-40, 40), histtype='stepfilled', hatch='\\', color='blue', label='May-Sep (1992-2000)')
    temp_summer_00.dropna().hist(ax=ax1, bins=np.arange(-40, 40), histtype='step', hatch='\\', color='black', label='May-Sep (2001-2010)')

    ax1.set_xlabel("$T$ $(^\circ C)$")
    ax1.set_ylabel("Counts")
    ax1.legend(fontsize='large')

# main
alpha_var_decr_20 = {'evergreen':0.0075, 'birch':0.00525, 'grassland':0.01375}
alpha_var_incr_20 = {'evergreen':0.005, 'birch':0.0035, 'grassland':0.0092}

# Set up the different species
# Evergreen
evergreen = JavisModel('evergreen', Tmin=0, Tmax=200, Topt=20, fmin=0.1, Dmax=0.8, Dmin=2.8, alpha=0.006, gmax=125)
evergreen_boreal = JavisModel('evergreen_boreal', Tmin=0, Tmax=200, Topt=10, fmin=0.1, Dmax=0.8, Dmin=2.8, alpha=0.006, gmax=125)
evergreen_cold = JavisModel('evergreen_cold', Tmin=0, Tmax=200, Topt=15, fmin=0.1, Dmax=0.8, Dmin=2.8, alpha=0.006, gmax=125)
# Birch
birch = JavisModel('birch', Tmin=5, Tmax=200, Topt=20, fmin=0.1, Dmax=0.5, Dmin=2.7, alpha=0.0042, gmax=240)
birch_boreal = JavisModel('birch_boreal', Tmin=5, Tmax=200, Topt=12, fmin=0.1, Dmax=0.5, Dmin=2.7, alpha=0.0042, gmax=240)
birch_cold = JavisModel('birch_cold', Tmin=5, Tmax=200, Topt=15, fmin=0.1, Dmax=0.5, Dmin=2.7, alpha=0.0042, gmax=240)
# Grassland
grassland = JavisModel('grassland', Tmin=10, Tmax=36, Topt=24, fmin=0.1, Dmax=1.75, Dmin=4.5, alpha=0.011, gmax=190)
grassland_boreal = JavisModel('grassland_boreal', Tmin=0, Tmax=24, Topt=12, fmin=0.1, Dmax=1.75, Dmin=4.5, alpha=0.011, gmax=190)
grassland_cold = JavisModel('grassland_cold', Tmin=0, Tmax=36, Topt=15, fmin=0.1, Dmax=1.75, Dmin=4.5, alpha=0.011, gmax=190)


# Load data
data_temp = import_data(src_svanvik)
data_rad = import_data(src_svanvik_rad)
# Comput vpd
vpd = VPD(data_temp.iloc[:,1], data_temp.iloc[:,0])/kilo

# Clean up
plt.close('all')

# Loop through species and plot them
#for species, i in zip((evergreen, birch, grassland, evergreen_boreal, birch_boreal, grassland_boreal, evergreen_cold, birch_cold, grassland_cold), np.arange(1,10)):
#    plot_f_functions(species, i)

varaiance = []
# Loop through species and print variance
for species in (evergreen, evergreen_boreal, evergreen_cold, birch_cold, birch_boreal, birch, grassland, grassland_boreal, grassland_cold):
    for kind in ('evergreen', 'birch', 'grassland'):
        if kind in species.name:
            print(species.name)
            result = get_f_function(species)
            varaiance.append(get_variance(result[-1], type='noon'))
            varaiance.append(get_variance(result[-1], type='morning'))
            
            temp = copy.deepcopy(species)
            temp.name = temp.name + "+20"
            temp.alpha_light = alpha_var_incr_20[kind]
            result = get_f_function(temp)
            print(temp.name)
            varaiance.append(get_variance(result[-1], type='noon'))
            varaiance.append(get_variance(result[-1], type='morning'))
            
            temp.name = temp.name[:-3] + "-20"
            temp.alpha_light = alpha_var_decr_20[kind]
            result = get_f_function(temp)
            print(temp.name)
            varaiance.append(get_variance(result[-1], type='noon'))
            varaiance.append(get_variance(result[-1], type='morning'))
    
plot_optimal(varaiance)   

# Show it
plt.show(block=False)
