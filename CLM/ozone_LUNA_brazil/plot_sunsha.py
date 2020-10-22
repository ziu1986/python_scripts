import os, sys, glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all, plot_error_bands, get_month_name, seconds_in_month
from mytools.ozone_tools import flunder
from mytools.clm_tools import *

# Clean up
plt.close("all")

# CLM simulations
# Reference - no ozone - ozone_luna_100
src = os.environ['CESM_RUN'] + '/archive/'
exp = ('brazil_2000_5.0.34',
       'brazil_2000_5.0.34_ozone_luna_40_wos',
       'brazil_2000_5.0.34_ozone_luna_20_wos')

data = {}

data.update({exp[0]:load_data(src + exp[0] + '/lnd/hist/*.nc', var=['PSNSHA', 'PSNSUN', 'GSSHA', 'GSSUN'])})

for iexp in exp[1:]:
    data.update({iexp:load_data(src + iexp + '/lnd/hist/*.nc', var=['PSNSHA', 'PSNSUN', 'GSSHA', 'GSSUN', 'O3UPTAKESHA', 'O3UPTAKESUN'])})

def plot(fig1, var, **karg):

    b_fit = karg.pop('fit', False)
    
    if b_fit:
        def funcP1d(x, slope, intercept):
            return(slope*x+1)

        from scipy.optimize import curve_fit
        par0 = (0,0)

        # Output
        f = open('fit_results_sunsha_%s.txt' % var, 'w')

    # Plot it
    
    fig1.canvas.set_window_title('sunsha_%s' % var)
    ax11 = plt.subplot(211)
    ax12 = plt.subplot(212)

    for i, marker, ls in zip((1,2), ('+', 'x'), ('--','-.')):
        for ax, mode in zip(fig1.axes, ('SHA', 'SUN')):
            ax.set_title(mode)
            test = (data[exp[i]]/data[exp[0]]).sel(time='2000')["%s%s" % (var, mode)].dropna(dim='time')
            test_time = test.time
            ax.plot(data[exp[i]].sel(time=test_time)['O3UPTAKE%s' % mode], test, ls='None', marker=marker, label="%s ppb" % exp[i][-6:-4])

            if b_fit:
                x = data[exp[i]].sel(time=test_time)['O3UPTAKE%s' % mode].values
                vals, covar = curve_fit(funcP1d, flunder(x), flunder(test.values), p0=par0)
                ax.plot(data[exp[i]].sel(time=test_time)['O3UPTAKE%s' % mode], funcP1d(flunder(x),*vals), ls=ls, color='black', label='fit')
                f.write("%s \n%s \n%s" % (exp[i], vals, covar))
        
    for ax in fig1.axes:
        ax.set_ylim(0.95,1.01)
        ax.set_xlim(0,135)
        ax.legend()
    if var == 'PSN':
        ax12.set_ylabel('$A_n^{O_3}/A_n^{ref}$', y=1)
    elif var == 'GS':
        ax12.set_ylabel('$g_{sto}^{O_3}/g_{sto}^{ref}$', y=1)
    
    ax12.set_xlabel('$\Sigma\Phi_{O_3}$ ($mmol\,m^{-2}$)')
    

    # Show it
    plt.show(block=False)
    
    
# Call
fig1 = plt.figure(1,figsize=(10,8))
plot(fig1, 'PSN', fit=True)
