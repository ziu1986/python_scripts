import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mytools.plot_tools import print_all

# Clean up
plt.close('all')

# Source
case = "svanvik"
var = 'Lai' #'Psn'
src = os.environ['DATA'] + "/astra_data/observations/MODIS/" + case +"/statistics*" + var + "*landcover*.csv"

try:
    data_list
except NameError:
    # Read data
    data_list = []

    for ifile in sorted(glob.glob(src)):
        data = pd.read_csv(ifile)
        data.index = pd.to_datetime(data['date'])
        data_list.append(data)

# Plot it
#fig1 = plt.figure(1)
#ax11 = plt.subplot()

#for each, imarker, icolor in zip(data_list, ('^', 'v'), ('red', 'blue')):
#    each['mean'].plot(ax=ax11, ls='none', marker=imarker, color=icolor, markersize=10)

#ax11.set_xlabel("Time (date)")
#ax11.set_ylabel("$A_{net}$ $(kgC\,m^{-2})$")

fig2 = plt.figure(2)
fig2.canvas.set_window_title("modis_%s" % (var))
ax21 = plt.subplot(211)
ax22 = plt.subplot(212)

from scipy.optimize import curve_fit

def poly2(x, *m):
    return(m[0]*(x-m[1])**2+m[2])
    
for iyear, iax in zip((2018, 2019), (ax21, ax22)):
    print("Year: "); print(iyear)
    for each, imarker, icolor in zip(data_list, ('^', 'v'), ('red', 'blue')):
        tmp = each['mean']['%d' % iyear].groupby(each['mean']['%d' % iyear].index.dayofyear).mean()
        tmp.plot(ax=iax, ls='none', marker=imarker, color=icolor, markersize=10)
        p0 = [1,150,0.06]
        popt, pcov = curve_fit(poly2, tmp.dropna().index.values, tmp.dropna().values, p0)
        roots = np.roots((popt[0], 2*popt[0]*popt[1], popt[0]*popt[1]**2+popt[2]))
        print('Unweighted fit parameters:', popt)
        print('Covariance matrix:'); print(pcov)
        print('Roots: '); print(np.abs(roots))
        fit_range = np.arange(0,366)
        yfit = poly2(fit_range, *popt)
        iax.plot(fit_range, yfit)
    
for ax in fig2.axes:
    if var=='Psn':
        ax.set_ylim(0,0.08)
    elif var=='Lai':
        ax.set_ylim(0,5)
    ax.set_xlabel("")
ax22.set_xlabel("Time (doy)")
if var=="Psn":
    ax22.set_ylabel("$A_{net}$ ($kgC\,m^{-2}$)", y=1)
elif var=="Lai":
    ax22.set_ylabel("$LAI$ $(m^{2}\,m^{-2})$", y=1)
# Show it
plt.show(block=False)
