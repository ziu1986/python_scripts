import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from scipy import fftpack
from read_sunspots import *

def test_smoothing_factor(x, y, **karg):
    weights = karg.pop('w', None)
    s_begin = karg.pop('s_begin', 1)
    s_end = karg.pop('s_end', len(weights))
    s_res = karg.pop('steps', 1)
    R2 = []
    spl_list = []
    for i in np.arange(s_begin, s_end, s_res):
        spl = UnivariateSpline(x,y, w=~np.isnan(y), s=i)
        spl_list.append(spl)
        R2.append(1-spl(x).var()/y.var())
    return(spl_list, np.arange(s_begin, s_end, s_res), R2)
    
    
# Clean up
plt.close('all')

# The data
nc_src = os.environ['DATA']+'/astra_data/observations/sunspots/SN_d_tot_V2.0.txt'
data = read_sunspots(nc_src)

# Spline fits
from scipy.interpolate import UnivariateSpline, LSQUnivariateSpline
x = np.arange(len(data.resample('1M').mean().index))
y = data['Ntot'].resample('1M').mean().copy()
y[np.isnan(y)] = 0.
stencils = (65,141,194,227,307,369,459,510,590,636,726,787,847,911,1005,1064,1136,1196,1263,1319,1385,1444,1514,1563,1633,1684,1756,1803,1908,1947,2025,2066,2144,2180,2289,2349)
stencils2 = (443,1011,1725)
spl = LSQUnivariateSpline(x,y, stencils, w=~np.isnan(y))
spl2 = LSQUnivariateSpline(x,y, stencils2, w=~np.isnan(y))

test = test_smoothing_factor(x,y, w=~np.isnan(y), s_begin=2.4e6, s_end=12.4e6, steps=0.5e6)

# Plot it
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("sunspot_freq_spectrum")

ax11 = plt.subplot(211)
data[['Ntot']].plot(ax=ax11, ls='None', marker='+', label='$N_{tot}^{daily}$')
#(data['Ntot'].resample('1M').mean().fillna(method='ffill')).plot(color='black', label='$N_{tot}^{monthly}$')
(data['Ntot'].fillna(method='ffill').resample('1M').mean()).interpolate(method='spline', order=3).plot(color='black')

ax11.set_xlabel("Time (year)")
ax11.set_ylabel("N$_{tot}^{ss}$")
ax11.legend()

ax12 = plt.subplot(212)
fft = fftpack.fft((data['Ntot'].resample('1M').mean().fillna(method='ffill')))
freqs = fftpack.fftfreq(len(fft))
ax12.stem(1/freqs/12, np.abs(fft))
    
ax12.set_ylabel("Amplitude")
ax12.set_xlabel("Frequency (years)")
ax12.set_xlim(0,205)

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("spline_manual_knots")
ax21 = plt.subplot(211)
ax21.plot(x,y, alpha=0.25, color='black', marker='x')
ax21.plot(x, spl(x), color='black', label='36k manual')
ax21.plot(x, test[0][0](x), color='blue', label='41k fit')
ax21.plot(x, spl2(x), color='red', label='3k manual')
ax21.set_xlabel("Time (months since begin of obs)")
ax21.set_ylabel("N$_{tot}^{ss}$")
ax21.legend()
ax22 = plt.subplot(212)
ax22.hist(np.sqrt((spl(x)-y)**2), histtype='step', bins=range(180), color='black')
ax22.hist(np.sqrt((test[0][0](x)-y)**2), histtype='step', bins=range(180), color='blue')
ax22.hist(np.sqrt((spl2(x)-y)**2), histtype='step', bins=range(180), color='red')
#ax22.scatter(x, np.sqrt((y-spl(x))**2), marker='x', label='residuals_manual')
#ax22.scatter(x, np.sqrt((y-test[0][0](x))**2), marker='+', label='residuals_49k')
#ax22.set_xlabel("Time (month since begin of obs)")
#ax22.set_ylabel("Residuals")
ax22.set_ylabel("Count")
ax22.set_xlabel("Residuals")

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("spline_smoothing_test")

ax31 = plt.subplot(211)
ax31.scatter(x,y, alpha=0.25, color='black', marker='x')
for each in test[0]:
    ax31.plot(x, each(x), label='N$_{knots}$=%d' % (len(each.get_knots())))

ax31.set_xlabel("Time (months since begin of obs)")
ax31.set_ylabel("N$_{tot}^{ss}$")
ax31.legend(ncol=5, fontsize='small')

#ax32 = plt.subplot(222)
#ax32.scatter(x,y, alpha=0.25, color='black', marker='x')
#ax32.plot(x,test[0][-4](x), color='red')

ax33 = plt.subplot(212)
ax33.plot(test[1], [len(each.get_knots()) for each in test[0]])
ax33.set_xlabel("s-factor")
ax33.set_ylabel("N$_{knots}$")

ax33t = ax33.twinx()
ax33t.set_ylabel("$R^2$", color='red')
ax33t.plot(test[1], test[2], color='red')
ax33t.tick_params(axis='y', labelcolor='red')


# Show it
plt.show(block=False)
