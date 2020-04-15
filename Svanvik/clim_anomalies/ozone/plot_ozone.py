import os, sys, glob
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.station_info import *

# Clean up
plt.close('all')

# Source
src_svanvik = os.environ['DATA']+'/astra_data/observations/ozone/Svanvik/NO0047R.*.ozone.1y.1h.xls'
src_svanvik_clim = os.environ['PY_SCRIPTS']+'/ozone_metrics/ozone_anomalies/obs_climatologies.pkl'

# Read data
data_list = []
try:
    data_svanvik
except NameError:

    for file in sorted(glob.glob(src_svanvik)):
        data_list.append(pd.read_excel(file, index_col=0, header=0))
    data_svanvik = pd.concat(data_list)
    data_svanvik = data_svanvik.where(data_svanvik['O3_mugm-3']>=0)
    data_svanvik.loc[:,'day'] = data_svanvik.index.day.values
    data_svanvik.loc[:,'month'] = data_svanvik.index.month.values
    data_svanvik.loc[:,'hour'] = data_svanvik.index.hour.values

    data_svanvik_ozone = data_svanvik['O3_mugm-3'].dropna()/2.

    # Load climatologies from observations
    import pickle
    with open(src_svanvik_clim, 'rb') as input:
        climatology_svanvik = pickle.load(input)
        yerr_mean_svanvik = pickle.load(input)

# Sampling from climatology
svanvik_ozone_clim =  {2018:climatology_svanvik(np.unique(data_svanvik_ozone['2018'].index.dayofyear)),
                       2019:climatology_svanvik(np.unique(data_svanvik_ozone['2019'].index.dayofyear))}
svanvik_ozone_clim_std =  {2018:yerr_mean_svanvik(np.unique(data_svanvik_ozone['2018'].index.dayofyear)),
                           2019:yerr_mean_svanvik(np.unique(data_svanvik_ozone['2019'].index.dayofyear))}

mapping = {2018:data_svanvik['2018'].dropna().groupby(data_svanvik['2018'].dropna().index.dayofyear).mean()['month'],
           2019:data_svanvik['2019'].dropna().groupby(data_svanvik['2019'].dropna().index.dayofyear).mean()['month']}

fig3 = plt.figure(3)
fig3.canvas.set_window_title("ozone_signific")
ax31 = plt.subplot(211)
ax31.set_title("(a)")
ax32 = plt.subplot(212)
ax32.set_title("(b)")
for iyear, iax, icolor in zip((2018, 2019),(ax31, ax32), ('violet', 'purple')):
    tmp = data_svanvik_ozone['%d' % iyear].groupby(data_svanvik_ozone['%d' % iyear].index.dayofyear).mean()
    
    test = pd.DataFrame({'O3': ((tmp-svanvik_ozone_clim[iyear])/(svanvik_ozone_clim_std[iyear])).values, 'month':mapping[iyear].values}, index=mapping[iyear].index.values)
    
    (test.where(test['O3']>0.25).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100).plot.bar(ax=iax, color=icolor)
    (test.where(test['O3']<-0.25).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*-100).plot.bar(ax=iax, color=icolor)
    iax.set_xlabel("")
    iax.set_ylabel("")
    #iax.set_xlim(-0.5,12.5)
    #iax.set_xticklabels(np.arange(1,13))
    iax.set_ylim(-60, 60)
    iax.axhline(0, color='grey', ls=':')
    iax.axhline(15.9, color='grey', ls='--')
    iax.axhline(-15.9, color='grey', ls='--')
    
    iax.text(9.25,53, '$>+1\,\sigma$: %3.2f %s' % (test.where(test>1).dropna().size/float(test.size)*100, "%"), size='x-large')
    iax.text(9.5,-55, '$<-1\,\sigma$: %3.2f %s' % (test.where(test<-1).dropna().size/float(test.size)*100, "%"), size='x-large')
    iax.tick_params(labelrotation=0)
    #iax.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)])
    
    print('%d %3.2f' % (iyear, test.where(test>1).dropna().size/float(test.size)*100))
    print('%d %3.2f' % (iyear, test.where(test<-1).dropna().size/float(test.size)*100))
   
    print('%d %s' % (iyear, test.where(test>1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    print('%d %s' % (iyear, test.where(test<-1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    
    
ax32.set_xlabel("Time (months)")
ax32.set_ylabel("#Days above $\pm 1\sigma_{clim}$ (%)", y=1)

# Show it
plt.show(block=False)
