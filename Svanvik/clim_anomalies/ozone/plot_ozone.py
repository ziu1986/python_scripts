import os, sys, glob
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.ozone_tools import *
from mytools.station_info import *

# Clean up
plt.close('all')

# Source
src_svanvik = os.environ['DATA']+'/astra_data/observations/ozone/Svanvik/NO0047R.*.ozone.1y.1h.xls'
#src_svanvik_clim = os.environ['PY_SCRIPTS']+'/ozone_metrics/ozone_anomalies/obs_climatologies.pkl'
src_svanvik_corr = os.environ['DATA']+'/astra_data/input_data/DO3SE_input/svanvik_ozone_20*.csv'
src = os.environ['DATA']+'/astra_data/observations/ozone/Svanvik/*.nas'

# Read data
data_list = []
try:
    data_svanvik
except NameError:
    data_svanvik_obs = load_data(src)
    
    for file in sorted(glob.glob(src_svanvik)):
        data_list.append(pd.read_excel(file, index_col=0, header=0))
    data_svanvik = pd.concat(data_list)
    data_svanvik = data_svanvik.where(data_svanvik['O3_mugm-3']>=0)
    data_svanvik.loc[:,'day'] = data_svanvik.index.day.values
    data_svanvik.loc[:,'month'] = data_svanvik.index.month.values
    data_svanvik.loc[:,'hour'] = data_svanvik.index.hour.values

    data_svanvik_ozone = data_svanvik['O3_mugm-3'].dropna()/2.

    # Read gap filled ozone data and convert file data
    data_svanvik_corr = []
    for file in sorted(glob.glob(src_svanvik_corr)):
        tmp_data_svanvik = pd.read_csv(file)
        data_svanvik_corr.append(tmp_data_svanvik)
    # Concat data Svanvik data
    data_svanvik_corr = pd.concat(data_svanvik_corr)
    data_svanvik_corr = data_svanvik_corr.set_index('Year_Mnth_Date_Time(NMT)')
    data_svanvik_corr.index = pd.to_datetime(data_svanvik_corr.index)

    data_svanvik_corr.loc[:,'day'] = data_svanvik_corr.index.day.values
    data_svanvik_corr.loc[:,'month'] = data_svanvik_corr.index.month.values
    data_svanvik_corr.loc[:,'hour'] = data_svanvik_corr.index.hour.values

    # Load climatologies from observations
    #import pickle
    #with open(src_svanvik_clim, 'rb') as input:
    #    climatology_svanvik = pickle.load(input)
    #    yerr_mean_svanvik = pickle.load(input)

mapping = {2018:data_svanvik['2018'].dropna().groupby(data_svanvik['2018'].dropna().index.dayofyear).mean()['month'],
           2019:data_svanvik['2019'].dropna().groupby(data_svanvik['2019'].dropna().index.dayofyear).mean()['month'],
           20181:data_svanvik_corr['2018'].dropna().groupby(data_svanvik_corr['2018'].dropna().index.dayofyear).mean()['month']} 
# Get climatology
yozone_svanvik, yerr_svanvik, yerr_mean_svanvik = compute_climatology(data_svanvik_obs)
# Compute spline fits
from scipy.interpolate import UnivariateSpline
# Weights for smoothing
w_svanvik = 1/yerr_mean_svanvik
fitSpl_dmean_svanvik = UnivariateSpline(yozone_svanvik.index.values, yozone_svanvik.values, w=w_svanvik)
# Bias correction
bias_corr = 1.2

# Sampling from climatology
svanvik_ozone_clim =  {2018:fitSpl_dmean_svanvik(np.unique(data_svanvik_ozone['2018'].index.dayofyear.values)),
                       2019:fitSpl_dmean_svanvik(np.unique(data_svanvik_ozone['2019'].index.dayofyear.values)),
                       20181:fitSpl_dmean_svanvik(np.unique(data_svanvik_corr['2018'].index.dayofyear.values))}
svanvik_ozone_clim_std =  {2018:yerr_svanvik.iloc[np.unique(data_svanvik_ozone['2018'].index.dayofyear.values)-1],
                           2019:yerr_svanvik.iloc[np.unique(data_svanvik_ozone['2019'].index.dayofyear.values)-1],
                           20181:yerr_svanvik.iloc[np.unique(data_svanvik_corr['2018'].index.dayofyear.values)-1]}
# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_time")
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)
for iyear, iax, icolor in zip((2018, 2019),(ax11, ax12), ('violet', 'purple')):
    tmp = data_svanvik_ozone['%d' % iyear].groupby(data_svanvik_ozone['%d' % iyear].index.dayofyear).mean()
    tmp.plot(ax=iax, color=icolor, ls='None', marker='o')
    plot_error_bands(iax, np.arange(1,367), fitSpl_dmean_svanvik(np.arange(1,367))+bias_corr, error=yerr_mean_svanvik.values, color='grey', label="spline fit")
    plot_error_bands(iax, np.arange(1,367), yozone_svanvik.values+bias_corr, error=yerr_mean_svanvik.values, color='blueviolet', label='clim. mean')
    iax.set_title("%d" % iyear)
    iax.set_ylim(0,60)
    iax.set_xlim(0,366)
    iax.set_xlabel("")
    #iax.legend()

data_svanvik_corr['2018'].iloc[:,0].groupby(data_svanvik_corr['2018'].index.dayofyear).mean().plot(ax=ax11, marker='x')
ax11.set_xlabel("")
ax11.set_ylabel("$[O_3]$ (ppb)", y=-0.1)
ax12.set_xlabel("Time (days of year)")


probe_p = []
probe_n = []

text_p = []
text_n = []

fig3 = plt.figure(3)
fig3.canvas.set_window_title("ozone_signific")
ax31 = plt.subplot(211)
ax31.set_title("(a)")
ax32 = plt.subplot(212)
ax32.set_title("(b)")
for iyear, iax, icolor in zip((2018, 2019),(ax31, ax32, ax31), ('violet', 'purple', 'grey')):
    tmp = data_svanvik_ozone['%d' % iyear].groupby(data_svanvik_ozone['%d' % iyear].index.dayofyear).mean()
    
    test = pd.DataFrame({'O3': ((tmp-svanvik_ozone_clim[iyear]-bias_corr)/(svanvik_ozone_clim_std[iyear])).values, 'month':mapping[iyear].values}, index=mapping[iyear].index.values)
    
    probe_p.append(((test.where(test['O3']>1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size)))*100))
    probe_p[-1].reindex(np.arange(1,13)).plot.bar(ax=iax, color=icolor)
    probe_n.append(((test.where(test['O3']<-1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size)))*-100))
    probe_n[-1].reindex(np.arange(1,13)).plot.bar(ax=iax, color=icolor)
    
    iax.set_xlabel("")
    iax.set_ylabel("")
    iax.set_ylim(-60, 60)
    iax.axhline(0, color='grey', ls=':')
    iax.axhline(15.9, color='grey', ls='--')
    iax.axhline(-15.9, color='grey', ls='--')

    text_p.append(test.where(test['O3']>1).dropna().size/float(test.size)*100)
    text_n.append(test.where(test['O3']<-1).dropna().size/float(test.size)*100)
    
    iax.text(9.25,53, '$>+1\,\sigma$: %3.2f %s' % (text_p[-1], "%"), size='x-large')
    iax.text(9.5,-55, '$<-1\,\sigma$: %3.2f %s' % (text_n[-1], "%"), size='x-large')
    iax.tick_params(labelrotation=0)
    iax.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)])
    
    print('%d %3.2f' % (iyear, text_p[-1]))
    print('%d %3.2f' % (iyear, text_n[-1]))
   
    print('%d %s' % (iyear, test.where(test['O3']>1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    print('%d %s' % (iyear, test.where(test['O3']<-1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    
    
ax32.set_xlabel("Time (months)")
ax32.set_ylabel("#Days above $\pm 2\sigma_{clim}$ (%)", y=1)

# Get evaluation for corrected 2018
tmp = data_svanvik_corr.iloc[:,0]['2018'].groupby(data_svanvik_corr.iloc[:,0]['2018'].index.dayofyear).mean()
test = pd.DataFrame({'O3': ((tmp-svanvik_ozone_clim[20181]-bias_corr)/(svanvik_ozone_clim_std[20181])).values, 'month':mapping[20181].values}, index=mapping[20181].index.values)
    
probe_p.append(((test.where(test['O3']>1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size)))*100))
probe_n.append(((test.where(test['O3']<-1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size)))*-100))
    
fig4 = plt.figure(4, figsize=(12,8))
fig4.canvas.set_window_title("ozone_signific_1")
ax41 = plt.subplot()
ax41.axhspan(-15.9, 15.9, color='linen')

plot_data_p = pd.DataFrame({"2018":probe_p[0].reindex(np.arange(1,13)), "2019":probe_p[1].reindex(np.arange(1,13))})#, "2018_corr":probe_p[2].reindex(np.arange(1,13))})
plot_data_n = pd.DataFrame({"2018":probe_n[0].reindex(np.arange(1,13)), "2019":probe_n[1].reindex(np.arange(1,13))}) #, "2018_corr":probe_n[2].reindex(np.arange(1,13)),})

plot_data_p.plot.bar(ax=ax41, width=0.85, color=('violet', 'purple'))
plot_data_n.plot.bar(ax=ax41, width=0.85, color=('violet', 'purple'))


ax41.legend(["_","2018","2019","_","_"], loc='upper left')

ax41.plot(6-0.2, probe_p[2][7], marker='*', markersize=12)
#ax41.bar(6-0.2, probe_p[2][7], width=0.425, facecolor='None', edgecolor='violet', hatch='--')
ax41.set_ylim(-60, 60)

ax41.axhline(0, color='black', ls=':')
ax41.axhline(15.9, color='grey', ls='--')
ax41.axhline(-15.9, color='grey', ls='--')
ax41.axhspan(-15.9, 15.9, facecolor='None', edgecolor='black', hatch='//', alpha=0.5)

ax41.text(9.75, 55, '$>+1\,\sigma$', size='x-large', color='black')
ax41.text(9.25, 51, '2018: %2.2f %s' % (text_p[0], "%"), size='x-large', color='violet')
ax41.text(9.25, 47, '2019: %2.2f %s' % (text_p[1], "%"), size='x-large', color='purple')

ax41.text(9.75, -47, '$>-1\,\sigma$', size='x-large', color='black')
ax41.text(9.25, -51, '2018: %2.2f %s' % (text_n[0], "%"), size='x-large', color='violet')
ax41.text(9.25, -55, '2019: %2.2f %s' % (text_n[1], "%"), size='x-large', color='purple')

    
    
ax41.tick_params(labelrotation=0)
ax41.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)])
    
ax41.set_xlabel("Time (months)")
ax41.set_ylabel("#Days above $\pm 1\sigma_{clim}$ (%)")


# Show it
plt.show(block=False)
