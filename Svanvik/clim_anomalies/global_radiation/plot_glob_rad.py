import os, sys
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mytools.plot_tools import *
from mytools.station_info import *


# Clean up
plt.close('all')

# Source
src_rad_svanvik_clim = os.environ['DATA']+'/astra_data/observations/temp_prec_rad/svanvik_glob_rad_2000_2017.csv'
src_rad_svanvik = os.environ['DATA']+'/astra_data/observations/temp_prec_rad/svanvik_glob_rad_2018_2019.csv'
# Load data
try:
    data_svanvik
except NameError:
    
    data_svanvik = pd.read_csv(src_rad_svanvik)
    data_svanvik_clim = pd.read_csv(src_rad_svanvik_clim)
    data_svanvik.index = pd.date_range("%s" % data_svanvik['Time measured'][0][:10], "%s" % data_svanvik['Time measured'].iloc[-1][:-3], freq='H')
    data_svanvik_clim.index = pd.date_range("%s" % data_svanvik_clim['Time measured'][0][:10], "%s" % data_svanvik_clim['Time measured'].iloc[-1][:-3], freq='H')
    data_svanvik = data_svanvik.drop(columns=['Time measured'])
    data_svanvik_clim = data_svanvik_clim.drop(columns=['Time measured'])

    data_svanvik_clim.loc[:,'hour'] = data_svanvik_clim.index.hour.values
    data_svanvik_clim.loc[:,'day'] = data_svanvik_clim.index.day.values
    data_svanvik_clim.loc[:,'month'] = data_svanvik_clim.index.month.values

    data_svanvik.loc[:,'hour'] = data_svanvik.index.hour.values
    data_svanvik.loc[:,'day'] = data_svanvik.index.day.values
    data_svanvik.loc[:,'month'] = data_svanvik.index.month.values
    
# Plot it
fig1 = plt.figure(1)
fig1.canvas.set_window_title("global_rad_clim")
ax11 = plt.subplot(211)
ax11.set_title("(a)")
xtime_hour = np.arange(1, data_svanvik_clim.groupby(['month','day','hour']).max().size+1)
ax11.plot(xtime_hour, data_svanvik_clim.groupby(['month','day','hour']).max(), color='red', label='max', alpha=0.5)
ax11.plot(xtime_hour, data_svanvik_clim.groupby(['month','day','hour']).mean(), color='black', label='mean', ls=':', alpha=0.75)
ax11.plot(xtime_hour, data_svanvik_clim.groupby(['month','day','hour']).min(), color='blue', label='min', ls='-', alpha=1)

ax12 = plt.subplot(223, sharex=ax11)
ax12.set_title("(b)")

ax12.plot(xtime_hour, (data_svanvik['2018'].groupby(['month','day','hour']).max()-data_svanvik_clim.groupby(['month','day','hour']).max()), color='red', ls='-', alpha=0.5)
ax12.plot(xtime_hour, (data_svanvik['2018'].groupby(['month','day','hour']).min()-data_svanvik_clim.groupby(['month','day','hour']).min()), color='blue', ls='-')
ax12.plot(xtime_hour, (data_svanvik['2018'].groupby(['month','day','hour']).mean()-data_svanvik_clim.groupby(['month','day','hour']).mean()), color='black', ls=':', alpha=0.75)

ax12.plot(xtime_hour[12::24*10], data_svanvik_clim.groupby(['month','day','hour']).std().groupby(['month','day']).max()[::10], color='cyan', ls='-', linewidth=3, label='$1\sigma_{clim}$', alpha=0.9)
#ax12.plot(xtime_hour[12::24*10], data_svanvik_clim.groupby(['month','day','hour']).std().groupby(['month','day']).min()[::10], color='grey', ls='-', linewidth=3)

ax12.plot(xtime_hour[12::24*10], -(data_svanvik_clim.groupby(['month','day','hour']).std().groupby(['month','day']).max())[::10], color='cyan', ls='-', linewidth=3, alpha=0.9)


ax13 = plt.subplot(224, sharex=ax11)
ax13.set_title("(c)")

ax13.plot(xtime_hour, (data_svanvik['2019'].groupby(['month','day','hour']).max()-data_svanvik_clim.groupby(['month','day','hour']).max()), color='red', ls='-', alpha=0.5)
ax13.plot(xtime_hour, (data_svanvik['2019'].groupby(['month','day','hour']).min()-data_svanvik_clim.groupby(['month','day','hour']).min()), color='blue', ls='-')
ax13.plot(xtime_hour, (data_svanvik['2019'].groupby(['month','day','hour']).mean()-data_svanvik_clim.groupby(['month','day','hour']).mean()), color='black', ls=':', alpha=0.75)

ax13.plot(xtime_hour[12::24*10], data_svanvik_clim.groupby(['month','day','hour']).std().groupby(['month','day']).max()[::10], color='cyan', ls='-', linewidth=3, label='$1\sigma_{clim}$', alpha=0.9)
ax13.plot(xtime_hour[12::24*10], -(data_svanvik_clim.groupby(['month','day','hour']).std().groupby(['month','day']).max())[::10], color='cyan', ls='-', linewidth=3, alpha=0.9)


for ax in fig1.axes:
    ax.legend()
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("")
    ybound = ax.get_ybound()[1]-50
    ax.text(1*24, ybound, "Jan")
    ax.text((31+1)*24, ybound,"Feb")
    ax.text((31+29+1)*24, ybound,"Mar")
    ax.text((31+29+31+1)*24, ybound,"Apr")
    ax.text((31+29+31+30+1)*24, ybound,"May")
    ax.text((31+29+31+30+31+1)*24, ybound, "Jun")
    ax.text((31+29+31+30+31+30+1)*24, ybound, "Jul")
    ax.text((31+29+31+30+31+30+31+1)*24, ybound, "Aug")
    ax.text((31+29+31+30+31+30+31+31+1)*24, ybound, "Sep")
    ax.text((31+29+31+30+31+30+31+31+30+1)*24, ybound, "Oct")
    ax.text((31+29+31+30+31+30+31+31+30+31+1)*24, ybound, "Nov")
    ax.text((31+29+31+30+31+30+31+31+30+31+30+1)*24, ybound, "Dec")
ax11.set_ylabel("Global radiation ($W\,m^{-2}$)")
#ax11.set_xlabel("")
ax12.set_ylabel("$\Delta$Global radiation ($W\,m^{-2}$)")


fig2 = plt.figure(2)
fig2.canvas.set_window_title("global_rad_clim_hist")
ax21 = plt.subplot(211)

((data_svanvik['2018'].groupby(['month','day','hour']).mean()-data_svanvik_clim.groupby(['month','day','hour']).mean())).iloc[:,0].hist(ax=ax21, bins=400, color='violet', label='2018', density=True)
((data_svanvik['2019'].groupby(['month','day','hour']).mean()-data_svanvik_clim.groupby(['month','day','hour']).mean())).iloc[:,0].hist(ax=ax21, bins=400, color='purple', histtype='step', label='2019', density=True)
#(data_svanvik_clim.groupby(['month','day']).max()/data_svanvik_clim.groupby(['month','day']).std()).iloc[:,0].plot(ax=ax21)
#(data_svanvik['2018'].groupby(['month','day']).max()/data_svanvik_clim.groupby(['month','day']).std()).iloc[:,0].plot(ax=ax21)
#(data_svanvik['2019'].groupby(['month','day']).max()/data_svanvik_clim.groupby(['month','day']).std()).iloc[:,0].plot(ax=ax21)

ax22 = plt.subplot(212)
((data_svanvik['2018'].groupby(['month','day','hour']).mean()-data_svanvik_clim.groupby(['month','day','hour']).mean())/data_svanvik_clim.groupby(['month','day','hour']).std()).iloc[:,0].hist(ax=ax22, bins=np.arange(-4.05,4.1,0.1), color='violet', label='2018', density=True)
((data_svanvik['2019'].groupby(['month','day','hour']).mean()-data_svanvik_clim.groupby(['month','day','hour']).mean())/data_svanvik_clim.groupby(['month','day','hour']).std()).iloc[:,0].hist(ax=ax22, bins=np.arange(-4.05,4.1,0.1), color='purple', histtype='step', label='2019', density=True)

for ax in fig2.axes:
    ax.legend()
ax21.set_xlim(-100, 100)
ax21.set_ylim(0,0.01)
ax21.set_ylabel("Probability density", y=0)
ax21.set_xlabel("$\Delta_{iyear-clim} Q_0$ ($W\,m^{-2}$)")
ax22.set_xlabel("$\Delta_{iyear-clim} Q_0/\sigma_{clim}$")

fig3 = plt.figure(3)
fig3.canvas.set_window_title("global_rad_signific")
ax31 = plt.subplot(211)
ax31.set_title("(a)")
ax32 = plt.subplot(212)
ax32.set_title("(b)")
for iyear, iax, icolor in zip((2018, 2019),(ax31, ax32), ('violet', 'purple')):
    test = ((data_svanvik['%d' % iyear].groupby(['month','day','hour']).max()-data_svanvik_clim.groupby(['month','day','hour']).mean())/data_svanvik_clim.groupby(['month','day','hour']).std())
    (test.where(test>1).iloc[:,0].dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100).plot.bar(ax=iax, color=icolor)
    (test.where(test<-1).iloc[:,0].dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*-100).plot.bar(ax=iax, color=icolor)
    iax.set_xlabel("")
    iax.set_ylabel("")
    iax.set_ylim(-60, 60)
    iax.axhline(0, color='grey', ls=':')
    iax.axhline(15.9, color='grey', ls='--')
    iax.axhline(-15.9, color='grey', ls='--')

    iax.text(9.25,53, '$>+1\,\sigma$: %3.2f %s' % (test.where(test>1).iloc[:,0].dropna().size/float(test.size)*100, "%"), size='x-large')
    iax.text(9.5,-53, '$<-1\,\sigma$: %3.2f %s' % (test.where(test<-1).iloc[:,0].dropna().size/float(test.size)*100, "%"), size='x-large')
    iax.tick_params(labelrotation=0)
    iax.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)])
    
    print('%d %3.2f' % (iyear, test.where(test>1).iloc[:,0].dropna().size/float(test.size)*100))
    print('%d %3.2f' % (iyear, test.where(test<-1).iloc[:,0].dropna().size/float(test.size)*100))
   
    print('%d %s' % (iyear, test.where(test>1).iloc[:,0].dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    print('%d %s' % (iyear, test.where(test<-1).iloc[:,0].dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    
    
ax32.set_xlabel("Time (months)")
ax32.set_ylabel("#Days above $\pm 1\sigma_{clim}$ (%)", y=1)


   

# Show it
plt.show(block=False)

