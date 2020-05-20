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

# Save data to file
#pd.DataFrame({'Q0 (Wm^2)':data_svanvik['2018'].groupby(['month','day','hour']).mean().values.flatten()}, index=data_svanvik['2018'].index.values).to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_global_irradiance_2018.csv")
#pd.DataFrame({'Q0 (Wm^2)':data_svanvik['2019'].groupby(['month','day','hour']).mean().values.flatten()}, index=data_svanvik['2019'].index.values).to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_global_irradiance_2019.csv")
# Drop Feb 29 from climatology and reindex it
#save_data = pd.DataFrame({'Q0 (Wm^2)':data_svanvik_clim.groupby(['month','day','hour']).mean().drop(data_svanvik_clim.groupby(['month','day','hour']).mean().loc[2,29,:].index).values.flatten()}, index=pd.date_range('2018-01-01', '2018-12-31 23:00', freq='1H'))
#save_data.to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_global_irradiance-climatology.csv")
    
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


probe_p = []
probe_n = []

text_p = []
text_n = []

fig3 = plt.figure(3)
fig3.canvas.set_window_title("global_rad_signific")
ax31 = plt.subplot(211)
ax31.set_title("(a)")
ax32 = plt.subplot(212)
ax32.set_title("(b)")
for iyear, iax, icolor in zip((2018, 2019),(ax31, ax32), ('violet', 'purple')):
    test = ((data_svanvik['%d' % iyear].groupby(['month','day','hour']).max()-data_svanvik_clim.groupby(['month','day','hour']).mean())/data_svanvik_clim.groupby(['month','day','hour']).std())

    probe_p.append((test.where(test>1).iloc[:,0].dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    probe_p[-1].plot.bar(ax=iax, color=icolor)
    probe_n.append((test.where(test<-1).iloc[:,0].dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*-100))
    probe_n[-1].plot.bar(ax=iax, color=icolor)
    iax.set_xlabel("")
    iax.set_ylabel("")
    iax.set_ylim(-60, 60)
    iax.axhline(0, color='grey', ls=':')
    iax.axhline(15.9, color='grey', ls='--')
    iax.axhline(-15.9, color='grey', ls='--')

    text_p.append(test.where(test>1).iloc[:,0].dropna().size/float(test.size)*100)
    text_n.append(test.where(test<-1).iloc[:,0].dropna().size/float(test.size)*100)
    
    iax.text(9.25,53, '$>+1\,\sigma$: %3.2f %s' % (text_p[-1], "%"), size='x-large')
    iax.text(9.5,-53, '$<-1\,\sigma$: %3.2f %s' % (text_n[-1], "%"), size='x-large')
    iax.tick_params(labelrotation=0)
    iax.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)])
    
    print('%d %3.2f' % (iyear, text_p[-1]))
    print('%d %3.2f' % (iyear, text_n[-1]))
   
    print('%d %s' % (iyear, test.where(test>1).iloc[:,0].dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    print('%d %s' % (iyear, test.where(test<-1).iloc[:,0].dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    
    
ax32.set_xlabel("Time (months)")
ax32.set_ylabel("#Days above $\pm 1\sigma_{clim}$ (%)", y=1)


fig4 = plt.figure(4, figsize=(12,8))
fig4.canvas.set_window_title("global_rad_signific_1")
ax41 = plt.subplot()
ax41.axhspan(-15.9, 15.9, color='linen')

plot_data_p = pd.DataFrame({"2018":probe_p[0], "2019":probe_p[1]})
plot_data_n = pd.DataFrame({"2018":probe_n[0], "2019":probe_n[1]})

plot_data_p.plot.bar(ax=ax41, width=0.85, color=('violet', 'purple'))
plot_data_n.plot.bar(ax=ax41, width=0.85, color=('violet', 'purple'))

ax41.legend(["_","2018","2019","_","_"], loc='upper left')

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

