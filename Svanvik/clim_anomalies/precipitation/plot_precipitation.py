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
src_precip_svanvik = os.environ['DATA']+'/astra_data/observations/temp_prec_rad/svanvik_accu_precip_2009-2020.xlsx'
src_svanvik_prec_clim = os.environ['DATA']+'/astra_data/observations/temp_prec_rad/precip_svanvik_1992_2020.csv'
#src_precip_cru_clim =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/SCA_precip_flux_WFDE5_CRU_climatology.nc'
#src_precip_cru_clim_std =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/SCA_precip_flux_WFDE5_CRU_climatology_std.nc'
#src_precip_cru =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/2017-2018/SCA_precip_flux_WFDE5_CRU_*.nc'
# Load data
try:
    data_svanvik
except NameError:
    
    data_svanvik = pd.read_excel(src_precip_svanvik).dropna()
    data_svanvik[u"sum(precipitation_amount P1D)"] = [float(each.replace(',','.')) for each in data_svanvik[u"sum(precipitation_amount P1D)"].values]
    data_svanvik.index = [np.datetime64("%s-%s-%s" % (itime[-4:], itime[-7:-5], itime[:2])) for itime in data_svanvik['time'].dropna().values ]

    #data_cru_clim = xr.open_dataset(src_precip_cru_clim)
    #data_cru_clim_std = xr.open_dataset(src_precip_cru_clim_std)


    svanvik_precip = data_svanvik[u"sum(precipitation_amount P1D)"]

    data_svanvik_precip_clim = pd.read_csv(src_svanvik_prec_clim)
    data_svanvik_precip_clim.index = pd.date_range("%s" % data_svanvik_precip_clim['Time Measured'][0][:-3], "%s" % data_svanvik_precip_clim['Time Measured'].iloc[-1][:-3], freq='D')
    data_svanvik_precip_clim = data_svanvik_precip_clim.drop(columns=['Time Measured'])[:'2013']

    data_svanvik_precip_clim.loc[:,'day'] = data_svanvik_precip_clim.index.day.values
    data_svanvik_precip_clim.loc[:,'month'] = data_svanvik_precip_clim.index.month.values
    data_svanvik_precip_clim.loc[:,'year'] = data_svanvik_precip_clim.index.year.values

    data_svanvik.loc[:,'day'] = data_svanvik.index.day.values
    data_svanvik.loc[:,'month'] = data_svanvik.index.month.values
    data_svanvik.loc[:,'year'] = data_svanvik.index.year.values

# Save data to file

data_svanvik['2018'][u"sum(precipitation_amount P1D)"].to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_precipitation_2018.csv")
data_svanvik['2019'][u"sum(precipitation_amount P1D)"].to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_precipitation_2019.csv")
# Drop Feb 29 from climatology and reindex it
save_data = pd.DataFrame({'Precipitation (mm)':data_svanvik_precip_clim.groupby(['month','day']).mean().drop(data_svanvik_precip_clim.groupby(['month','day']).mean().loc[2,29,:].index)['RR'].values}, index=pd.date_range('2018-01-01', '2018-12-31 23:00', freq='1D'))
save_data.to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_precipitation-climatology.csv")

    #svanvik_precip_cru = data_cru_clim.sel(lat=station_location['Svanvik'].lat, lon=station_location['Svanvik'].lon, method='nearest')['Rainf']
    #svanvik_precip_cru_std = data_cru_clim_std.sel(lat=station_location['Svanvik'].lat, lon=station_location['Svanvik'].lon, method='nearest')['Rainf']

    #svanvik_ap_cru = (svanvik_precip_cru*60**2).resample(time='1D').sum()
    #svanvik_ap_cru_std = np.sqrt(((svanvik_precip_cru_std*60**2)**2).resample(time='1D').sum())

    #svanvik_ap_cru_test_list = []
    #svanvik_ap_cru_test_std_list = []
    #for iyear in range(2014,2020):
    #    if iyear != 2016:
    #        #print(iyear)
    #        svanvik_ap_cru_test = svanvik_ap_cru.drop(pd.date_range('2016-02-29',freq='D', periods=1), dim='time').copy()
    #        svanvik_ap_cru_test['time'] = pd.date_range('%d-01-01' % iyear, periods=365, freq='D')
    #        svanvik_ap_cru_test_std = svanvik_ap_cru_std.drop(pd.date_range('2016-02-29',freq='D', periods=1), dim='time').copy()
    #        svanvik_ap_cru_test_std['time'] = pd.date_range('%d-01-01' % iyear, periods=365, freq='D')
    #    else:
            #print(iyear, "ja")
    #        svanvik_ap_cru_test = svanvik_ap_cru.copy()
    #        svanvik_ap_cru_test['time'] = pd.date_range('%d-01-01' % iyear, periods=366, freq='D')
    #        svanvik_ap_cru_test_std = svanvik_ap_cru_std.copy()
    #        svanvik_ap_cru_test_std['time'] = pd.date_range('%d-01-01' % iyear, periods=366, freq='D')
        
    #    svanvik_ap_cru_test_list.append(svanvik_ap_cru_test)
    #    svanvik_ap_cru_test_std_list.append(svanvik_ap_cru_test_std)

    #svanvik_ap_cru_test_list = xr.concat(svanvik_ap_cru_test_list, dim='time')
    #svanvik_ap_cru_test_std_list = xr.concat(svanvik_ap_cru_test_std_list, dim='time')
    #svanvik_pull = ((svanvik_precip['2014':'2019']-svanvik_ap_cru_test_list.sel(time=svanvik_precip['2014':'2019'].index.values))/svanvik_ap_cru_test_std_list.sel(time=svanvik_precip['2014':'2019'].index.values))



#svanvik_precip_clim = data_svanvik_precip_clim.groupby(['year','month']).sum().groupby(['month']).mean()['RR']
#svanvik_precip_clim_std = data_svanvik_precip_clim.groupby(['year','month']).sum().groupby(['month']).std()['RR']

svanvik_precip_clim = data_svanvik_precip_clim.groupby(['month','day']).mean()['RR']
svanvik_precip_clim_std = data_svanvik_precip_clim.groupby(['month','day']).std()['RR']

# Plot it
#fig1 = plt.figure(1)
#fig1.canvas.set_window_title("plot_precipitation_anomalies_svanvik")
#ax11 = plt.subplot(211)
#ax12 = plt.subplot(212)

#svanvik_ap_cru_test_list.plot(ax=ax11, color='black', label="CRU reanalysis", ls='--')
#svanvik_precip['2014':].plot(ax=ax11, color='blueviolet', label='Svanvik', ls='-')
#(svanvik_ap_cru_test_list+svanvik_ap_cru_test_std_list).plot(ax=ax11, color='black', alpha=0.45)
#(svanvik_ap_cru_test_list-svanvik_ap_cru_test_std_list).plot(ax=ax11, color='black', alpha=0.45)

#ax11.legend(ncol=2, loc='upper left')
#ax11.set_xlabel("")
#ax11.set_ylabel("$Precip$ (mm day)")

#pull_plot = svanvik_pull.where((svanvik_pull.index.month>6)&(svanvik_pull.index.month<10)).dropna()
#pull_plot.plot(ax=ax12, ls='None', marker='x')

#ax12.axhline(2, color='red', ls=':')
#ax12.axhline(1, color='orange', ls=':')
#ax12.axhline(-2, color='red', ls=':')
#ax12.axhline(-1, color='orange', ls=':')

#for iyear in range(2014,2020):
#    print(iyear, pull_plot.where(pull_plot['%d' % iyear]<-2).dropna().count(), pull_plot.where(pull_plot['%d' % iyear]<-1).dropna().count())

#ax12.set_xlabel("Time (years)")
#ax12.set_ylabel("$\Delta_{clim} Precip / \sigma Precip_{clim}$")
#ax12.set_ylim(-2.5, 2.5)

#fig2 = plt.figure("precipitation_hist")
#ax21 = plt.subplot()

#svanvik_precip.hist(ax=ax21, bins=np.arange(0,50,0.1), density=True, histtype='step', color='blueviolet')
#ax21.hist((svanvik_precip_cru*60**2).resample(time='1D').sum(), density=True, bins=np.arange(0,50,0.1), histtype='step')

fig3 = plt.figure(3)
fig3.canvas.set_window_title("precipitation_signific")
ax31 = plt.subplot(211)
ax31.set_title("(a)")
ax32 = plt.subplot(212)
ax32.set_title("(b)")
for iyear, iax, icolor in zip((2018, 2019),(ax31, ax32), ('violet', 'purple')):
    tmp = data_svanvik['%d' % iyear].groupby(['month', 'day']).mean()
    test = ((tmp[u"sum(precipitation_amount P1D)"]-svanvik_precip_clim)/svanvik_precip_clim_std)
    (test.where(test>0.5).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100).plot.bar(ax=iax, color=icolor)
    (test.where(test<-0.5).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*-100).plot.bar(ax=iax, color=icolor)
    iax.set_xlabel("")
    iax.set_ylabel("")
    iax.set_ylim(-60, 60)
    iax.axhline(0, color='grey', ls=':')
    iax.axhline(30.9, color='grey', ls='--')
    iax.axhline(-30.9, color='grey', ls='--')
    
    iax.text(9.15,53, '$>+1/2\,\sigma$: %3.2f %s' % (test.where(test>0.5).dropna().size/float(test.size)*100, "%"), size='x-large')
    iax.text(9.25,-55, '$<-1/2\,\sigma$: %3.2f %s' % (test.where(test<-0.5).dropna().size/float(test.size)*100, "%"), size='x-large')
    iax.tick_params(labelrotation=0)
    iax.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)])
    
    print('%d %3.2f' % (iyear, test.where(test>0.5).dropna().size/float(test.size)*100))
    print('%d %3.2f' % (iyear, test.where(test<-0.5).dropna().size/float(test.size)*100))
   
    print('%d %s' % (iyear, test.where(test>0.5).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    print('%d %s' % (iyear, test.where(test<-0.5).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    
    
ax32.set_xlabel("Time (months)")
ax32.set_ylabel("#Days above $\pm 1/2\sigma_{clim}$ (%)", y=1)


# Show it
plt.show(block=False)
