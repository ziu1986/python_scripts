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
src_temp_svanvik = os.environ['DATA']+'/astra_data/observations/temp_prec_rad/svanvik_temp_deg_2013-2019.xls'
src_temp_svanvik_clim = os.environ['DATA']+'/astra_data/observations/temp_prec_rad/temp_svanvik*'
#src_temp_cru_clim =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/SCA_Tair_WFDE5_CRU_climatology.nc'
#src_temp_cru_clim_std =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/SCA_Tair_WFDE5_CRU_climatology_std.nc'
#src_temp_cru =  os.environ['DATA']+'/nird_data/reanalysis/Copernicus/CRU/2017-2018/SCA_Tair_WFDE5_CRU_*.nc'
# Load data
data_list = []
try:
    data_svanvik
except NameError:
    for file in sorted(glob.glob(src_temp_svanvik_clim)):
        print(file)
        data_temp_svanvik_clim = pd.read_csv(file)
        data_temp_svanvik_clim.index = pd.date_range("%s" % data_temp_svanvik_clim['Time measured'][0][:-3], "%s" % data_temp_svanvik_clim['Time measured'].iloc[-1][:-3], freq='H')
        data_temp_svanvik_clim = data_temp_svanvik_clim.drop(columns=['Time measured'])
        data_list.append(data_temp_svanvik_clim)
        
    data_temp_svanvik_clim = pd.concat(data_list)
    
    data_svanvik = pd.read_excel(src_temp_svanvik)
    data_svanvik.index = data_svanvik[u'Fra-tid']
    data_svanvik.loc[:,'day'] = data_svanvik.index.day.values
    data_svanvik.loc[:,'month'] = data_svanvik.index.month.values
    data_svanvik.loc[:,'hour'] = data_svanvik.index.hour.values

    #data_cru_clim = xr.open_dataset(src_temp_cru_clim)
    #data_cru_clim_std = xr.open_dataset(src_temp_cru_clim_std)

    data_temp_svanvik_clim.loc[:,'day'] = data_temp_svanvik_clim.index.day.values
    data_temp_svanvik_clim.loc[:,'month'] = data_temp_svanvik_clim.index.month.values
    data_temp_svanvik_clim.loc[:,'hour'] = data_temp_svanvik_clim.index.hour.values

    svanvik_temp = (273.15+data_svanvik[u"Svanvik | Ambient Temperature | degC"].where((data_svanvik[u"Svanvik | Ambient Temperature | degC"]>-50) & (data_svanvik[u"Svanvik | Ambient Temperature | degC"]<50)).dropna())
    #svanvik_temp_cru = data_cru_clim.sel(lat=station_location['Svanvik'].lat, lon=station_location['Svanvik'].lon, method='nearest')['Tair']
    #svanvik_temp_cru_std = data_cru_clim_std.sel(lat=station_location['Svanvik'].lat, lon=station_location['Svanvik'].lon, method='nearest')['Tair']

    #svanvik_temp_cru_test_list = []
    #svanvik_temp_cru_test_std_list = []
    #for iyear in range(2014,2020):
    #    if iyear != 2016:
    #        #print(iyear)
    #        svanvik_temp_cru_test = svanvik_temp_cru.drop(pd.date_range('2016-02-29',freq='H', periods=24), dim='time').copy()
    #        svanvik_temp_cru_test['time'] = pd.date_range('%d-01-01' % iyear, periods=365*24, freq='H')
    #        svanvik_temp_cru_test_std = svanvik_temp_cru_std.drop(pd.date_range('2016-02-29',freq='H', periods=24), dim='time').copy()
    #        svanvik_temp_cru_test_std['time'] = pd.date_range('%d-01-01' % iyear, periods=365*24, freq='H')
    #    else:
    #        #print(iyear, "ja")
    #        svanvik_temp_cru_test = svanvik_temp_cru.copy()
    #        svanvik_temp_cru_test['time'] = pd.date_range('%d-01-01' % iyear, periods=366*24, freq='H')
    #        svanvik_temp_cru_test_std = svanvik_temp_cru_std.copy()
    #        svanvik_temp_cru_test_std['time'] = pd.date_range('%d-01-01' % iyear, periods=366*24, freq='H')
        
    #    svanvik_temp_cru_test_list.append(svanvik_temp_cru_test)
    #    svanvik_temp_cru_test_std_list.append(svanvik_temp_cru_test_std)

    #    svanvik_temp_cru_test_list = xr.concat(svanvik_temp_cru_test_list, dim='time')
    #    svanvik_temp_cru_test_std_list = xr.concat(svanvik_temp_cru_test_std_list, dim='time')
    #    svanvik_pull = ((svanvik_temp['2014':'2019']-svanvik_temp_cru_test_list.sel(time=svanvik_temp['2014':'2019'].index.values))/svanvik_temp_cru_test_std_list.sel(time=svanvik_temp['2014':'2019'].index.values))

    svanvik_temp_clim = data_temp_svanvik_clim.groupby(['month','day','hour']).mean().iloc[:,0]+273.15
    svanvik_temp_clim_std = data_temp_svanvik_clim.groupby(['month','day','hour']).std().iloc[:,0]

svanvik_temp_clim_deg = data_temp_svanvik_clim.groupby(['month','day','hour']).mean().iloc[:,0]

# Save data to file
#data_svanvik['2018']['Svanvik | Ambient Temperature | degC'].to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_temperature_2018.csv")
#data_svanvik['2019']['Svanvik | Ambient Temperature | degC'].to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_temperature_2019.csv")
# Drop Feb 29 from climatology and reindex it
#save_data = pd.DataFrame({'Temperature (degC)':svanvik_temp_clim_deg.drop(svanvik_temp_clim_deg.loc[2,29,:].index).values}, index=pd.date_range('2018-01-01', '2018-12-31 23:00', freq='1H'))
#save_data.to_csv(os.environ['DATA']+'/DO3SE_input/'+"svanvik_temperature-climatology.csv")


# Plot it
#fig1 = plt.figure(1)
#fig1.canvas.set_window_title("plot_temperature_anomalies_svanvik")
#ax11 = plt.subplot(211)
#ax12 = plt.subplot(212)

#svanvik_temp.plot(ax=ax11, color='blueviolet', label='Svanvik')
#svanvik_temp_cru_test_list.plot(ax=ax11, color='black', label="CRU reanalysis")
#(svanvik_temp_cru_test_list+svanvik_temp_cru_test_std_list).plot(ax=ax11, color='black', alpha=0.45)
#(svanvik_temp_cru_test_list-svanvik_temp_cru_test_std_list).plot(ax=ax11, color='black', alpha=0.45)

#ax11.legend(ncol=2, loc='upper left')

#svanvik_pull.plot(ax=ax12)

#ax12.axhline(4, color='red', ls=':')
#ax12.axhline(2, color='orange', ls=':')
#ax12.axhline(-4, color='red', ls=':')
#ax12.axhline(-2, color='orange', ls=':')

#for iyear in range(2014,2020):
#    print(iyear, svanvik_pull.where(svanvik_pull['%d' % iyear]>4).dropna().count(), svanvik_pull.where(svanvik_pull['%d' % iyear]>2).dropna().count())

#ax11.set_xlabel("Time (years)")
#ax11.set_ylabel("$T_{2m}$ (K)")
#ax12.set_xlabel("Time (months)")
#ax12.set_ylabel("$\Delta_{clim} T_{2m} / \sigma T_{2m, clim}$")

probe_p = []
probe_n = []

text_p = []
text_n = []

fig3 = plt.figure(3)
fig3.canvas.set_window_title("temperature_signific")
ax31 = plt.subplot(211)
ax31.set_title("(a)")
ax32 = plt.subplot(212)
ax32.set_title("(b)")
for iyear, iax, icolor in zip((2018, 2019),(ax31, ax32), ('violet', 'purple')):
    tmp = (data_svanvik['%d' % iyear].where((data_svanvik['%d' % iyear][u"Svanvik | Ambient Temperature | degC"]>-50) & (data_svanvik['%d' % iyear][u"Svanvik | Ambient Temperature | degC"]<50)).dropna()).groupby(['month','day','hour']).max()[u"Svanvik | Ambient Temperature | degC"]+273.15
    test = ((tmp-svanvik_temp_clim)/svanvik_temp_clim_std)
    
    probe_p.append((test.where(test>1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    probe_n.append((test.where(test<-1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*-100))
    probe_p[-1].plot.bar(ax=iax, color=icolor)
    probe_n[-1].plot.bar(ax=iax, color=icolor)
    iax.set_xlabel("")
    iax.set_ylabel("")
    iax.set_ylim(-60, 60)
    iax.axhline(0, color='grey', ls=':')
    iax.axhline(15.9, color='grey', ls='--')
    iax.axhline(-15.9, color='grey', ls='--')

    text_p.append(test.where(test>1).dropna().size/float(test.size)*100)
    text_n.append(test.where(test<-1).dropna().size/float(test.size)*100)
              
    iax.text(9.25,53, '$>+1\,\sigma$: %3.2f %s' % (text_p[-1], "%"), size='x-large')
    iax.text(9.5,-55, '$<-1\,\sigma$: %3.2f %s' % (text_n[-1], "%"), size='x-large')
    iax.tick_params(labelrotation=0)
    iax.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)])
    
    print('%d %3.2f' % (iyear, text_p[-1]))
    print('%d %3.2f' % (iyear, text_n[-1]))
   
    print('%d %s' % (iyear, test.where(test>1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    print('%d %s' % (iyear, test.where(test<-1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    
    
ax32.set_xlabel("Time (months)")
ax32.set_ylabel("#Days above $\pm 1\sigma_{clim}$ (%)", y=1)

fig4 = plt.figure(4, figsize=(12,8))
fig4.canvas.set_window_title("temperature_signific_1")
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
