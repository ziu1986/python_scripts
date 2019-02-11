execfile('read_ebas.py')
def plot_month_span(ax):
    ax.axvspan(1,31, color='linen')
    ax.axvspan(31+28, 31+28+31, color='linen')
    ax.axvspan(31+28+31+30, 31+28+31+30+31, color='linen')
    ax.axvspan(31+28+31+30+31+30, 31+28+31+30+31+30+31, color='linen')
    ax.axvspan(31+28+31+30+31+30+31+31, 31+28+31+30+31+30+31+31+30, color='linen')
    ax.axvspan(31+28+31+30+31+30+31+31+30+31, 31+28+31+30+31+30+31+31+30+31+30, color='linen')
def plot_month_name(ax, ypos):
    xpos = (1, 31, 31+28, 31+28+31, 31+28+31+30, 31+28+31+30+31,
            31+28+31+30+31+30, 31+28+31+30+31+30+31, 31+28+31+30+31+30+31+31,
            31+28+31+30+31+30+31+31+30, 31+28+31+30+31+30+31+31+30+31, 31+28+31+30+31+30+31+31+30+31+30)
    for i in range(1,13):
        ax.text(xpos[i-1], ypos, get_month_name(i, length=3))
    
    
# Close the previous plots
plt.close('all')

# Directories of data
nc_src_o3 = os.environ['DATA']+'/astra_data/ctm_results/CTM3_oivind/osloctm_ozone*.nc'
nc_src_so2 = os.environ['DATA']+'/astra_data/ctm_results/CTM3_oivind/osloctm_so2*.nc'
nc_src_no2 = os.environ['DATA']+'/astra_data/ctm_results/CTM3_oivind/osloctm_no2*.nc'
src_karasjok_o3 = os.environ['DATA']+'/processed_data/observations/Karasjok/NO0055R.*ozone*.nas'
src_karasjok_so2 = os.environ['DATA']+'/processed_data/observations/Karasjok/NO0055R.*pack*.nas'
src_karasjok_no2 = os.environ['DATA']+'/processed_data/observations/Karasjok/NO0055R.*nitrogen_dioxide*.nas'


try:
    data_karasjok
except NameError:
    data_karasjok_o3 = []
    data_karasjok_so2 = []
    data_karasjok_no2 = []
    
    for file in sorted(glob.glob(src_karasjok_o3)):
        tmp = read_station_data(file)
        data_karasjok_o3.append(pd.Series(tmp['O3'],index=tmp['time']))
    data_karasjok_o3 = pd.concat(data_karasjok_o3)
    for file in sorted(glob.glob(src_karasjok_so2)):
        tmp = read_station_data(file, tracer=('SO2',))
        data_karasjok_so2.append(pd.Series(tmp['SO2'],index=tmp['time']))
    data_karasjok_so2 = pd.concat(data_karasjok_so2)
    for file in sorted(glob.glob(src_karasjok_no2)):
        tmp = read_station_data(file, tracer=('NO2',))
        data_karasjok_no2.append(pd.Series(tmp['NO2'],index=tmp['time']))
    data_karasjok_no2 = pd.concat(data_karasjok_no2)
try:
    data_o3
except NameError:
    data_o3 = read_data(nc_src_o3,var='O3',lev=(1,1))
    data_so2 = read_data(nc_src_so2,var='SO2',lev=(1,1))
    data_no2 = read_data(nc_src_no2,var='NO2',lev=(1,1))

# Data selection
sel_data_o3 = (data_o3.sel(lat=69,method='nearest').sel(lon=24,method='nearest')/1e-9)
sel_data_so2 =  (data_so2.sel(lat=69,method='nearest').sel(lon=24,method='nearest')/1e-9)
sel_data_no2 = (data_no2.sel(lat=69,method='nearest').sel(lon=24,method='nearest')/1e-9)

clim_o3 = sel_data_o3.groupby('time.dayofyear').mean()
clim_so2 = sel_data_so2.groupby('time.dayofyear').mean()
clim_no2 = sel_data_no2.groupby('time.dayofyear').mean()
climerr_o3 = sel_data_o3.groupby('time.dayofyear').std()/np.sqrt(sel_data_o3.groupby('time.dayofyear').sum()/clim_o3)
climerr_so2 = sel_data_so2.groupby('time.dayofyear').std()/np.sqrt(sel_data_so2.groupby('time.dayofyear').sum()/clim_so2)
climerr_no2 = sel_data_no2.groupby('time.dayofyear').std()/np.sqrt(sel_data_no2.groupby('time.dayofyear').sum()/clim_no2)
anom_o3 = sel_data_o3.groupby('time.dayofyear')-clim_o3
anom_so2 = sel_data_so2.groupby('time.dayofyear')-clim_so2
anom_no2 = sel_data_no2.groupby('time.dayofyear')-clim_no2

clim_karasjok_o3 = data_karasjok_o3.groupby(data_karasjok_o3.index.dayofyear).mean()
clim_karasjok_so2 = data_karasjok_so2.groupby(data_karasjok_so2.index.dayofyear).mean()
clim_karasjok_no2 = data_karasjok_no2.groupby(data_karasjok_no2.index.dayofyear).mean()
climerr_karasjok_o3 = data_karasjok_o3.groupby(data_karasjok_o3.index.dayofyear).std()/np.sqrt(data_karasjok_o3.groupby(data_karasjok_o3.index.dayofyear).sum()/clim_karasjok_o3)
climerr_karasjok_so2 = data_karasjok_so2.groupby(data_karasjok_so2.index.dayofyear).std()/np.sqrt(data_karasjok_so2.groupby(data_karasjok_so2.index.dayofyear).sum()/clim_karasjok_so2)
climerr_karasjok_no2 = data_karasjok_no2.groupby(data_karasjok_no2.index.dayofyear).std()/np.sqrt(data_karasjok_no2.groupby(data_karasjok_no2.index.dayofyear).sum()/clim_karasjok_no2)
#anom_karasjok_o3 = data_karasjok_o3.groupby(data_karasjok_o3.index.dayofyear).apply(lambda x: x - clim_karasjok_o3)
#anom_karasjok_so2 = data_karasjok_so2.apply(lambda x: x - clim_karasjok_so2)
#anom_karasjok_no2 = data_karasjok_no2.apply(lambda x: x - clim_karasjok_no2)

# Plotting
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("ebas_timeseries")
ax11 = plt.subplot(311)
ax12 = plt.subplot(312)
ax13 = plt.subplot(313)
sel_data_o3.plot(ax=ax11, alpha=0.15, color='blue', label='OsloCTM3')
sel_data_so2.plot(ax=ax12, alpha=0.15, color='blue', label='OsloCTM3')
sel_data_no2.plot(ax=ax13, alpha=0.15, color='blue', label='OsloCTM3')
data_karasjok_o3['1998-01-01':].plot(ax=ax11, marker='+', ls='none', color='grey', alpha=0.15, label='Karasjok')
data_karasjok_so2['1998-01-01':].plot(ax=ax12, marker='+', ls='none', color='grey', alpha=0.15, label='Karasjok')
data_karasjok_no2['1998-01-01':].plot(ax=ax13, marker='+', ls='none', color='grey', alpha=0.15, label='Karasjok')

for ax in fig1.axes[:-1]:
    ax.set_xlabel('')
for ax in fig1.axes:   
    ax.set_title('')
    ax.set_ylabel('%s (ppb)' % ax.get_ylabel())
    ax.legend()
    
fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("ebas_climatology_profile")
ax21 = plt.subplot(311)
ax22 = plt.subplot(312)
ax23 = plt.subplot(313)

clim_o3.plot(ax=ax21, color='blue', label='OsloCTM3')
clim_so2.plot(ax=ax22, color='blue', label='OsloCTM3')
clim_no2.plot(ax=ax23, color='blue', label='OsloCTM3')
plot_error_bands(ax21,
                 clim_o3.dayofyear,
                 clim_o3.data,
                 climerr_o3)
plot_error_bands(ax22,
                 clim_so2.dayofyear,
                 clim_so2.data,
                 climerr_so2)
plot_error_bands(ax23,
                 clim_no2.dayofyear,
                 clim_no2.data,
                 climerr_no2)

clim_karasjok_o3.plot(ax=ax21,
                      yerr=climerr_karasjok_o3,
                      marker='.', ls='none', color='grey', label='Karasjok')
clim_karasjok_so2.plot(ax=ax22,
                       yerr=climerr_karasjok_so2,
                       marker='.', ls='none', color='grey', label='Karasjok')
clim_karasjok_no2.plot(ax=ax23,
                       yerr=climerr_karasjok_no2,
                       marker='.', ls='none', color='grey', label='Karasjok')

for ax in fig2.axes[:-1]:
    ax.set_xlabel('')
for ax in fig2.axes:   
    ax.set_title('')
    ax.set_ylabel('%s (ppb)' % ax.get_ylabel())
    ax.set_xlim(0,367)
    plot_month_span(ax)
    ax.legend()
plot_month_name(ax21, 55)
ax21.set_ylim(0,60)
ax22.set_ylim(0,3)
ax23.set_ylim(0,3)

fig3 = plt.figure(3, figsize=(16,9))
fig3.canvas.set_window_title("ebas_anomalies")
ax31 = plt.subplot(311)
ax32 = plt.subplot(312)
ax33 = plt.subplot(313)
anom_o3.plot(ax=ax31, color='blue', label='OsloCTM3')
anom_so2.plot(ax=ax32, color='blue', label='OsloCTM3')
anom_no2.plot(ax=ax33, color='blue', label='OsloCTM3')
#anom_karasjok_o3.plot(ax=ax31, marker='.', ls='none', color='grey', label='Karasjok')
#anom_karasjok_so2.plot(ax=ax32, marker='.', ls='none', color='grey', label='Karasjok')
#anom_karasjok_no2.plot(ax=ax33, marker='.', ls='none', color='grey', label='Karasjok')

# Show it
plt.show(block=False)
