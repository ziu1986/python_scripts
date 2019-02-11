# Load modules and EBAS reading routine
execfile('read_ebas')
# Close the previous plots
plt.close('all')

# Directories of data
src_jergul = os.environ['DATA']+'/astra_data/ctm_results/observations/Jergul/NO0030R.*ozone*.nas'
src_karasjok = os.environ['DATA']+'/astra_data/ctm_results/observations/Karasjok/NO0055R.*ozone*.nas'
src_svanvik = os.environ['DATA']+'/astra_data/ctm_results/observations/Svanvik/NO0047R.*ozone*.nas'
nc_src = os.environ['DATA']+'/astra_data/ctm_results/CTM3_oivind/osloctm_ozone*.nc'
# Loop through EBAS data and transform them to pandas timeseries
data_jergul = []
data_karasjok = []
data_svanvik = []
for file in sorted(glob.glob(src_jergul)):
    tmp = read_station_data(file)
    data_jergul.append(pd.Series(tmp['O3'],index=tmp['time']))
for file in sorted(glob.glob(src_karasjok)):
    tmp = read_station_data(file)
    data_karasjok.append(pd.Series(tmp['O3'],index=tmp['time']))
for file in sorted(glob.glob(src_svanvik)):
    tmp = read_station_data(file)
    data_svanvik.append(pd.Series(tmp['O3'],index=tmp['time']))
try:
    data
except NameError:
    data = read_data(nc_src,var='O3',lev=(1,1))
# Concatenate the lists
data_jergul = pd.concat(data_jergul)
data_karasjok = pd.concat(data_karasjok)
data_svanvik = pd.concat(data_svanvik)
# Select the data

# Plotting
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("ozone_ebas_1997")
ax11 = plt.subplot()
#ax12 = plt.subplot(212)
(data.sel(lat=69,method='nearest').sel(lon=24,method='nearest')/1e-9).sel(time=slice('1997-01-01','1998-01-01')).plot(label='OsloCTM3')
data_jergul['1997-01-01':'1998-01-01'].plot(marker='x', ls='none', color='red', alpha=0.75, label='Jergul')
data_karasjok['1997-01-01':'1998-01-01'].plot(marker='+', ls='none', color='orange', alpha=0.75, label='Karasjok')
#data_svanvik['1997-01-01':'1998-01-01'].plot(marker='v', ls='none', color='blueviolet', alpha=0.75, label='Svanvik')
ax11.set_xlabel('Time')
ax11.set_ylabel('$O_3$ (ppb)')
ax11.legend()

fig2 = plt.figure(2,figsize=(16,9))
b_sc = True
b_mm = True
b_sel = True
ncols = 4
fig2.canvas.set_window_title("ozone_ebas_timeseries_sampled")
ax21 = plt.subplot()
selection_jergul = data_jergul.where(
    (data_jergul.index.month>=5)&
    (data_jergul.index.month<8)&
    (data_jergul.data>=40))
selection_karasjok = data_karasjok.where(
    (data_karasjok.index.month>=5)&
    (data_karasjok.index.month<8)&
    (data_karasjok.data>=40))
selection_svanvik = data_svanvik.where(
    (data_svanvik.index.month>=5)&
    (data_svanvik.index.month<8)&
    (data_svanvik.data>=40))
selection_ctm = (data/1e-9).where(
    (data.time.dt.month>=5)&
    (data.time.dt.month<8)&
    ((data/1e-9)>=40))

selection2_jergul = data_jergul.where(
    (data_jergul.index.month>=5)&
    (data_jergul.index.month<8)&
    (data_jergul.index.hour>=8)&
    (data_jergul.index.hour<=20)&
    (data_jergul.data>=40))
selection2_karasjok = data_karasjok.where(
    (data_karasjok.index.month>=5)&
    (data_karasjok.index.month<8)&
    (data_karasjok.index.hour>=8)&
    (data_karasjok.index.hour<=20)&
    (data_karasjok.data>=40))
selection2_svanvik = data_svanvik.where(
    (data_svanvik.index.month>=5)&
    (data_svanvik.index.month<8)&
    (data_svanvik.index.hour>=8)&
    (data_svanvik.index.hour<=20)&
    (data_svanvik.data>=40))
selection2_ctm = (data/1e-9).where(
    (data.time.dt.month>=5)&
    (data.time.dt.month<8)&
    (data.time.dt.hour>=8)&
    (data.time.dt.hour<=20)&
    ((data/1e-9)>=40))

data_jergul.plot(marker='x', ls='none', color='grey', alpha=0.15, label='Jergul')
data_karasjok.plot(marker='+', ls='none', color='grey', alpha=0.15, label='Karasjok')
data_svanvik.plot(marker='', ls='-', color='blueviolet', alpha=0.15, label='Svanvik')
(data.sel(lat=69,method='nearest').sel(lon=24,method='nearest')/1e-9).plot(color='blue', alpha=0.15,label='OsloCTM3')

if b_mm:
    data_jergul[:'1996-12-31'].resample('1M', convention='start').mean().plot(marker='', ls='-', color='red', label='Jergul mm')
    data_karasjok.resample('1M').mean().plot(marker='', ls='-', color='orange', label='Karasjok mm')
    data_svanvik.resample('1M').mean().plot(marker='', ls='-', color='blueviolet', label='Svanvik mm')
    (data.sel(lat=69,method='nearest').sel(lon=24,method='nearest')/1e-9).resample(time='1M').mean().plot(color='blue',label='OsloCTM3 mm')
    ncols = 5
if b_sel:
    selection_jergul.plot(marker='.', ls='none', color='red', alpha=0.75, label='Jergul $\mathrm{(\geq 40\,ppb)}$')
    selection_karasjok.plot(marker='.', ls='none', color='orange', alpha=0.75, label='Karasjok $\mathrm{(\geq 40\,ppb)}$')
    selection_svanvik.plot(marker='.', ls='none', color='blueviolet', alpha=0.15, label='Svanvik $\mathrm{(\geq 40\,ppb)}$')
    selection_ctm.sel(lat=69,method='nearest').sel(lon=24,method='nearest').plot(marker='.', ls='none', color='blue', alpha=0.75, label='OsloCTM3 $\mathrm{(\geq 40\,ppb)}$')
if b_sc:
    ax21.axvspan(date2num(dt.datetime.strptime('1986-09','%Y-%m')),date2num(dt.datetime.strptime('1996-08','%Y-%m')),color='linen')
    ax21.axvspan(date2num(dt.datetime.strptime('2008-12','%Y-%m')),date2num(dt.datetime.strptime('2011-01','%Y-%m')),color='linen')
    ax21.plot(date2num(dt.datetime.strptime('1989-11','%Y-%m')),(-5),ls='',marker='s')
    ax21.text(date2num(dt.datetime.strptime('1989-11','%Y-%m')),(-3),"maximum solar cycle 22")
    ax21.plot(date2num(dt.datetime.strptime('2001-11','%Y-%m')),(-5),ls='',marker='s')
    ax21.text(date2num(dt.datetime.strptime('2001-11','%Y-%m')),(-3),"maximum solar cycle 23")

ax21.set_xlabel('Time')
ax21.set_ylabel('$O_3$ (ppb)')
ax21.legend(ncol=ncols)

fig3 = plt.figure(3,figsize=(16,9))
fig3.canvas.set_window_title("ozone_ebas_aot40vstime")
ax31 = plt.subplot()
cumsum_ctm = (selection_ctm.sel(lat=69,method='nearest').sel(lon=24,method='nearest').resample(time='1H').interpolate()-40).resample(time='D').reduce(np.nansum)
cumsum_jergul = (selection_jergul-40).resample('D').sum()
cumsum_karasjok = (selection_karasjok-40).resample('D').sum()
cumsum2_ctm = (selection2_ctm.sel(lat=69,method='nearest').sel(lon=24,method='nearest').resample(time='1H').interpolate()-40).resample(time='D').reduce(np.nansum)
cumsum2_jergul = (selection2_jergul-40).resample('D').sum()
cumsum2_karasjok = (selection2_karasjok-40).resample('D').sum()
cumsum_ctm.plot(zorder=3, color='blue',label='OsloCTM3')
cumsum_jergul.plot(zorder=2, marker='x', ls='none', color='red',label='Jergul')
cumsum_karasjok.plot(zorder=2, marker='+', ls='none', color='orange',label='Karasjok')

ax31.set_xlabel("Time (year)")
ax31.set_ylabel("AOT40 (ppb)")

ax31.legend()

fig4 = plt.figure(4, figsize=(16,9))
fig4.canvas.set_window_title("ozone_ebas_occurance_hist")
ax41 = plt.subplot()

filtered_data_jergul = []
filtered_data_karasjok = []
filtered_data_sim = []
filtered_data2_jergul = []
filtered_data2_karasjok = []
filtered_data2_sim = []

from pandas.tseries.offsets import *
for aot, date in zip(cumsum_jergul, cumsum_jergul.index):
    if (date.month > 4) & (date.month < 8):
        if pd.Timestamp('%s'%(date.year)).is_leap_year:
            filtered_data_jergul.append((aot, date.dayofyear-121))
        else:
            filtered_data_jergul.append((aot, date.dayofyear-120))
        
for aot, date in zip(cumsum_karasjok, cumsum_karasjok.index):
    if (date.month > 4) & (date.month < 8):
        if pd.Timestamp('%s'%(date.year)).is_leap_year:
            filtered_data_karasjok.append((aot, date.dayofyear-121))
        else:
            filtered_data_karasjok.append((aot, date.dayofyear-120))
            
for aot, date in zip(cumsum_ctm.data, cumsum_ctm.time):
    if (date.dt.month > 4) & (date.dt.month < 8):
        if pd.Timestamp('%s'%(date.dt.year.data)).is_leap_year:
            filtered_data_sim.append((aot,date.dt.dayofyear.data-121))
        else:
            filtered_data_sim.append((aot,date.dt.dayofyear.data-120))

for aot, date in zip(cumsum2_jergul, cumsum2_jergul.index):
    if (date.month > 4) & (date.month < 8):
        if pd.Timestamp('%s'%(date.year)).is_leap_year:
            filtered_data2_jergul.append((aot, date.dayofyear-121))
        else:
            filtered_data2_jergul.append((aot, date.dayofyear-120))
        
for aot, date in zip(cumsum2_karasjok, cumsum2_karasjok.index):
    if (date.month > 4) & (date.month < 8):
        if pd.Timestamp('%s'%(date.year)).is_leap_year:
            filtered_data2_karasjok.append((aot, date.dayofyear-121))
        else:
            filtered_data2_karasjok.append((aot, date.dayofyear-120))
            
for aot, date in zip(cumsum2_ctm.data, cumsum2_ctm.time):
    if (date.dt.month > 4) & (date.dt.month < 8):
        if pd.Timestamp('%s'%(date.dt.year.data)).is_leap_year:
            filtered_data2_sim.append((aot,date.dt.dayofyear.data-121))
        else:
            filtered_data2_sim.append((aot,date.dt.dayofyear.data-120))
            
hist_data_jergul = np.array(zip(*filtered_data_jergul)[1])[np.where(np.array(zip(*filtered_data_jergul)[0])>0)]
violin_data_jergul = np.array(zip(*filtered_data_jergul)[0]).reshape(len(filtered_data_jergul)/92,92)
hist_data_karasjok = np.array(zip(*filtered_data_karasjok)[1])[np.where(np.array(zip(*filtered_data_karasjok)[0])>0)]
violin_data_karasjok = np.array(zip(*filtered_data_karasjok)[0]).reshape(len(filtered_data_karasjok)/92,92)
hist_data_sim = np.array(zip(*filtered_data_sim)[1])[np.where(np.array(zip(*filtered_data_sim)[0])>0)]
violin_data_sim = np.array(zip(*filtered_data_sim)[0]).reshape(len(filtered_data_sim)/92,92)

hist_data2_jergul = np.array(zip(*filtered_data2_jergul)[1])[np.where(np.array(zip(*filtered_data2_jergul)[0])>0)]
violin_data2_jergul = np.array(zip(*filtered_data2_jergul)[0]).reshape(len(filtered_data2_jergul)/92,92)
hist_data2_karasjok = np.array(zip(*filtered_data2_karasjok)[1])[np.where(np.array(zip(*filtered_data2_karasjok)[0])>0)]
violin_data2_karasjok = np.array(zip(*filtered_data2_karasjok)[0]).reshape(len(filtered_data2_karasjok)/92,92)
hist_data2_sim = np.array(zip(*filtered_data2_sim)[1])[np.where(np.array(zip(*filtered_data2_sim)[0])>0)]
violin_data2_sim = np.array(zip(*filtered_data2_sim)[0]).reshape(len(filtered_data2_sim)/92,92)
ax41.axvspan(1,31,alpha=0.25,color='grey')
ax41.axvspan(62,92,alpha=0.25,color='grey')

hist_sim = ax41.hist(hist_data_sim,
                     bins=np.arange(1,93),
                     weights=np.repeat(1/13.,hist_data_sim.size),
                     label='OsloCTM3 (1997-01-01 $\mathrm{-}$ 2010-12-31)')

hist_jergul = ax41.hist(hist_data_jergul,
                        bins=np.arange(1,93),
                        histtype='step',
                        #color='red',
                        linewidth='2',
                        edgecolor='red', hatch="/",
                        weights=np.repeat(1/8.,hist_data_jergul.size),
                        label="Jergul (*1988-04-01$\mathrm{ - \dagger}$1997-01-07)")
hist_karasjok = ax41.hist(hist_data_karasjok,
                          bins=np.arange(1,93),
                          histtype='step',
                          #color='orange',
                          linewidth='1.75',
                          edgecolor='orange', hatch="\\",
                          weights=np.repeat(1/12.,hist_data_karasjok.size),
                          label="Karasjok (*1997-02-01$\mathrm{ - \dagger}$2010-02-28)")

ax41.set_ylabel("Count/N$_{years}$")
ax41.set_xlabel("Time (Days starting from May 1st)")
ax41.set_xlim(0,94)
ax41.set_ylim(0,1)
ax41.text(2,0.95, get_month_name(5), size='large')
ax41.text(32,0.95, get_month_name(6), size='large')
ax41.text(63,0.95, get_month_name(7), size='large')
ax41.legend(loc='upper right')

fig5 = plt.figure(5, figsize=(16,9))
fig5.canvas.set_window_title("ozone_ebas_violins")
ax51 = plt.subplot()
ax51.axvspan(1,31,alpha=0.25,color='grey')
ax51.axvspan(62,92,alpha=0.25,color='grey')
viol_jergul = ax51.violinplot(violin_data_jergul,
                positions=np.arange(0.75,92),
                points=92,
                showmeans=False,
                showmedians=False)
for pc in viol_jergul['bodies']:
    pc.set_facecolor('red')
    pc.set_edgecolor('red')
viol_jergul['cbars'].set_color('red')

viol_karasjok = ax51.violinplot(violin_data_karasjok,
                positions=np.arange(1,93),
                points=92,
                showmeans=False,
                showmedians=False)
for pc in viol_karasjok['bodies']:
    pc.set_facecolor('orange')
    pc.set_edgecolor('orange')
viol_karasjok['cbars'].set_color('orange')
viol_sim = ax51.violinplot(violin_data_sim,
                positions=np.arange(1.25,93),
                points=92,
                showmeans=False,
                showmedians=False)
for pc in viol_sim['bodies']:
    pc.set_facecolor('blue')
    pc.set_edgecolor('blue')
viol_sim['cbars'].set_color('blue')
ax51.plot(np.arange(0.75,92),violin_data_jergul.mean(axis=0),color='red',marker='x', ls='', label='Jergul mean')
ax51.plot(np.arange(1.,93),violin_data_karasjok.mean(axis=0),color='orange',marker='x', ls='', label='Karasjok mean')
ax51.plot(np.arange(1.25,93),violin_data_sim.mean(axis=0),color='blue',marker='x', ls='', label='OsloCTM3 mean')
ax51.plot(np.arange(0.75,92),np.median(violin_data_jergul,axis=0),color='mistyrose',marker='.', ls='', label='Jergul median')
ax51.plot(np.arange(1.,93),np.median(violin_data_karasjok,axis=0),color='peachpuff',marker='.', ls='', label='Karasjok median')
ax51.plot(np.arange(1.25,93),np.median(violin_data_sim,axis=0),color='lightblue',marker='.', ls='', label='OsloCTM3 median')
ax51.set_xlabel("Time (Days starting from May 1st)")
ax51.set_ylabel("$AOT_{40}$ (ppb)")
ax51.text(2,715, get_month_name(5), size='large')
ax51.text(32,715, get_month_name(6), size='large')
ax51.text(63,715, get_month_name(7), size='large')
ax51.set_xlim(0,93)
ax51.legend(loc='upper right')

fig6 = plt.figure(6, figsize=(16,9))
fig6.canvas.set_window_title("ozone_ebas_histAOT40")
ax61 = plt.subplot()
bins = [0,]
for i in np.arange(1,650,5):
    bins.append(i)
hist2_sim = ax61.hist(violin_data_sim.flat,
                      bins=bins,
                      #normed=True,
                      #weights=np.repeat(1/13.,hist_data_sim.size),
                      label='OsloCTM3 (1997-01-01 $\mathrm{-}$ 2010-12-31)')
hist2_jergul = ax61.hist(violin_data_jergul.flat,
                         bins=bins,
                         #normed=True,
                         histtype='step',
                         linewidth='2',
                         edgecolor='red', hatch="/",
                         #weights=np.repeat(1/8.,hist_data_jergul.size),
                         label="Jergul (*1988-04-01$\mathrm{ - \dagger}$1997-01-07)")
hist2_karasjok = ax61.hist(violin_data_karasjok.flat,
                           bins=bins,
                           #normed=True,
                           histtype='step',
                           linewidth='1.75',
                           edgecolor='orange', hatch="\\",
                           #weights=np.repeat(1/12.,hist_data_karasjok.size),
                           label="Karasjok (*1997-02-01$\mathrm{ - \dagger}$2010-02-28)")
ax61.set_ylabel("Counts")
ax61.set_xlabel("$AOT_{40}$ (ppb)")
ax61.set_ylim(0,90)
ax61.legend(loc='upper right')

fig7 = plt.figure(7,figsize=(16,9))
fig7.canvas.set_window_title("ozone_ebas_hist_time")
ax71 = plt.subplot(321)
ax71.set_title("OsloCTM3",x=0.9,y=0.85)
ax72 = plt.subplot(323)
ax72.set_title("Jergul",x=0.9,y=0.85)
ax73 = plt.subplot(325)
ax73.set_title("Karasjok",x=0.9,y=0.85)
ax74 = plt.subplot(322)
ax75 = plt.subplot(324)
ax76 = plt.subplot(326)
for ax in fig7.axes[:3]:
    ax.axvspan(1,31,alpha=0.25,color='grey')
    ax.axvspan(62,92,alpha=0.25,color='grey')
    ax.set_ylim(0,1)
    ax.set_xlim(1,93)
for ax in fig7.axes[3:]:
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    ax.plot((0,1),(0,1),ls='--',color='black')
ax71.hist(hist_data_sim,
          bins=np.arange(1,93),
          color='grey',
          weights=np.repeat(1/13.,hist_data_sim.size),
          label='[0-24]')
ax72.hist(hist_data_jergul,
          bins=np.arange(1,93),
          color='grey',
          weights=np.repeat(1/8.,hist_data_jergul.size),
          label="[0-24]")
ax73.hist(hist_data_karasjok,
          bins=np.arange(1,93),
          color='grey',
          weights=np.repeat(1/12.,hist_data_karasjok.size),
          label="[0-24]")

hist12_sim = ax71.hist(hist_data2_sim,
                       bins=np.arange(1,93),
                       histtype='step',
                       linewidth='2',
                       edgecolor='blue', hatch="/",
                       weights=np.repeat(1/13.,hist_data2_sim.size),
                       label='[8-20]')
hist12_jergul = ax72.hist(hist_data2_jergul,
                          bins=np.arange(1,93),
                          histtype='step',
                          linewidth='2',
                          edgecolor='red', hatch="/",
                          weights=np.repeat(1/8.,hist_data2_jergul.size),
                          label="[8-20]")
hist12_karasjok = ax73.hist(hist_data2_karasjok,
                            bins=np.arange(1,93),
                            histtype='step',
                            linewidth='1.75',
                            edgecolor='orange', hatch="/",
                            weights=np.repeat(1/12.,hist_data2_karasjok.size),
                            label="[8-20]")
ax74.plot(hist_sim[0], hist12_sim[0],ls='', marker='x')
ax75.plot(hist_jergul[0], hist12_jergul[0],ls='', marker='x',color='red')
ax76.plot(hist_karasjok[0], hist12_karasjok[0],ls='', marker='x',color='orange')
ax71.text(2,0.9, get_month_name(5), size='large')
ax71.text(32,0.9, get_month_name(6), size='large')
ax71.text(63,0.9, get_month_name(7), size='large')
ax72.set_ylabel("Count/N$_{years}$")
ax73.set_xlabel("Time (Days starting from May 1st)")
ax75.set_ylabel("(Count/N$_{years}$)$_{[0-24]}$")
ax76.set_xlabel("(Count/N$_{years}$)$_{[8-20]}$")
for ax in fig7.axes[:3]:
    ax.legend(loc='center right')

fig8 = plt.figure(8,figsize=(16,9))
fig8.canvas.set_window_title("ozone_ebas_hist_timeintegral_corr")
ax81 = plt.subplot(321)
ax81.set_title("OsloCTM3",x=0.9,y=0.85)
ax82 = plt.subplot(323)
ax82.set_title("Jergul",x=0.9,y=0.85)
ax83 = plt.subplot(325)
ax83.set_title("Karasjok",x=0.9,y=0.85)
ax84 = plt.subplot(322)
ax85 = plt.subplot(324)
ax86 = plt.subplot(326)
for ax in fig8.axes[:3]:
    ax.set_ylim(0,90)
    ax.set_xlim(0,650)
for ax in fig8.axes[3:]:
    ax.set_ylim(0,650)
    ax.set_xlim(0,650)
    ax.plot((0,650),(0,600),ls=':',color='black')
ax81.hist(violin_data_sim.flat,
          bins=bins,
          color='grey',
          #weights=np.repeat(1/13.,hist_data_sim.size),
          label='[0-24]')
ax82.hist(violin_data_jergul.flat,
          bins=bins,
          color='grey',
          #weights=np.repeat(1/8.,hist_data_jergul.size),
          label="[0-24]")
ax83.hist(violin_data_karasjok.flat,
          bins=bins,
          color='grey',
          #weights=np.repeat(1/12.,hist_data_karasjok.size),
          label="[0-24]")

hist22_sim = ax81.hist(violin_data2_sim.flat,
                       bins=bins,
                       histtype='step',
                       linewidth='2',
                       edgecolor='blue', hatch="/",
                       #weights=np.repeat(1/13.,hist_data2_sim.size),
                       label='[8-20]')
hist22_jergul = ax82.hist(violin_data2_jergul.flat,
                          bins=bins,
                          histtype='step',
                          linewidth='2',
                          edgecolor='red', hatch="/",
                          #weights=np.repeat(1/8.,hist_data2_jergul.size),
                          label="[8-20]")
hist22_karasjok = ax83.hist(violin_data2_karasjok.flat,
                            bins=bins,
                            histtype='step',
                            linewidth='1.75',
                            edgecolor='orange', hatch="/",
                            #weights=np.repeat(1/12.,hist_data2_karasjok.size),
                            label="[8-20]")
ax84.plot(np.array(zip(*filtered_data2_sim)[0]), np.array(zip(*filtered_data_sim)[0]), ls='', marker='x')
ax85.plot(np.array(zip(*filtered_data2_jergul)[0]), np.array(zip(*filtered_data_jergul)[0]), ls='', marker='x',color='red')
ax86.plot(np.array(zip(*filtered_data2_karasjok)[0]), np.array(zip(*filtered_data_karasjok)[0]), ls='', marker='x',color='orange')

fit2_sim = linregress(np.array(zip(*filtered_data2_sim)[0]), np.array(zip(*filtered_data_sim)[0]))
ffit2_sim = np.poly1d((fit2_sim[0],fit2_sim[1]))
ax84.plot((0,650),ffit2_sim((0,650)),color='blue',ls='--')
ax84.text(500,600,"$R^2 = %1.2f$" % (fit2_sim[2]**2))
fit2_jergul= linregress(np.array(zip(*filtered_data2_jergul)[0]), np.array(zip(*filtered_data_jergul)[0]))
ffit2_jergul = np.poly1d((fit2_jergul[0],fit2_jergul[1]))
ax85.plot((0,650),ffit2_jergul((0,650)),color='red',ls='--')
ax85.text(500,600,"$R^2 = %1.2f$" % (fit2_jergul[2]**2))
fit2_karasjok = linregress(np.array(zip(*filtered_data2_karasjok)[0]), np.array(zip(*filtered_data_karasjok)[0]))
ffit2_karasjok = np.poly1d((fit2_karasjok[0],fit2_karasjok[1]))
ax86.plot((0,650),ffit2_karasjok((0,650)),color='orange',ls='--')
ax86.text(500,600,"$R^2 = %1.2f$" % (fit2_karasjok[2]**2))

ax82.set_ylabel("Count")
ax83.set_xlabel("$AOT_{40}$ (ppb)")
ax85.set_ylabel("$AOT_{40}^{[0-24]}$ (ppb)")
ax86.set_xlabel("$AOT_{40}^{[8-20]}$ (ppb)")

for ax in fig8.axes[:3]:
    ax.legend(loc='center right')

#from scipy import fftpack
fig9 = plt.figure(9,figsize=(16,9))
fig9.canvas.set_window_title("ozone_ebas_resampling_test")
ax91 = plt.subplot()
data_karasjok['2006-06-13':'2006-06-14'].plot(color='orange',label='Karasjok')
data_karasjok.resample('1H').interpolate()[5::6]['2006-06-13':'2006-06-14'].plot(ls='',marker='s',label='Karasjok 6h sample')
data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate()['2006-06-13':'2006-06-14'].plot(ls=':',label='Karasjok resampled')
ax91.fill_between(data_karasjok['2006-06-13':'2006-06-14'].index,
                  data_karasjok['2006-06-13':'2006-06-14'].data,
                  data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate()['2006-06-13':'2006-06-14'],
                  where=data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate()['2006-06-13':'2006-06-14']>=data_karasjok['2006-06-13':'2006-06-14'].data,
                  facecolor='green', interpolate=True)
ax91.fill_between(data_karasjok['2006-06-13':'2006-06-14'].index,
                  data_karasjok['2006-06-13':'2006-06-14'].data,
                  data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate()['2006-06-13':'2006-06-14'],
                  where=data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate()['2006-06-13':'2006-06-14']<=data_karasjok['2006-06-13':'2006-06-14'].data,
                  facecolor='red', interpolate=True)
(data.sel(lat=69,method='nearest').sel(lon=24,method='nearest').sel(time=slice('2006-06-13','2006-06-14'))/1e-9).plot(color='blue',label='OsloCTM3')
ax91.set_ylim()
ax91.set_ylabel("$O_3$ (ppb)")
ax91.set_xlabel("Time")
ax91.legend()

fig10 = plt.figure(10,figsize=(16,9))
fig10.canvas.set_window_title("ozone_ebas_resampling_test_AOT40diff")
ax101 = plt.subplot(221)
ax102 = plt.subplot(222)
ax103 = plt.subplot(223)
ax104 = plt.subplot(224)
selection_test = data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate().where(
    (data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate().index.month>=5)&
    (data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate().index.month<8)&
    (data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate()>=40))
cumsum_test = (selection_test-40).resample('D').sum()
selection2_test = data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate().where(
    (data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate().index.month>=5)&
    (data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate().index.month<8)&
    (data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate().index.hour>=8)&
    (data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate().index.hour<=20)&
    (data_karasjok.resample('1H').interpolate()[5::6].resample('1H').interpolate()>=40))
cumsum2_test = (selection2_test-40).resample('D').sum()

ax101.plot(cumsum_karasjok, cumsum_test, ls='', marker='x', color='orange', label='[0-24]')
ax103.plot(cumsum2_karasjok, cumsum2_test, ls='', marker='x', color='orange', label='[0-8]')
hist102 = ax102.hist((cumsum_test-cumsum_karasjok.where((cumsum_test.index.month >=5) & (cumsum_test.index.month<8))).dropna(), bins=np.arange(-40,41,1))
hist104 = ax104.hist((cumsum2_test-cumsum2_karasjok.where((cumsum_test.index.month >=5) & (cumsum_test.index.month<8))).dropna(), bins=np.arange(-40,41,1))

for ax in fig10.axes[::2]:
    ax.set_xlim(-5,705)
    ax.set_ylim(-5,705)
    ax.plot((0,700),(0,700), ls='--', color='black')
    ax.legend()
for ax in fig10.axes[1::2]:
    ax.set_ylim(0,100)
    ax.set_xlim(-40,40)
ax101.set_xlabel("")
ax103.set_xlabel("AOT$_{40}^{origres}$ (ppb)")
ax101.set_ylabel("")
ax103.set_ylabel("AOT$_{40}^{resampled}$ (ppb)", y=1)
ax102.set_xlabel("")
ax104.set_xlabel("$\Delta$AOT$_{40}$ (ppb)")
ax102.set_ylabel("")
ax104.set_ylabel("count", y=1)
 
#ax101.set_ylabel("$\Delta AOT_{40}$ (ppb)")
#ax101.set_xlabel("Time (years)")
#ax101.set_ylim(-40,40)
#ax101.legend()

from scipy import fftpack
fig11 = plt.figure(11,figsize=(16,9))
fig11.canvas.set_window_title("ozone_ebas_freq_spectrum")
ax111 = plt.subplot(211)
ax112 = plt.subplot(212)
fft_sim = fftpack.fft((data.sel(lat=69,method='nearest').sel(lon=24,method='nearest')/1e-9).resample(time='1M').mean())
freqs_sim = fftpack.fftfreq(len(fft_sim))
fft_karasjok = fftpack.fft(pd.concat((data_jergul,data_karasjok)).resample('1M').mean().fillna(method='ffill'))
freqs_karasjok = fftpack.fftfreq(len(fft_karasjok))
ax111.stem(1/freqs_sim/12, np.abs(fft_sim))
ax111.set_title("OsloCTM3",x=0.9,y=0.85)
ax112.stem(1/freqs_karasjok/12, np.abs(fft_karasjok))
ax112.set_title("Jergul $\mathrm{\cap}$ Karasjok",x=0.9,y=0.85)
for ax in fig11.axes:
    ax.set_xlim(0,25)
    
ax111.set_ylabel("")
ax112.set_ylabel("Amplitude",y=1)
ax111.set_xlabel("")
ax112.set_xlabel("Time (years)")
    
# Show it
if b_plot:
    plt.show(block=False)
