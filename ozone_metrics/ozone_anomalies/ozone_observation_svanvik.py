plt.close('all')

fig1 = plt.figure(1)
fig1.canvas.set_window_title("ozone_obs_2018_07")
ax11 = plt.subplot(211)
data_svanvik_OzoNorClim.loc['2018-07'].plot(ax=ax11, color='blueviolet', ls='None', marker='x', label='Svanvik (2018)')
data['Esrange'].loc['2018-07'].plot(ls='None', marker='o', color='blue', label="Esrange (2018)")
data['Pallas'].loc['2018-07'].plot(ls='None', marker='^', color='black', label="Pallas (2018)")
ax11.plot(pd.date_range("2018-07","2018-07-31 23:00", freq='H'), clim_hourly.loc[7].values, color='red', label="Hourly clim.")
ax11.plot(pd.date_range("2018-07","2018-07-31 23:00", freq='H'), clim_hourly_svanvik.loc[7].values, color='magenta', label="Hourly clim. Svanvik")

plt.plot((scaling_max_pallas*sample_clim_hourly_svanvik[0]['2018-07']+sample_clim_hourly_svanvik[0]['2018-07']), color='grey', label='corr. Svanvik')
data_svanvik_rra.plot(color='orange', label='rra Svanvik')

ax11.set_ylabel("$[O_3] (ppb)$")
ax11.set_ylim(0,80)
ax11.legend(ncol=2)

ax12 = plt.subplot(212)
ax12.hist((data['Pallas'].loc['2018-07'].values-clim_hourly.loc[7].values.flatten())/clim_hourly_err.loc[7].values.flatten(), density=True, histtype='step', label='Pallas 2018')
ax12.hist((data['Esrange'].loc['2018-07'].values-clim_hourly.loc[7].values.flatten())/clim_hourly_err.loc[7].values.flatten(), density=True, histtype='step', color='blue', label='Esrange 2018')
ax12.hist((data_svanvik_OzoNorClim.loc['2018-07']-sample_clim_hourly_svanvik[0][data_svanvik_OzoNorClim.loc['2018-07'].index])/sample_clim_hourly_err_svanvik[0][data_svanvik_OzoNorClim.loc['2018-07'].index], density=True, histtype='step', color='violet', label='Svanvik 2018')
ax12.hist((data_svanvik_rra.loc['2018-07']-sample_clim_hourly_svanvik[0][data_svanvik_rra.time.values])/sample_clim_hourly_err_svanvik[0][data_svanvik_rra.time.values], density=True, histtype='step', color='orange', label='Svanvik 2018 rra')


ax12.plot(x_test_pallas_hourly_2018, pdf_test_pallas_hourly_2018, color='black', ls='-', label='Skew normal fit: Pallas')
stats_text(ax12, stat_test_pallas_hourly_2018, fit_test_pallas_hourly_2018, name="Pallas 2018", ypos=0.7)

ax12.plot(x_test_esrange_hourly_2018, pdf_test_esrange_hourly_2018, color='black', ls='--', label='Skew normal fit: Esrange')
stats_text(ax12, stat_test_esrange_hourly_2018, fit_test_esrange_hourly_2018, name="Esrange 2018", ypos=0.4)

ax12.plot(x_test_svanvik_hourly_2018, pdf_test_svanvik_hourly_2018, color='black', ls='-.', label='Skew normal fit: Svanvik')
stats_text(ax12, stat_test_svanvik_hourly_2018, fit_test_svanvik_hourly_2018, name="Svanvik 2018", ypos=0.1)

ax12.plot(x_test_svanvik_rra_hourly_2018, pdf_test_svanvik_rra_hourly_2018, color='black', ls=':', label='Skew normal fit: Svanvik_Rra')
stats_text(ax12, stat_test_svanvik_rra_hourly_2018, fit_test_svanvik_rra_hourly_2018, name="Svanvik rra 2018", ypos=0.7, xpos=0.18)

ax12.set_ylabel("Probability density")
ax12.set_xlabel("$[O_3]_{2018}-[O_3]_{clim}/\sigma [O_3]_{clim}$")
ax12.legend(loc='upper right')


# Monthly anomalies
data_svanvik = data_svanvik_OzoNorClim.copy()
data_svanvik.index = pd.to_datetime(data_svanvik.index.values.astype(str))
data_svanvik.loc[:,'day'] = data_svanvik.index.day.values
data_svanvik.loc[:,'month'] = data_svanvik.index.month.values
data_svanvik.loc[:,'hour'] = data_svanvik.index.hour.values

data_svanvik_clim = data['Svanvik'].copy()
data_svanvik_clim.loc[:,'day'] = data_svanvik_clim.index.day.values
data_svanvik_clim.loc[:,'month'] = data_svanvik_clim.index.month.values
data_svanvik_clim.loc[:,'hour'] = data_svanvik_clim.index.hour.values

svanvik_ozone_clim = data_svanvik_clim.groupby(['month','day','hour']).mean().iloc[:,0]
svanvik_ozone_clim_std = data_svanvik_clim.groupby(['month','day','hour']).std().iloc[:,0]

fig3 = plt.figure(3)
fig3.canvas.set_window_title("ozone_signific")
ax31 = plt.subplot(211)
ax31.set_title("(a)")
ax32 = plt.subplot(212)
ax32.set_title("(b)")
for iyear, iax, icolor in zip((2018, 2019),(ax31, ax32), ('violet', 'purple')):
    tmp = (data_svanvik['%d' % iyear].dropna()).groupby(['month','day','hour']).mean()
    test = ((tmp-svanvik_ozone_clim)/svanvik_ozone_clim_std)
    (test.where(test>1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100).plot.bar(ax=iax, color=icolor)
    (test.where(test<-1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*-100).plot.bar(ax=iax, color=icolor)
    iax.set_xlabel("")
    iax.set_ylabel("")
    iax.set_ylim(-60, 60)
    iax.axhline(0, color='grey', ls=':')
    iax.axhline(15.9, color='grey', ls='--')
    iax.axhline(-15.9, color='grey', ls='--')
    
    iax.text(9.25,53, '$>+1\,\sigma$: %3.2f %s' % (test.where(test>1).dropna().size/float(test.size)*100, "%"), size='x-large')
    iax.text(9.5,-55, '$<-1\,\sigma$: %3.2f %s' % (test.where(test<-1).dropna().size/float(test.size)*100, "%"), size='x-large')
    iax.tick_params(labelrotation=0)
    iax.set_xticklabels([get_month_name(imonth, length=3) for imonth in np.arange(1,13)])
    
    print('%d %3.2f' % (iyear, test.where(test>1).dropna().size/float(test.size)*100))
    print('%d %3.2f' % (iyear, test.where(test<-1).dropna().size/float(test.size)*100))
   
    print('%d %s' % (iyear, test.where(test>1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    print('%d %s' % (iyear, test.where(test<-1).dropna().groupby(['month']).apply(np.size).astype(float)/(test.groupby(['month']).apply(np.size))*100))
    
    
ax32.set_xlabel("Time (months)")
ax32.set_ylabel("#Days above $\pm 1\sigma_{clim}$ (%)", y=1)


plt.show(block=False)
