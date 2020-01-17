plt.close('all')

fig1 = plt.figure(1)
fig1.canvas.set_window_title("ozone_obs_2018")
ax11 = plt.subplot(211)
data_svanvik_OzoNorClim.loc['2018-07'].plot(ax=ax11, color='blueviolet', ls='None', marker='x', label='Svanvik (2018)')
data['Esrange'].loc['2018-07'].plot(ls='None', marker='o', color='blue', label="Esrange (2018)")
data['Pallas'].loc['2018-07'].plot(ls='None', marker='^', color='black', label="Pallas (2018)")
ax11.plot(pd.date_range("2018-07","2018-07-31 23:00", freq='H'), clim_hourly.loc[7].values, color='red', label="Hourly clim.")
ax11.plot(pd.date_range("2018-07","2018-07-31 23:00", freq='H'), clim_hourly_svanvik.loc[7].values, color='magenta', label="Hourly clim. Svanvik")

plt.plot((scaling_max_pallas*sample_clim_hourly_svanvik[0]['2018-07']+sample_clim_hourly_svanvik[0]['2018-07']), color='grey', label='corr. Svanvik')

ax11.set_ylabel("$[O_3] (ppb)$")
ax11.set_ylim(0,80)
ax11.legend(ncol=2)

ax12 = plt.subplot(212)
ax12.hist((data['Pallas'].loc['2018-07'].values-clim_hourly.loc[7].values.flatten())/clim_hourly_err.loc[7].values.flatten(), density=True, histtype='step', label='Pallas 2018')
ax12.hist((data['Esrange'].loc['2018-07'].values-clim_hourly.loc[7].values.flatten())/clim_hourly_err.loc[7].values.flatten(), density=True, histtype='step', color='blue', label='Esrange 2018')
ax12.hist((data_svanvik_OzoNorClim.loc['2018-07']-sample_clim_hourly_svanvik[0][data_svanvik_OzoNorClim.loc['2018-07'].index])/sample_clim_hourly_err_svanvik[0][data_svanvik_OzoNorClim.loc['2018-07'].index], density=True, histtype='step', color='violet', label='Svanvik 2018')


ax12.plot(x_test_pallas_hourly_2018, pdf_test_pallas_hourly_2018, color='black', ls='-', label='Skew normal fit: Pallas')
stats_text(ax12, stat_test_pallas_hourly_2018, fit_test_pallas_hourly_2018, name="Pallas 2018", ypos=0.7)

ax12.plot(x_test_esrange_hourly_2018, pdf_test_esrange_hourly_2018, color='black', ls='--', label='Skew normal fit: Esrange')
stats_text(ax12, stat_test_esrange_hourly_2018, fit_test_esrange_hourly_2018, name="Esrange 2018", ypos=0.4)

ax12.plot(x_test_svanvik_hourly_2018, pdf_test_svanvik_hourly_2018, color='black', ls='-.', label='Skew normal fit: Svanvik')
stats_text(ax12, stat_test_svanvik_hourly_2018, fit_test_svanvik_hourly_2018, name="Svanvik 2018", ypos=0.1)

ax12.set_ylabel("Probability density")
ax12.set_xlabel("$[O_3]_{2018}-[O_3]_{clim}/\sigma [O_3]_{clim}$")
ax12.legend(loc='upper right')





plt.show(block=False)
