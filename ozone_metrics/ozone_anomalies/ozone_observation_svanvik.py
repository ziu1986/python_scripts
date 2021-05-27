plt.close('all')

# Bias correction for historical climatology to present day
bias_corr = 1.2
# Time lag corection (same for Esrange and Pallas)
time_lag_corr = lag_max['svanvik_esrange']

# Scaling factor
scaling = sample_clim_hourly_svanvik/sample_clim_hourly.shift(-time_lag_corr)
anomaly_pallas = data['Pallas']['07-2018']-sample_clim_hourly['07-2018'][0]
anomaly_esrange = data['Esrange']['07-2018']-sample_clim_hourly['07-2018'][0]
anomaly_svanvik = data_svanvik_OzoNorClim['2018-07']-sample_clim_hourly_svanvik[0][data_svanvik_OzoNorClim['2018-07'].index]-bias_corr

reco_anomaly_svanvik = anomaly_pallas.shift(-time_lag_corr)*scaling['07-2018'][0]
reco_svanvik = reco_anomaly_svanvik+sample_clim_hourly_svanvik['2018-07'][0]+bias_corr


fig1 = plt.figure(1, figsize=(10,12))
fig1.canvas.set_window_title("ozone_reconstruction_2018_07")
ax11 = plt.subplot(311)
ax11.set_title('(a)')

data['Esrange']['2018-07'].plot(ax=ax11, ls='None', marker='o', fillstyle='none', color='blue', label="Esrange")
data['Pallas']['2018-07'].plot(ax=ax11, ls='None', marker='^', fillstyle='none', color='black', label="Pallas")
data_svanvik_OzoNorClim['2018-07'].plot(ax=ax11, color='blueviolet', ls='None', marker='d', label='Svanvik')
sample_clim_hourly['07-2018'][0].plot(ax=ax11, color='red', label="Hourly clim.")
sample_clim_hourly.shift(-time_lag_corr)['2018-07'][0].plot(ax=ax11, color='red', ls='--', label="Hourly clim. + time lag corr.")
sample_clim_hourly_svanvik['07-2018'][0].plot(ax=ax11, color='red', ls='-.', label="Hourly clim. Svanvik")

ax11.set_ylabel("$[O_3] (ppb)$")
ax11.set_ylim(0,75)
ax11.set_xticklabels("")
ax11.set_xlabel('')
ax11.legend(ncol=2)
#
ax12 = plt.subplot(312)
ax12.set_title('(b)')
anomaly_pallas.plot(ax=ax12, ls='None', marker='^', fillstyle='none', color='black', label="Pallas")
anomaly_esrange.plot(ax=ax12, ls='None', marker='o', fillstyle='none', color='blue', label="Esrange")

anomaly_svanvik.plot(ax=ax12, ls='None', color='blueviolet', label='Svanvik', marker='d')
reco_anomaly_svanvik.plot(ax=ax12, color='magenta', label='Reco. Svanvik')


ax12.set_ylabel("$\Delta [O_3]$ (ppb)")
#ax12.set_xlabel("Time (days)")
ax12.set_ylim(-30, 30)
ax12.legend(ncol=2)

#fig2 = plt.figure(2, figsize=(16,9))
#fig2.canvas.set_window_title("ozone_reconstruction_svanvik")
ax13 = plt.subplot(313)
ax13.set_title('(c)')

reco_svanvik.plot(ax=ax13, ls='-', color='magenta', marker='*', label='Reco. Svanvik')
data_svanvik_OzoNorClim['2018-07'].plot(ax=ax13, color='blueviolet', fillstyle='none', ls='None', marker='d', label='Svanvik')
data_svanvik_rra.to_pandas().plot(ax=ax13, color='orange', fillstyle='none', ls=':', label='Regional Model Reanalysis')

ax13.set_ylabel("$[O_3] (ppb)$")
ax13.set_ylim(0,75)
ax13.set_xlabel('Time (days)')
ax13.legend(ncol=3)

print("RMSE reconstruction: %1.2f" % np.sqrt(((reco_svanvik-data_svanvik_OzoNorClim['2018-07'])**2).sum()/reco_svanvik.size))
print("RMSE regional reanalysis: %1.2f" % np.sqrt(((reco_svanvik-data_svanvik_rra.to_pandas())**2).sum()/reco_svanvik.size))

plt.show(block=False)

"""
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
"""

# Save reconstruction to csv file
save_data = pd.concat((data_svanvik_OzoNorClim['2018-01':'2018-07-09 9:00'], reco_svanvik['2018-07-09 9:00':'2018-07-24 8:00'].dropna(), data_svanvik_OzoNorClim['2018-07-24 9:00':]))

save_data.to_csv("svanvik_ozone_2018.csv")
data_svanvik_OzoNorClim['2019'].to_csv("svanvik_ozone_2019.csv")
(sample_clim_hourly_svanvik[0]+bias_corr).to_csv("svanvik_ozone_corr-climatology.csv")
