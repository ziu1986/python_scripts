# Clean up
plt.close('all')

from mpl_toolkits.mplot3d import Axes3D # Register 3d projection

fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("svanvik_krekling_conductance")
ax1 = plt.subplot()
ax1.hist(g_sto_o3, bins=100, density=True, color='orange', label='all')
g_sto_o3_month['Jun'].hist(ax=ax1, bins=100, histtype='step', density=True, label='Jun')
g_sto_o3_month['Aug'].hist(ax=ax1, bins=100, histtype='step', density=True, color='red', label='Aug' )
g_sto_o3_month['Sep'].hist(ax=ax1, bins=100, histtype='step', density=True, color='blue', label='Sep')


ax1.plot(x_sample, pdf, label="fit all months", color='orange')
stats_text(ax1, stat, fit, name="Cond. all months", ypos=0.7)

#for imonth, icolor in zip(('Jun','Aug', 'Sep'), ('black', 'red', 'blue')):
#    ax1.plot(fit_results_month[imonth]['x_sample'], fit_results_month[imonth]['pdf'], label="fit %s" % imonth, color=icolor)
#stats_text(ax1, stat, fit, name="g_sto", ypos=0.7)

ax1.set_xlabel("$g_{sto}$ $(mmol\,m^{-2}\,s^{-1}\,g^{-1})$")
ax1.set_ylabel("Probability Density")

ax1.legend(ncol=2)

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("svanvik_krekling-contuctance_photosyth")
ax21 = plt.subplot(221)
ax21.plot(data_krekling['HHMMSS'], (data_krekling['June Photo']), ls='None', marker='x', label='Jun')
ax21.plot(data_krekling['HHMMSS.1'], (data_krekling['Aug Photo']), ls='None', marker='+', color='red', label='Aug')
ax21.plot(data_krekling['HHMMSS.2'], (data_krekling['Sept Photo']), ls='None', marker='.', color='blue', label='Sep')

ax21.plot(data_krekling['HHMMSS'], A_net_o3_month['Jun'], ls='None', marker='o', fillstyle='none', markersize=10, label='_')
ax21.plot(data_krekling['HHMMSS.1'], A_net_o3_month['Aug'], ls='None', marker='o', fillstyle='none', markersize=10, color='red', label='_')
ax21.plot(data_krekling['HHMMSS.2'], A_net_o3_month['Sep'], ls='None', marker='o', fillstyle='none', markersize=10, color='blue', label='_')


ax21.set_ylabel("$A_{net}$ $(mol\,m^{-2}\,s^{-1}\,g^{-1})$")
ax21.set_xlabel("Time")

ax22 = plt.subplot(222)
ax22.plot(data_krekling['HHMMSS'], (1e3*k_O3*data_krekling['June Cond']), ls='None', marker='x', label='Jun')
ax22.plot(data_krekling['HHMMSS.1'], (1e3*k_O3*data_krekling['Aug Cond']), ls='None', marker='+', color='red', label='Aug')
ax22.plot(data_krekling['HHMMSS.2'], (1e3*k_O3*data_krekling['Sept Cond']), ls='None', marker='.', color='blue', label='Sep')

ax22.plot(data_krekling['HHMMSS'], g_sto_o3_month['Jun'], ls='None', marker='o', fillstyle='none', markersize=10, label='_')
ax22.plot(data_krekling['HHMMSS.1'], g_sto_o3_month['Aug'], ls='None', marker='o', fillstyle='none', markersize=10, color='red', label='_')
ax22.plot(data_krekling['HHMMSS.2'], g_sto_o3_month['Sep'], ls='None', marker='o', fillstyle='none', markersize=10, color='blue', label='_')


ax22.set_ylabel("$g_{sto}$ $(mmol\,m^{-2}\,s^{-1}\,g^{-1})$")
ax22.set_xlabel("Time")
ax22.legend()

ax23 = plt.subplot(223)
ax23.plot(np.log(data_krekling['PARo']), (data_krekling['June Photo']), ls='None', marker='x', label='Jun')
ax23.plot(np.log(data_krekling['PARo.1']), (data_krekling['Aug Photo']), ls='None', marker='+', color='red', label='Aug')
ax23.plot(np.log(data_krekling['PARo.2']), (data_krekling['Sept Photo']), ls='None', marker='.', color='blue', label='Sep')

ax23.set_ylabel("$A_{net}$ $(mol\,m^{-2}\,s^{-1}\,g^{-1})$")
#ax23.set_xlabel("PAR $(\mu\,mol\,m^{-2}\,s^{-1})$")
ax23.set_xlabel("log(PAR)")

ax24 = plt.subplot(224)
ax24.plot(np.log(data_krekling['PARo']), (1e3*k_O3*data_krekling['June Cond']), ls='None', marker='x', label='Jun')
ax24.plot(np.log(data_krekling['PARo.1']), (1e3*k_O3*data_krekling['Aug Cond']), ls='None', marker='+', color='red', label='Aug')
ax24.plot(np.log(data_krekling['PARo.2']), (1e3*k_O3*data_krekling['Sept Cond']), ls='None', marker='.', color='blue', label='Sep')

ax24.set_ylabel("$g_{sto}$ $(mmol\,m^{-2}\,s^{-1}\,g^{-1})$")
#ax24.set_xlabel("PAR $(\mu\,mol\,m^{-2}\,s^{-1})$")
ax24.set_xlabel("log(PAR)")
ax24.legend()

fig3 = plt.figure(3)
fig3.canvas.set_window_title("svanvik_krekling_3d_ppfd_temp")
ax31 = plt.subplot(projection='3d')
ax31.scatter((data_krekling['PARo']),
             (data_krekling['Tleaf (air)']),
             g_sto_o3_month['Jun']/gmax,
             marker='x', label='Jun')
ax31.scatter((data_krekling['PARo.1']),
             (data_krekling['Tleaf (air).1']),
             g_sto_o3_month['Aug']/gmax,
             marker='+', color='red', label='Aug')
ax31.scatter((data_krekling['PARo.2']),
             (data_krekling['Tleaf (air).2']),
             g_sto_o3_month['Sep']/gmax,
             marker='o', color='blue', label='Sep')

ax31.set_xlabel("PAR $(\mu\,mol\,m^{-2}\,s^{-1})$")
ax31.set_ylabel("T $(^\circ C)$")
ax31.set_zlabel("$g_{sto}^{rel}$")

ax31.legend()

fig4 = plt.figure(4)
fig4.canvas.set_window_title("svanvik_krekling_3d_ppfd_vpd")
ax41 = plt.subplot(projection='3d')
ax41.scatter((data_krekling['PARo']),
             VPD(data_krekling['RH_R'], data_krekling['Tleaf (air)'], version='buck')/kilo,
             g_sto_o3_month['Jun']/gmax,
             marker='x', label='Jun')
ax41.scatter((data_krekling['PARo.1']),
             VPD(data_krekling['RH_R.1'], data_krekling['Tleaf (air).1'], version='buck')/kilo,
             g_sto_o3_month['Aug']/gmax,
             marker='+', color='red', label='Aug')
ax41.scatter((data_krekling['PARo.2']),
             VPD(data_krekling['RH_R.2'], data_krekling['Tleaf (air).2'], version='buck')/kilo,
             g_sto_o3_month['Sep']/gmax,
             marker='o', color='blue', label='Sep')

ax41.set_xlabel("PAR $(\mu\,mol\,m^{-2}\,s^{-1})$")
ax41.set_ylabel("VPD $(kPa)$")
ax41.set_zlabel("$g_{sto}^{rel}$")

ax41.legend()

fig5 = plt.figure(5)
fig5.canvas.set_window_title("svanvik_krekling_3d_ppfd_relhum")
ax51 = plt.subplot(projection='3d')
ax51.scatter((data_krekling['PARo']),
             data_krekling['RH_R'],
             g_sto_o3_month['Jun']/gmax,
             marker='x', label='Jun')
ax51.scatter((data_krekling['PARo.1']),
             data_krekling['RH_R.1'],
             g_sto_o3_month['Aug']/gmax,
             marker='+', color='red', label='Aug')
ax51.scatter((data_krekling['PARo.2']),
             data_krekling['RH_R.2'],
             g_sto_o3_month['Sep']/gmax,
             marker='o', color='blue', label='Sep')

ax51.set_xlabel("PAR $(\mu\,mol\,m^{-2}\,s^{-1})$")
ax51.set_ylabel("RelHum $(\%)$")
ax51.set_zlabel("$g_{sto}^{rel}$")

ax51.legend()

##

fig6 = plt.figure(6)
fig6.canvas.set_window_title("svanvik_krekling_f_light")
ax61 = plt.subplot()
ax61.plot((data_krekling['PARo']),
          g_sto_o3_month['Jun']/gmax,
          ls='None', marker='x', label='Jun')
ax61.plot((data_krekling['PARo.1']),
          g_sto_o3_month['Aug']/gmax,
          ls='None', marker='+', color='red', label='Aug')
ax61.plot((data_krekling['PARo.2']),
          g_sto_o3_month['Sep']/gmax,
          ls='None', marker='.', color='blue', label='Sep')

# Plot fit
from decimal import Decimal

sample_ppfd = np.arange(0,2000)
ax61.plot(sample_ppfd, f_light(sample_ppfd, fit_params[0]), color='orange', ls='--', lw=5)
ax61.text(0.5, 0.9, "Fit\n alpha = (%0.2E +/- %0.2E) $m^2 s\,\mu mol^{-1}$" % (Decimal(fit_params[0]), Decimal(cov_mat[0][0])), size='large', color='orange', transform=ax61.transAxes)

ax61.set_xlabel("PPFD $(\mu\,mol\,m^{-2}\,s^{-1})$")
ax61.set_ylabel("$g_{sto}^{rel}$")
ax61.legend()


fig7 = plt.figure(7)
fig7.canvas.set_window_title("svanvik_krekling_f_temp")
ax71 = plt.subplot()
ax71.plot((data_krekling['Tleaf (air)']),
          g_sto_o3_month['Jun']/gmax,
          ls='None', marker='x', label='Jun')
ax71.plot((data_krekling['Tleaf (air).1']),
          g_sto_o3_month['Aug']/gmax,
          ls='None', marker='+', color='red', label='Aug')
ax71.plot((data_krekling['Tleaf (air).2']),
          g_sto_o3_month['Sep']/gmax,
          ls='None', marker='.', color='blue', label='Sep')

sample_temp = np.arange(0,50)
ax71.plot(sample_temp, f_temp((sample_temp, fmin), 12., 0., 25.), color='grey', ls='--', lw=5)

ax71.plot(sample_temp, f_temp((sample_temp, fmin), fit_params_temp[0], fit_params_temp[1], fit_params_temp[2]), color='orange', ls='--', lw=5)
ax71.text(0.5, 0.88, "Fit\n T_opt = (%0.2E+/-%0.2E) deg C\n T_min = (%0.2E+/-%0.2E) deg C\n T_max = (%0.2E+/-%0.2E) deg C" % (Decimal(fit_params_temp[0]),Decimal(cov_mat_temp[0][0]), Decimal(fit_params_temp[1]),Decimal(cov_mat_temp[1][1]), Decimal(fit_params_temp[2]),Decimal(cov_mat_temp[2][2])), size='large', color='orange', transform=ax71.transAxes)

ax71.text(0.5, 0.76, "Guess\n T_opt = (%d) deg C\n T_min = (%d) deg C\n T_max = (%d) deg C" % (12, 0, 25), size='large', color='grey', transform=ax71.transAxes)


ax71.set_xlabel("T $(^\circ C)$")
ax71.set_ylabel("$g_{sto}^{rel}$")
ax71.legend()

ax71.set_xlim(0, 50)


fig8 = plt.figure(8)
fig8.canvas.set_window_title("svanvik_krekling_f_vpd")
ax81 = plt.subplot()
ax81.plot(VPD(data_krekling['RH_R'], data_krekling['Tleaf (air)'], version='buck')/kilo,
          g_sto_o3_month['Jun']/gmax,
           ls='None', marker='x', label='Jun')
ax81.plot(VPD(data_krekling['RH_R.1'], data_krekling['Tleaf (air).1'], version='buck')/kilo,
          g_sto_o3_month['Aug']/gmax,
          ls='None', marker='+', color='red', label='Aug')
ax81.plot(VPD(data_krekling['RH_R.2'], data_krekling['Tleaf (air).2'], version='buck')/kilo,
          g_sto_o3_month['Sep']/gmax,
          ls='None', marker='.', color='blue', label='Sep')

sample_vpd = np.arange(0,6)
ax81.plot(sample_vpd, f_vpd((sample_vpd, fmin), 0.5, 2.), color='grey', ls='--', lw=5)

ax81.plot(sample_vpd, f_vpd((sample_vpd, fmin), fit_params_vpd[0], fit_params_vpd[1]), color='orange', ls='--', lw=5)
ax81.text(0.5, 0.88, "Fit\n VPD_max = (%0.2E+/-%0.2E) kPa\n VPD_min = (%0.2E+/-%0.2E) kPa" % (Decimal(fit_params_vpd[0]),Decimal(cov_mat_vpd[0][0]), Decimal(fit_params_vpd[1]),Decimal(cov_mat_vpd[1][1])), size='large', color='orange', transform=ax81.transAxes)

ax81.text(0.5, 0.76, "Guess\n VPD_max = (%1.1f) kPa\n VPD_min = (%1.1f) kPa" % (0.5, 2), size='large', color='grey', transform=ax81.transAxes)


ax81.set_xlabel("VPD $(kPa)$")
ax81.set_ylabel("$g_{sto}^{rel}$")
ax81.legend()

ax81.set_xlim(0, 5)

fig9 = plt.figure(9)
ax91 = plt.subplot()

ax91.hist(temp-dew_point(relHum, temp))
ax91.set_xlabel("$T_{dp} - T_{leaf}$")
ax91.set_ylabel("Counts")

# Show it
plt.show(block=False)
