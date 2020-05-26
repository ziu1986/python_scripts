# Clean up
plt.close('all')

fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("svanvik_krekling_conductance")
ax1 = plt.subplot()
(1e3*k_O3*data_krekling['June Cond']).where(np.log(data_krekling['PARo'])>=4).hist(ax=ax1, bins=100, density=True, label='Jun')
(1e3*k_O3*data_krekling['Aug Cond']).where(np.log(data_krekling['PARo.1'])>=4).hist(ax=ax1, bins=100, histtype='step', density=True, color='red', label='Aug' )
(1e3*k_O3*data_krekling['Sept Cond']).where(np.log(data_krekling['PARo.2'])>=4).hist(ax=ax1, bins=100, histtype='step', density=True, color='blue', label='Sep')

ax1.plot(x_sample, pdf, label="fit")
stats_text(ax1, stat, fit, name="Conductance", ypos=0.7)

ax1.set_xlabel("$g_{sto}$ $(mmol\,m^{-2}\,s^{-1}\,g^{-1})$")
ax1.set_ylabel("Probability Density")

ax1.legend()

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("svanvik_krekling-contuctance_photosyth")
ax21 = plt.subplot(221)
ax21.plot(data_krekling['HHMMSS'], (data_krekling['June Photo']), ls='None', marker='x', label='Jun')
ax21.plot(data_krekling['HHMMSS.1'], (data_krekling['Aug Photo']), ls='None', marker='+', color='red', label='Aug')
ax21.plot(data_krekling['HHMMSS.2'], (data_krekling['Sept Photo']), ls='None', marker='.', color='blue', label='Sep')

ax21.set_ylabel("$A_{net}$ $(mol\,m^{-2}\,s^{-1}\,g^{-1})$")
ax21.set_xlabel("Time")

ax22 = plt.subplot(222)
ax22.plot(data_krekling['HHMMSS'], (1e3*k_O3*data_krekling['June Cond']), ls='None', marker='x', label='Jun')
ax22.plot(data_krekling['HHMMSS.1'], (1e3*k_O3*data_krekling['Aug Cond']), ls='None', marker='+', color='red', label='Aug')
ax22.plot(data_krekling['HHMMSS.2'], (1e3*k_O3*data_krekling['Sept Cond']), ls='None', marker='.', color='blue', label='Sep')

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

#fig3 = plt.figure(3)
#ax31 = plt.subplot()

#ax31.plot(data_krekling['HHMMSS'], data_krekling['PARo'], ls='None', marker='x', label='Jun')
#ax31.plot(data_krekling['HHMMSS.1'], data_krekling['PARo.1'], ls='None', marker='+', color='red', label='Aug')
#ax31.plot(data_krekling['HHMMSS.2'], data_krekling['PARo.2'], ls='None', marker='o', color='blue', label='Sept')

#ax31.set_xlabel("Time (hours)")
#ax31.set_ylabel("PAR $(\mu\,mol\,m^{-2}\,s^{-1})$")

##

fig6 = plt.figure(6)
fig6.canvas.set_window_title("svanvik_krekling_f_light")
ax61 = plt.subplot()
ax61.plot((data_krekling['PARo']), (1e3*k_O3*data_krekling['June Cond'].where((np.log(data_krekling['PARo'])>photosynth_loglimit)))/gmax, ls='None', marker='x', label='Jun')
ax61.plot((data_krekling['PARo.1']), (1e3*k_O3*data_krekling['Aug Cond'].where((np.log(data_krekling['PARo.1'])>photosynth_loglimit)))/gmax, ls='None', marker='+', color='red', label='Aug')
ax61.plot((data_krekling['PARo.2']), (1e3*k_O3*data_krekling['Sept Cond'].where((np.log(data_krekling['PARo.2'])>photosynth_loglimit)))/gmax, ls='None', marker='.', color='blue', label='Sep')

# Plot fit
from decimal import Decimal

sample_ppfd = np.arange(0,2000)
ax61.plot(sample_ppfd, f_light(sample_ppfd, fit_params[0]), color='orange', ls='--', lw=5)
ax61.text(0.5, 0.9, "Fit\n alpha = %0.2E+/-%0.2E" % (Decimal(fit_params[0]),Decimal(cov_mat[0][0])), size='large', color='orange', transform=ax61.transAxes)

ax61.set_xlabel("PPFD $(\mu\,mol\,m^{-2}\,s^{-1})$")
ax61.set_ylabel("$g_{sto}^{rel}$")
ax61.legend()


fig7 = plt.figure(7)
fig7.canvas.set_window_title("svanvik_krekling_f_temp")
ax71 = plt.subplot()
ax71.plot((data_krekling['Tleaf (air)']), (1e3*k_O3*data_krekling['June Cond'].where((np.log(data_krekling['PARo'])>photosynth_loglimit)))/gmax, ls='None', marker='x', label='Jun')
ax71.plot((data_krekling['Tleaf (air).1']), (1e3*k_O3*data_krekling['Aug Cond'].where((np.log(data_krekling['PARo.1'])>photosynth_loglimit)))/gmax, ls='None', marker='+', color='red', label='Aug')
ax71.plot((data_krekling['Tleaf (air).2']), (1e3*k_O3*data_krekling['Sept Cond'].where((np.log(data_krekling['PARo.2'])>photosynth_loglimit)))/gmax, ls='None', marker='.', color='blue', label='Sep')

sample_temp = np.arange(0,50)
ax71.plot(sample_temp, f_temp(sample_temp, 12, 0, 25), color='grey', ls='--', lw=5)

ax71.plot(sample_temp, f_temp(sample_temp, fit_params_temp[0], fit_params_temp[1], fit_params_temp[2]), color='orange', ls='--', lw=5)
ax61.text(0.5, 0.9, "Fit\n alpha = %0.2E+/-%0.2E" % (Decimal(fit_params[0]),Decimal(cov_mat[0][0])), size='large', color='orange', transform=ax61.transAxes)


ax71.set_xlabel("T $(^\circ C)$")
ax71.set_ylabel("$g_{sto}^{rel}$")
ax71.legend()

ax71.set_xlim(0, 50)

# Show it
plt.show(block=False)
