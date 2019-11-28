# Execute plot_ozone_observations.py first
plt.close('all')
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(311)
ax12 = plt.subplot(312)
ax13 = plt.subplot(313)
spli_deriv = []
spli_deriv2 = []
for smoothing in np.linspace(0,1000):
    fitSpl_dmean = UnivariateSpline(np.unique(doys), climatology.groupby(climatology.index.dayofyear).apply(np.nanmean), w=w)
    fitSpl_dmean.set_smoothing_factor(smoothing)
    ax11.plot(np.linspace(1,367), fitSpl_dmean(np.linspace(1,367)))
    spli_deriv.append(fitSpl_dmean.derivative(1)(np.linspace(1,367))[0]-fitSpl_dmean.derivative(1)(np.linspace(1,367))[-1])
    spli_deriv2.append(fitSpl_dmean.derivative(2)(np.linspace(1,367))[0]-fitSpl_dmean.derivative(1)(np.linspace(1,367))[-1])
    print(smoothing, len(fitSpl_dmean.get_knots()), fitSpl_dmean.derivative(1)(np.linspace(1,367))[0]-fitSpl_dmean.derivative(1)(np.linspace(1,367))[-1])


ax12.plot(np.linspace(0,1000), spli_deriv, label='derivative_1')
ax12.plot(np.linspace(0,1000), spli_deriv2, ls='--', label='derivative_2')
minimum = np.linspace(0,1000)[np.where((np.array(spli_deriv)<=0.01) & (np.array(spli_deriv)>=-0.01))[0]]
ax12.axvline(minimum, color='black', ls=':')

fitSpl_dmean = UnivariateSpline(np.unique(doys), climatology.groupby(climatology.index.dayofyear).apply(np.nanmean), w=w)
fitSpl_dmean.set_smoothing_factor(minimum)
ax13.plot(np.linspace(1,367), fitSpl_dmean(np.linspace(1,367)))
ax13.errorbar(xtime[::10], yozone[::10], yerr=yerr[::10], color='black', marker='o', ls='None', label='monthly mean (1 $\sigma$)')
ax13.errorbar(xtime[::10], yozone[::10], yerr=yerr_mean[::10], color='blue', marker='None', ls='None', label='std_err_mean')

ax11.set_xlabel("Day of year")
ax11.set_ylabel("$[O_3]$ (ppb)")
ax12.set_xlabel("Smoothing factor")
ax12.set_ylabel("Deviation of endpoints")
ax12.legend()
ax13.set_xlabel("Day of year")
ax13.set_ylabel("$[O_3]$ (ppb)")

fig2 = plt.figure(2, figsize=(16,9))

plt.show(block=False)
