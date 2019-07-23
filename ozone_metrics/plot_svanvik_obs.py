# To be executed after plot_ozone_observations.py

clim_o3 = data['Svanvik'].groupby(data['Svanvik'].index.dayofyear).mean()
climerr_o3 = data['Svanvik'].groupby(data['Svanvik'].index.dayofyear).std()/np.sqrt(data['Svanvik'].groupby(data['Svanvik'].index.dayofyear).sum()/clim_o3)

# Plotting
plt.close('all')
fig1 = plt.figure(1,figsize=(16,3))
fig1.canvas.set_window_title("ebas_climatology_svanvik")

ax11 = plt.subplot()
clim_o3.plot(ax=ax11,
             yerr=climerr_o3,
             marker='.', ls='none', color='grey', label='Svanvik (1986-1996)')
data_svanvik_2018.groupby(data_svanvik_2018.index.dayofyear).mean().plot(ax=ax11, ls='None', marker='x', label='Svanvik (2018)', color='blueviolet')

ax11.set_title('')
ax11.set_ylabel('%s (ppb)' % ('O3'))
ax11.set_xlabel("dayofyear")
ax11.set_xlim(0,367)
ax11.set_ylim(0,60)
ax11.legend()
plot_month_span(ax11)
plot_month_name(ax11, 55)

# Show it
plt.show(block=False)
