# Clean up
plt.close('all')

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)

clim_karasjok_o3 = data_finnmark.groupby('time.hour').mean()
climerr_karasjok_o3 = data_finnmark.groupby('time.hour').std()/np.sqrt(data_finnmark.size)

clim_karasjok_o3_winter = data_finnmark_winter.groupby('time.hour').mean()
climerr_karasjok_o3_winter = data_finnmark_winter.groupby('time.hour').std()/np.sqrt(data_finnmark_winter.size)

clim_karasjok_o3_summer = data_finnmark_summer.groupby('time.hour').mean()
climerr_karasjok_o3_summer = data_finnmark_summer.groupby('time.hour').std()/np.sqrt(data_finnmark_summer.size)

obs_karasjok_winter = data_karasjok_winter.groupby(data_karasjok_winter.index.hour).apply(np.nanmean)
obserr_karasjok_winter = data_karasjok_winter.groupby(data_karasjok_winter.index.hour).apply(np.nanstd)/np.sqrt(data_karasjok_winter.size)

obs_karasjok_summer = data_karasjok_summer.groupby(data_karasjok_summer.index.hour).apply(np.nanmean)
obserr_karasjok_summer = data_karasjok_summer.groupby(data_karasjok_summer.index.hour).apply(np.nanstd)/np.sqrt(data_karasjok_summer.size)

# OsloCTM3 v1.0
#clim_karasjok_o3.plot(ax=ax31, ls="None", marker='o', fillstyle='none',color='black', label='OsloCTM3 v1.0')
#plot_error_bands(ax11,
#                 clim_karasjok_o3.hour,
#                 clim_karasjok_o3.data,
#                 climerr_karasjok_o3)

# Split into summer and winter
clim_karasjok_o3_winter.plot(ax=ax11, ls="None", marker='o', color='black', alpha=0.5, label='OsloCTM3 v1.0')
clim_karasjok_o3_summer.plot(ax=ax12, ls="None", marker='o', color='black', alpha=0.5, label='OsloCTM3 v1.0 ')

plot_error_bands(ax11,
                 clim_karasjok_o3_winter.hour,
                 clim_karasjok_o3_winter.data,
                 climerr_karasjok_o3_winter)

clim_karasjok_o3_summer.plot(ax=ax12)
plot_error_bands(ax12,
                 clim_karasjok_o3_summer.hour,
                 clim_karasjok_o3_summer.data,
                 climerr_karasjok_o3_summer)

# Split into summer and winter
#obs_karasjok_winter.plot(ax=ax11, ls="--", color='orange', alpha=0.5, label='Karasjok')
#obs_karasjok_summer.plot(ax=ax12, ls="--", color='orange', alpha=0.5, label='Karasjok')

#plot_error_bands(ax11,
#                 obs_karasjok_winter.hour,
#                 obs_karasjok_winter.data,
#                 obserr_karasjok_winterclimerr_karasjok_o3_winter)

#clim_karasjok_o3_summer.plot(ax=ax12)
#plot_error_bands(ax12,
#                 clim_karasjok_o3_summer.hour,
#                 clim_karasjok_o3_summer.data,
#                 climerr_karasjok_o3_summer)

# Decorate the plots
for ax in fig1.axes:
    ax.set_ylim(20,45)
    ax.set_xlabel("")
    ax.set_ylabel("[$O_3$] (ppb)")
    
ax11.set_title("Winter")
ax12.set_title("Summer")
ax.set_xlabel("Time (hour)")
# Show it
plt.show(block=False)
