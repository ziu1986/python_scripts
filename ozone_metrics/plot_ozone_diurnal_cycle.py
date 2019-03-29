# Clean up
plt.close('all')

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(211)
ax12 = plt.subplot(212)

try:
    clim_karasjok_o3 = data_finnmark.groupby('time.hour').mean()
    climerr_karasjok_o3 = data_finnmark.groupby('time.hour').std()/np.sqrt(data_finnmark.size)
    
    clim_karasjok_o3_winter = data_finnmark_winter.groupby('time.hour').mean()
    climerr_karasjok_o3_winter = data_finnmark_winter.groupby('time.hour').std()/np.sqrt(data_finnmark_winter.size)

    clim_karasjok_o3_summer = data_finnmark_summer.groupby('time.hour').mean()
    climerr_karasjok_o3_summer = data_finnmark_summer.groupby('time.hour').std()/np.sqrt(data_finnmark_summer.size)
    b_model = True
except NameError:
    print("Model data not loaded!")
    b_model = False

try:
    data_jergul_winter = data_jergul.loc[(data_jergul.index.month<6) | (data_jergul.index.month>=9)]
    data_jergul_summer = data_jergul.loc[(data_jergul.index.month>=6) & (data_jergul.index.month<9)]
    data_karasjok_winter = data_karasjok.loc[(data_karasjok.index.month<6) | (data_karasjok.index.month>=9)]
    data_karasjok_summer = data_karasjok.loc[(data_karasjok.index.month>=6) & (data_karasjok.index.month<9)]
    data_svanvik_winter = data_svanvik.loc[(data_svanvik.index.month<6) | (data_svanvik.index.month>=9)]
    data_svanvik_summer = data_svanvik.loc[(data_svanvik.index.month>=6) & (data_svanvik.index.month<9)]
except NameError:
    data_jergul_winter = data['Jergul'].loc[(data['Jergul'].index.month<6) | (data['Jergul'].index.month>=9)]
    data_jergul_summer = data['Jergul'].loc[(data['Jergul'].index.month>=6) & (data['Jergul'].index.month<9)]
    data_karasjok_winter = data['Karasjok'].loc[(data['Karasjok'].index.month<6) | (data['Karasjok'].index.month>=9)]
    data_karasjok_summer = data['Karasjok'].loc[(data['Karasjok'].index.month>=6) & (data['Karasjok'].index.month<9)]
    data_svanvik_winter = data['Svanvik'].loc[(data['Svanvik'].index.month<6) | (data['Svanvik'].index.month>=9)]
    data_svanvik_summer = data['Svanvik'].loc[(data['Svanvik'].index.month>=6) & (data['Svanvik'].index.month<9)]

    obs_karasjok_winter = data_karasjok_winter.groupby(data_karasjok_winter.index.hour).apply(np.nanmean)
    obserr_karasjok_winter = data_karasjok_winter.groupby(data_karasjok_winter.index.hour).apply(np.nanstd)/np.sqrt(data_karasjok_winter.size)

    obs_karasjok_summer = data_karasjok_summer.groupby(data_karasjok_summer.index.hour).apply(np.nanmean)
    obserr_karasjok_summer = data_karasjok_summer.groupby(data_karasjok_summer.index.hour).apply(np.nanstd)/np.sqrt(data_karasjok_summer.size)

    obs_jergul_winter = data_jergul_winter.groupby(data_jergul_winter.index.hour).apply(np.nanmean)
    obserr_jergul_winter = data_jergul_winter.groupby(data_jergul_winter.index.hour).apply(np.nanstd)/np.sqrt(data_jergul_winter.size)

    obs_jergul_summer = data_jergul_summer.groupby(data_jergul_summer.index.hour).apply(np.nanmean)
    obserr_jergul_summer = data_jergul_summer.groupby(data_jergul_summer.index.hour).apply(np.nanstd)/np.sqrt(data_jergul_summer.size)

    obs_svanvik_winter = data_svanvik_winter.groupby(data_svanvik_winter.index.hour).apply(np.nanmean)
    obserr_svanvik_winter = data_svanvik_winter.groupby(data_svanvik_winter.index.hour).apply(np.nanstd)/np.sqrt(data_svanvik_winter.size)

    obs_svanvik_summer = data_svanvik_summer.groupby(data_svanvik_summer.index.hour).apply(np.nanmean)
    obserr_svanvik_summer = data_svanvik_summer.groupby(data_svanvik_summer.index.hour).apply(np.nanstd)/np.sqrt(data_svanvik_summer.size)

    data_esrange_winter = data['Esrange'].loc[(data['Esrange'].index.month<6) | (data['Esrange'].index.month>=9)]
    data_esrange_summer = data['Esrange'].loc[(data['Esrange'].index.month>=6) & (data['Esrange'].index.month<9)]
    data_pallas_winter = data['Pallas'].loc[(data['Pallas'].index.month<6) | (data['Pallas'].index.month>=9)]
    data_pallas_summer = data['Pallas'].loc[(data['Pallas'].index.month>=6) & (data['Pallas'].index.month<9)]

    data_janiskoski_winter = data['Janiskoski'].loc[(data['Janiskoski'].index.month<6) | (data['Janiskoski'].index.month>=9)]
    data_janiskoski_summer = data['Janiskoski'].loc[(data['Janiskoski'].index.month>=6) & (data['Janiskoski'].index.month<9)]

    obs_esrange_winter = data_esrange_winter.groupby(data_esrange_winter.index.hour).apply(np.nanmean)
    obserr_esrange_winter = data_esrange_winter.groupby(data_esrange_winter.index.hour).apply(np.nanstd)/np.sqrt(data_esrange_winter.size)

    obs_esrange_summer = data_esrange_summer.groupby(data_esrange_summer.index.hour).apply(np.nanmean)
    obserr_esrange_summer = data_esrange_summer.groupby(data_esrange_summer.index.hour).apply(np.nanstd)/np.sqrt(data_esrange_summer.size)

    obs_pallas_winter = data_pallas_winter.groupby(data_pallas_winter.index.hour).apply(np.nanmean)
    obserr_pallas_winter = data_pallas_winter.groupby(data_pallas_winter.index.hour).apply(np.nanstd)/np.sqrt(data_pallas_winter.size)

    obs_pallas_summer = data_pallas_summer.groupby(data_pallas_summer.index.hour).apply(np.nanmean)
    obserr_pallas_summer = data_pallas_summer.groupby(data_pallas_summer.index.hour).apply(np.nanstd)/np.sqrt(data_pallas_summer.size)

    obs_janiskoski_winter = data_janiskoski_winter.groupby(data_janiskoski_winter.index.hour).apply(np.nanmean)
    obserr_janiskoski_winter = data_janiskoski_winter.groupby(data_janiskoski_winter.index.hour).apply(np.nanstd)/np.sqrt(data_janiskoski_winter.size)

    obs_janiskoski_summer = data_janiskoski_summer.groupby(data_janiskoski_summer.index.hour).apply(np.nanmean)
    obserr_janiskoski_summer = data_janiskoski_summer.groupby(data_janiskoski_summer.index.hour).apply(np.nanstd)/np.sqrt(data_janiskoski_summer.size)
if b_model:
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
obs_karasjok_winter.plot(ax=ax11, ls="None", marker='+', color='orange', alpha=0.5, label='Karasjok')
obs_karasjok_summer.plot(ax=ax12, ls="None", marker='+', color='orange', alpha=0.5, label='Karasjok')

plot_error_bands(ax11,
                 obs_karasjok_winter.index,
                 obs_karasjok_winter.values,
                 obserr_karasjok_winter,
                 color='orange')

plot_error_bands(ax12,
                 obs_karasjok_summer.index,
                 obs_karasjok_summer.values,
                 obserr_karasjok_summer,
                 color='orange')

obs_jergul_winter.plot(ax=ax11, ls="None", marker='x', color='red', alpha=0.5, label='Jergul')
obs_jergul_summer.plot(ax=ax12, ls="None", marker='x', color='red', alpha=0.5, label='Jergul')

plot_error_bands(ax11,
                 obs_jergul_winter.index,
                 obs_jergul_winter.values,
                 obserr_jergul_winter,
                 color='red')

plot_error_bands(ax12,
                 obs_jergul_summer.index,
                 obs_jergul_summer.values,
                 obserr_jergul_summer,
                 color='red')

obs_svanvik_winter.plot(ax=ax11, ls="None", marker='o', color='blueviolet', alpha=0.5, label='Svanvik')
obs_svanvik_summer.plot(ax=ax12, ls="None", marker='o', color='blueviolet', alpha=0.5, label='Svanvik')

plot_error_bands(ax11,
                 obs_svanvik_winter.index,
                 obs_svanvik_winter.values,
                 obserr_svanvik_winter,
                 color='blueviolet')

plot_error_bands(ax12,
                 obs_svanvik_summer.index,
                 obs_svanvik_summer.values,
                 obserr_svanvik_summer,
                 color='blueviolet')

obs_esrange_winter.plot(ax=ax11, ls="None", marker='^', color='blue', alpha=0.5, label='Esrange')
obs_esrange_summer.plot(ax=ax12, ls="None", marker='^', color='blue', alpha=0.5, label='Esrange')

plot_error_bands(ax11,
                 obs_esrange_winter.index,
                 obs_esrange_winter.values,
                 obserr_esrange_winter,
                 color='blue')

plot_error_bands(ax12,
                 obs_esrange_summer.index,
                 obs_esrange_summer.values,
                 obserr_esrange_summer,
                 color='blue')

obs_pallas_winter.plot(ax=ax11, ls="None", marker='v', color='black', alpha=0.5, label='Pallas')
obs_pallas_summer.plot(ax=ax12, ls="None", marker='v', color='black', alpha=0.5, label='Pallas')

plot_error_bands(ax11,
                 obs_pallas_winter.index,
                 obs_pallas_winter.values,
                 obserr_pallas_winter,
                 color='black')

plot_error_bands(ax12,
                 obs_pallas_summer.index,
                 obs_pallas_summer.values,
                 obserr_pallas_summer,
                 color='black')

obs_janiskoski_winter.plot(ax=ax11, ls="None", marker='d', color='grey', alpha=0.5, label='Janiskoski')
obs_janiskoski_summer.plot(ax=ax12, ls="None", marker='d', color='grey', alpha=0.5, label='Janiskoski')

plot_error_bands(ax11,
                 obs_janiskoski_winter.index,
                 obs_janiskoski_winter.values,
                 obserr_janiskoski_winter,
                 color='grey')

plot_error_bands(ax12,
                 obs_janiskoski_summer.index,
                 obs_janiskoski_summer.values,
                 obserr_janiskoski_summer,
                 color='grey')



# Decorate the plots
for ax in fig1.axes:
    ax.set_ylim(10,45)
    ax.set_xlabel("")
    ax.set_ylabel("[$O_3$] (ppb)")
    ax.legend(ncol=6)
    
ax11.set_title("Winter")
ax12.set_title("Summer")
ax12.set_xlabel("Time (hour)")


# Show it
plt.show(block=False)
