# Plot it
## Clean up
plt.close('all')

if plot_timeseries:
    fig1 = plt.figure(1, figsize=(16,9))
    fig1.canvas.set_window_title("ozone_timeseries_ltobs")
    ax11 = plt.subplot(311)
    ax12 = plt.subplot(312, sharex=ax11)
    ax13 = plt.subplot(313, sharex=ax11)
    data["Barrow"].plot(ax=ax11, ls='None', marker='x', label='Utqiagvik (USA)')
    data['Prestebakke'].plot(ax=ax12, zorder=2, ls='None', marker='.', label='Prestebakke (NOR)', color='red')
    data_jergkara.plot(ax=ax13, zorder=2, ls='None', marker='+', label='Jergul/Karasjok (NOR)', color='orange')

    ax12.axvspan(date2num(data['Prestebakke'].index[0].date()), date2num(dt.datetime.strptime('1996-12','%Y-%m')), zorder=3, facecolor='None', edgecolor='black', hatch='//')
    ax13.axvspan(date2num(data_jergkara.index[0].date()), date2num(dt.datetime.strptime('1996-12','%Y-%m')), zorder=3, facecolor='None', edgecolor='black', hatch='//')
    ax13.set_xlabel("Time (year)")
    ax12.set_ylabel("[$O_3$] (ppb)")

    for ax in fig1.axes:
        ax.set_ylim(0,100)
        ax.legend()

    #    
    fig2 = plt.figure(2, figsize=(16,9))
    fig2.canvas.set_window_title("ozone_timeseries_fenoscandic_obs")
    ax21 = plt.subplot(511)
    ax22 = plt.subplot(512, sharex=ax21)
    ax23 = plt.subplot(513, sharex=ax21)
    ax24 = plt.subplot(514, sharex=ax21)
    ax25 = plt.subplot(515, sharex=ax21)
    ax21.plot(data['Esrange'].index, data['Esrange'], ls='None', marker='+', label='Esrange (SWE)', color='blue')
    ax22.plot(data['Pallas'].index, data['Pallas'], ls='None', marker='+', label='Pallas (FIN)', color='black')
    data_jergkara.plot(ax=ax23, ls='None', marker='+', label='Jergul/Karasjok (NOR)', color='orange')
    ax24.plot(data['Svanvik'].index, data['Svanvik'], ls='None', marker='+', label='Svanvik (NOR)', color='blueviolet')
    ax24.plot(data_svanvik_OzoNorClim.index, data_svanvik_OzoNorClim, ls='None', marker='x', label='Svanvik, 2018/19 (NOR)', color='blueviolet')
    ax25.plot(data['Janiskoski'].index, data['Janiskoski'], ls='None', marker='+', label='Janiskoski (RUS)', color='grey')

    ax23.axvspan(date2num(data_jergkara.index[0].date()), date2num(dt.datetime.strptime('1996-12','%Y-%m')), zorder=3, facecolor='None', edgecolor='black', hatch='//')
    ax24.axvspan(date2num(data['Svanvik'].index[0].date()), date2num(dt.datetime.strptime('1996-12','%Y-%m')), zorder=3, facecolor='None', edgecolor='black', hatch='//')

    ax25.set_xlabel("Time (year)")
    ax23.set_ylabel("[$O_3$] (ppb)", y=1)
    for ax in fig2.axes:
        ax.set_ylim(0,100)
        ax.legend(ncol=2)

#
if plot_climatology:
    fig3 = plt.figure(3, figsize=(16,9))
    fig3.canvas.set_window_title("ozone_climatology_fenoscandic_obs")
    ax31 = plt.subplot()

    data['Prestebakke'].groupby(data['Prestebakke'].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Prestebakke (NOR)', color='red')
    data['Esrange'].groupby(data['Esrange'].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Esrange (SWE)', color='blue', ls='none', marker="v")
    data['Pallas'].groupby(data['Pallas'].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Pallas (FIN)', color='black', ls='none', marker="^")
    data_jergkara.groupby(data_jergkara.index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Jergul/Karasjok (NOR)', color='orange', ls='none', marker="x")
    data['Svanvik'].groupby(data['Svanvik'].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Svanvik (NOR)', color='blueviolet', ls='none', marker="+")
    data['Janiskoski'].groupby(data['Janiskoski'].index.dayofyear).apply(np.nanmean).plot(ax=ax31, label='Janiskoski (RUS)', color='grey', ls='none', marker="d")

    ax31.set_xlabel("Time (day of year)")
    ax31.set_ylabel("[$O_3$] (ppb)")
    for ax in fig3.axes:
        ax.set_ylim(0,60)
        ax.legend(loc='lower left')

    plot_month_span(ax31)
    plot_month_name(ax31, ypos=58)

if plot_spectrum:
    fig20 = plt.figure(20, figsize=(16,9))
    fig20.canvas.set_window_title("fequency_spectrum")
    ax201 = plt.subplot()
    ax201.stem(1/freqs_barrow/12, np.abs(fft_barrow)/np.abs(fft_barrow).max(), label='Utqiagvik (USA)')
    markerline, stemlines, baseline = ax201.stem(1/freqs_prestebakke/12, np.abs(fft_prestebakke)/np.abs(fft_prestebakke).max(), label='Prestebakke (NOR)')
    plt.setp(stemlines, color='red')
    plt.setp(markerline, color='red')
    markerline, stemlines, baseline = ax201.stem(1/freqs_jergkara/12, np.abs(fft_jergkara)/np.abs(fft_jergkara).max(), label='Jergul/Karasjok (NOR)')
    plt.setp(stemlines, color='orange')
    plt.setp(markerline, color='orange')
    ax201.set_xlim(0,40)
    ax201.set_ylabel("Normalized Amplitude")
    ax201.set_xlabel("Frequency (years)")

    ax201.legend()

#
if plot_map:
    # Show the stations
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cartopy.io.img_tiles as cimgt
    #import cartopy.geodesic as cgeo
    from matplotlib.transforms import offset_copy

    fig4 = plt.figure(4)
    fig4.canvas.set_window_title("station_map_fennoscandia")
    stamen_terrain = cimgt.Stamen('terrain-background')
    ax41 = fig4.add_subplot(1,1,1, projection=stamen_terrain.crs) #ccrs.PlateCarree()
    ax41.set_extent([19.4,31.4,67.6,71.4], crs=ccrs.PlateCarree())
    ax41.add_image(stamen_terrain, 8)
    
    ax41.plot(station_location['Jergul'].lon, station_location['Jergul'].lat, fillstyle='left', marker='o', color='orange', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax41.plot(station_location['Karasjok'].lon, station_location['Karasjok'].lat, fillstyle='right', marker='o', color='orange', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax41.plot(station_location['Svanvik'].lon, station_location['Svanvik'].lat, marker='o', color='blueviolet', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax41.plot(station_location['Esrange'].lon, station_location['Esrange'].lat, marker='o', color='blue', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax41.plot(station_location['Pallas'].lon, station_location['Pallas'].lat, marker='o', color='black', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    ax41.plot(station_location['Janiskoski'].lon, station_location['Janiskoski'].lat, marker='o', color='grey', markersize=12, alpha=0.7, transform=ccrs.Geodetic())
    geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax41)
    text_transform = offset_copy(geodetic_transform, units='dots', x=-25)
    text_transform_2 = offset_copy(geodetic_transform, units='dots', y=25)

    # Add text.
    ax41.text(station_location['Jergul'].lon, station_location['Jergul'].lat, u'Jergul',
              verticalalignment='center', horizontalalignment='right',
              transform=text_transform,
              bbox=dict(facecolor='orange', alpha=0.5, boxstyle='round'))
    ax41.text(station_location['Karasjok'].lon, station_location['Karasjok'].lat, u'Karasjok',
              verticalalignment='center', horizontalalignment='right',
              transform=text_transform_2,
              bbox=dict(facecolor='orange', alpha=0.5, boxstyle='round'))
    ax41.text(station_location['Svanvik'].lon, station_location['Svanvik'].lat, u'Svanvik',
              verticalalignment='center', horizontalalignment='right',
              transform=text_transform,
              bbox=dict(facecolor='blueviolet', alpha=0.5, boxstyle='round'))
    ax41.text(station_location['Esrange'].lon, station_location['Esrange'].lat, u'Esrange',
              verticalalignment='center', horizontalalignment='right',
              transform=text_transform,
              bbox=dict(facecolor='blue', alpha=0.5, boxstyle='round'))
    ax41.text(station_location['Pallas'].lon, station_location['Pallas'].lat, u'Pallas',
              verticalalignment='center', horizontalalignment='right',
              transform=text_transform,
              bbox=dict(facecolor='black', alpha=0.5, boxstyle='round'))
    ax41.text(station_location['Janiskoski'].lon, station_location['Janiskoski'].lat, u'Janiskoski',
              verticalalignment='center', horizontalalignment='right',
              transform=text_transform,
              bbox=dict(facecolor='grey', alpha=0.5, boxstyle='round'))


#
if plot_correlation:
    fig5 = plt.figure(5, figsize=(16,9))
    fig5.canvas.set_window_title("density_distribution")
    ax51 = plt.subplot(221)
    ax52 = plt.subplot(222)
    ax53 = plt.subplot(223)
    ax54 = plt.subplot(224)
    bins = range(81)
    hist2d_1 = ax51.hist2d(data['Esrange'].dropna(), data_jergkara[data['Esrange'].dropna().index], bins=(bins), cmap=plt.cm.hot_r)
    hist2d_2 = ax52.hist2d(data['Pallas'].dropna(), data_jergkara[data['Pallas'].dropna().index], bins=(bins), cmap=plt.cm.hot_r)
    hist2d_3 = ax53.hist2d(data['Esrange'].dropna(), data['Svanvik'][data['Esrange'].dropna().index], bins=(bins), cmap=plt.cm.hot_r)
    hist2d_4 = ax54.hist2d(data['Pallas'].dropna(), data['Svanvik'][data['Pallas'].dropna().index], bins=(bins), cmap=plt.cm.hot_r)

    for ax in fig5.axes:
        ax.plot(bins, bins, color='grey',ls='--')
    
    # Placing the colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    for ax,hist in zip(fig5.axes, (hist2d_1,hist2d_2,hist2d_3,hist2d_4)):
        axins = inset_axes(ax,
                           width="2.5%",  # width = 50% of parent_bbox width
                           height="70%",  # height : 5%
                           loc='lower right',
                           #bbox_to_anchor=(1, 0., 1, 1),
                           #bbox_transform=ax.transAxes,
                       #borderpad=0,
        )
        fig5.colorbar(hist[-1], cax=axins, orientation='vertical')
        axins.set_ylabel("Counts")

    ax51.set_ylabel("$[O_3]_{Jergul/Karasjok}$ (ppb)")
    ax53.set_ylabel("$[O_3]_{Svanvik}$ (ppb)")
    ax53.set_xlabel("$[O_3]_{Esrange}$ (ppb)")
    ax54.set_xlabel("$[O_3]_{Pallas}$ (ppb)")
    
    ax51.text(60,75, "$r^2 = %1.2f$" % (data['Esrange'].corr(data_jergkara)), size='large')
    ax52.text(60,75, "$r^2 = %1.2f$" % (data['Pallas'].corr(data_jergkara)), size='large')
    ax53.text(60,75, "$r^2 = %1.2f$" % (data['Esrange'].corr(data['Svanvik'])), size='large')
    ax54.text(60,75, "$r^2 = %1.2f$" % (data['Pallas'].corr(data['Svanvik'])), size='large')

#
if plot_timelag:
    fig6 = plt.figure(figsize=(16,9))
    fig6.canvas.set_window_title("ozone_observation_timelag")
  
    ax61 = plt.subplot(121)
    ax62 = plt.subplot(122)
    #ax63 = plt.subplot(133)
    ax61.set_title("Jergul/Karasjok")
    ax62.set_title("Svanvik")
    #ax63.set_title("Svanvik")
    ax61.plot(time_lag, lag_jergkara_esrange, color='blue', label='Esrange')
    ax61.plot(time_lag, lag_jergkara_pallas, color='black', label='Pallas')
    ax61.plot(np.array(time_lag)*(-1), lag_svanvik_jergkara, ls='--', color='blueviolet', label='Svanvik')
    ax61.plot(time_lag, lag_jergkara_janiskoski, color='grey', ls='--', label='Janiskoski')

    ax62.plot(time_lag, lag_svanvik_esrange, color='blue', ls='--', label='Esrange')
    ax62.plot(time_lag, lag_svanvik_pallas, color='black', ls='--', label='Pallas')
    ax62.plot(time_lag, lag_svanvik_jergkara, color='orange', ls='--', label='Jergul/Karasjok')
    ax62.plot(time_lag, lag_svanvik_janiskoski, color='grey', label='Janiskoski')

    for ax in fig6.axes:
        ax.set_xlabel('Lag (hours)')
        #ax.set_ylim(0,1)
        ax.legend(ncol=2)
    ax61.set_ylabel('Correlation Coefficient')
#
if plot_climatology:
    fig7 = plt.figure(7, figsize=(16,9))
    fig7.canvas.set_window_title("ozone_climatology_fenoscandic_obs_norm")
    ax71 = plt.subplot()

    (data['Prestebakke'].groupby(data['Prestebakke'].index.dayofyear).apply(np.nanmean)-data['Prestebakke'].mean()).plot(ax=ax71, label='Prestebakke (NOR)', color='red')
    (data['Esrange'].groupby(data['Esrange'].index.dayofyear).apply(np.nanmean)-data['Esrange'].mean()).plot(ax=ax71, label='Esrange (SWE)', color='blue', ls='none', marker="v")
    (data['Pallas'].groupby(data['Pallas'].index.dayofyear).apply(np.nanmean)-data['Pallas'].mean()).plot(ax=ax71, label='Pallas (FIN)', color='black', ls='none', marker="^")
    (data_jergkara.groupby(data_jergkara.index.dayofyear).apply(np.nanmean)-data_jergkara.mean()).plot(ax=ax71, label='Jergul/Karasjok (NOR)', color='orange', ls='none', marker="x")
    (data['Svanvik'].groupby(data['Svanvik'].index.dayofyear).apply(np.nanmean)-data['Svanvik'].mean()).plot(ax=ax71, label='Svanvik (NOR)', color='blueviolet', ls='none', marker="+")
    (data['Janiskoski'].groupby(data['Janiskoski'].index.dayofyear).apply(np.nanmean)-data['Janiskoski'].mean()).plot(ax=ax71, label='Janiskoski (RUS)', color='grey', ls='none', marker="d")

    ax71.set_xlabel("Time (day of year)")
    ax71.set_ylabel("$\Delta$[$O_3$] (ppb)")
    for ax in fig7.axes:
        #ax.set_ylim(0,60)
        ax.legend(loc='lower left')

    plot_month_span(ax71)
    plot_month_name(ax71, ypos=20)

#
if plot_splines:
    fig8 = plt.figure(8, figsize=(10,12))
    fig8.canvas.set_window_title("ozone_climatology_fenoscandic_obs_spline")
    ax81 = plt.subplot(211)
    ax82 = plt.subplot(212)

    xtime = np.arange(1,367)

    hist81 = ax81.hist2d(doys, ozone_days, bins=(np.arange(0.5,367),np.linspace(0,81,160)), cmap=plt.cm.hot_r)
    ax81.errorbar(xtime[::10], yozone[::10], yerr=yerr[::10], color='black', marker='o', fillstyle='none', ls='None', label='monthly mean (1 $\sigma$)')
    ax81.errorbar(xtime[::10], yozone[::10], yerr=yerr_mean[::10], color='cornflowerblue', marker='None', ls='None', label='std err mean')
    ax81.errorbar(xtime[::10]-1, yozone_max[::10], yerr=yerr_max[::10], color='darkred', marker='v', ls='None', label='monthly mean max (1 $\sigma$)')
    ax81.errorbar(xtime[::10]+1, yozone_min[::10], yerr=yerr_min[::10], color='darkblue', marker='^', ls='None', label='monthly mean min (1 $\sigma$)')
    ax81.plot(np.linspace(1,367), fitSpl_dmean(np.linspace(1,367)), ls='--', label='spline fit: daily mean')
    ax81.plot(np.linspace(1,367), fitSpl_dmax(np.linspace(1,367)), ls=':', label='spline fit: mean daily max', color='red')
   

    hist82 = ax82.hist2d(doys_svanvik, ozone_days_svanvik, bins=(np.arange(0.5,367),np.linspace(0,81,160)), cmap=plt.cm.hot_r)
    ax82.plot(np.linspace(1,367), fitSpl_dmean_svanvik(np.linspace(1,367)), ls='--', label='spline fit: daily mean')
    ax82.plot(np.linspace(1,367), fitSpl_dmax_svanvik(np.linspace(1,367)), ls=':', label='spline fit: mean daily max', color='red')
    ax82.errorbar(xtime[::10], yozone_svanvik[::10], yerr=yerr_svanvik[::10], color='black', marker='o', ls='None', label='monthly mean (1 $\sigma$)')
    ax82.errorbar(xtime[::10], yozone_svanvik[::10], yerr=yerr_mean_svanvik[::10], color='cornflowerblue', marker='None', ls='None', label='std err mean')
    ax82.errorbar(xtime[::10]-1, yozone_max_svanvik[::10], yerr=yerr_max_svanvik[::10], color='darkred', marker='v', ls='None', label='monthly mean max (1 $\sigma$)')
    ax82.errorbar(xtime[::10]+1, yozone_min_svanvik[::10], yerr=yerr_min_svanvik[::10], color='darkblue', marker='^', ls='None', label='monthly mean min (1 $\sigma$)')
    

    # Setting colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    for ax,hist in zip(fig8.axes, (hist81, hist82)):
        axins = inset_axes(ax,
                           width="30%",  # width = 50% of parent_bbox width
                           height="2.5%",  # height : 5%
                           loc='upper left',
                           #bbox_to_anchor=(0.01, 0.95, 1, 1),
                           #bbox_transform=ax.transAxes,
                           #borderpad=0,
        )
        fig8.colorbar(hist[-1], cax=axins, orientation='horizontal')
        axins.set_xlabel("Counts")
    plot_month_name(ax81, ypos=82)

    ax82.set_xlabel("Time (day of year)")
    ax81.set_ylabel("[$O_3$] (ppb)", y=-0.25)

    ax81.legend()
    ax82.legend()
    ax82.set_title("Svanvik")

#
if plot_residuals:
    fig9 = plt.figure(9, figsize=(10,12))
    fig9.canvas.set_window_title("ozone_climatology_fenoscandic_obs_residuals-Svanvik")
    ax91 = plt.subplot(211)
    ax92 = plt.subplot(212)

    hist91 = ax91.hist((svanvik_daily.dropna()-sample).values, density=True, color='blueviolet', label='Svanvik clim.')
    hist91_2018 = ax91.hist((svanvik_daily_2018.dropna()-sample_2018).values, density=True, histtype='step', color='violet', label='Svanvik 2018')
    hist91_2019 = ax91.hist((svanvik_daily_2019.dropna()-sample_2019).values, density=True, histtype='step', color='purple', label='Svanvik 2019')

    hist92 = ax92.hist((svanvik_daily_2018.dropna()-sample_2018_svanvik).values, density=True, color='violet', label='Svanvik 2018')
    hist92 = ax92.hist((svanvik_daily_2019.dropna()-sample_2019_svanvik).values, density=True, histtype='step', color='purple', label='Svanvik 2019')

    ax91.set_xlabel("$<O_3>_{daily}^{Svanvik}-O_3^{clim}$ (ppb)")
    ax92.set_xlabel("$<O_3>_{daily}^{Svanvik}-{O_3}^{clim}_{Svanvik}$ (ppb)")
    ax91.set_ylabel("Probability density", y=-0.25)

    ax91.plot(x_sample, pdf, color='black', label='Skew normal fit')
    stats_text(ax91, stat, fit, name="Svanvik clim.", ypos=0.7)

    ax91.plot(x_sample_2018, pdf_2018, color='black', ls='--', label='Skew normal fit: 2018')
    stats_text(ax91, stat_2018, fit_2018, name="Svanvik 2018", ypos=0.4)

    ax91.plot(x_sample_2019, pdf_2019, color='black', ls='-.', label='Skew normal fit: 2019')
    stats_text(ax91, stat_2019, fit_2019, name="Svanvik 2019", ypos=0.1)

    ax92.plot(x_sample_svanvik, pdf_svanvik, color='black', label='Skew normal fit: 2018')
    stats_text(ax92, stat_svanvik, fit_svanvik, ypos=0.7)

    ax92.plot(x_sample_svanvik_2019, pdf_svanvik_2019, color='black', label='Skew normal fit: 2019')
    stats_text(ax92, stat_svanvik_2019, fit_svanvik_2019, ypos=0.4)

    ax91.legend()
    ax92.legend()

    fig10 = plt.figure(10, figsize=(10,12))
    fig10.canvas.set_window_title("ozone_climatology_fenoscandic_obs_residuals")
    ax101 = plt.subplot(311)
    ax102 = plt.subplot(312, sharex=ax101)
    ax103 = plt.subplot(313, sharex=ax101)

    ax101.set_title("Esrange", y=0.85, x=0.05)
    ax102.set_title("Pallas", y=0.85, x=0.05)
    ax103.set_title("Prestebakke", y=0.85, x=0.05)

    hist101_2018 = ax101.hist((esrange_daily_2018.dropna()-sample_2018_esrange).values, density=True, histtype='step', color='blue', label='Esrange 2018')
    hist102_2018 = ax102.hist((pallas_daily_2018.dropna()-sample_2018_pallas).values, density=True, histtype='step', color='black', label='Pallas 2018')
    hist103_2018 = ax103.hist((prestebakke_daily_2018.dropna()-sample_2018_prestebakke).values, density=True, histtype='step', color='red', label='Prestebakke 2018')

    ax101.plot(x_sample_esrange, pdf_esrange, color='black', label='Skew normal fit')
    stats_text(ax101, stat_esrange, fit_esrange, ypos=0.4)

    ax102.plot(x_sample_pallas, pdf_pallas, color='black', label='Skew normal fit')
    stats_text(ax102, stat_pallas, fit_pallas, ypos=0.4)

    ax103.plot(x_sample_prestebakke, pdf_prestebakke, color='black', label='Skew normal fit')
    stats_text(ax103, stat_prestebakke, fit_prestebakke, ypos=0.4)

    ax103.set_xlabel("$<O_3>_{daily}-O_3^{clim}$ (ppb)")

if plot_aot:
    fig11 = plt.figure(11, figsize=(16,9))
    fig11.canvas.set_window_title("ozone_fenoscandic_obs_aot40")
    ax111 = plt.subplot()
    compute_aot(data['Prestebakke'], month_start=6, month_end=8).plot(label='Prestebakke (NOR)', color='red')
    compute_aot(data_jergkara, month_start=6, month_end=8).plot(label='Jergul/Karasjok (NOR)', color='orange')
    compute_aot(data['Pallas'], month_start=6, month_end=8).plot(label='Pallas (FIN)', color='black')
    compute_aot(data['Esrange'], month_start=6, month_end=8).plot(label='Esrange (SWE)', color='blue')
    compute_aot(data['Svanvik'], month_start=6, month_end=8).plot(label='Svanvik (NOR)', color='blueviolet')
    compute_aot(data_svanvik_OzoNorClim, month_start=6, month_end=8).plot(label='', color='blueviolet')

    compute_aot(data_jergkara, time_start=1, time_end=23, month_start=6, month_end=8).plot(ls='--', color='orange', label='')
    compute_aot(data['Pallas'], time_start=1, time_end=23, month_start=6, month_end=8).plot(ls='--', color='black', label='')
    compute_aot(data['Esrange'], time_start=1, time_end=23, month_start=6, month_end=8).plot(ls='--', color='blue', label='')
    compute_aot(data['Svanvik'], time_start=1, time_end=23, month_start=6, month_end=8).plot(ls='--', color='blueviolet', label='')
    compute_aot(data_svanvik_OzoNorClim, time_start=1, time_end=23, month_start=6, month_end=8).plot(ls='--', color='blueviolet', label='')


    ax111.set_xlabel("Time (years)")
    ax111.set_ylabel("AOT40 (ppb h)")
    ax111.legend()
    ax111.axhline(3000, ls=':', color='black')

if plot_rollingsum:
    fig12 = plt.figure(12, figsize=(16,14))
    fig12.canvas.set_window_title("ozone_fenoscandic_obs_rolling_sum40")
    ax121 = plt.subplot(211)
    ax122 = plt.subplot(212)

    try:
        AOT40_820
    except NameError:
        AOT40_820 = {}
        AOT40_820['Prestebakke'] = compute_aot(data['Prestebakke'], month_start=5, month_end=9, rolling=True)
        AOT40_820['Jergkara'] = compute_aot(data_jergkara, month_start=5, month_end=9, rolling=True)
        AOT40_820['Pallas'] = compute_aot(data['Pallas'], month_start=5, month_end=9, rolling=True)
        AOT40_820['Esrange'] = compute_aot(data['Esrange'], month_start=5, month_end=9, rolling=True)
        AOT40_820['Svanvik'] = compute_aot(data['Svanvik'], month_start=5, month_end=9, rolling=True)
        #compute_aot(data_svanvik_OzoNorClim, month_start=5, month_end=9, rolling=True)
    
        AOT40_123 = {}
        AOT40_123['Prestebakke'] = compute_aot(data['Prestebakke'], time_start=1, time_end=23, month_start=5, month_end=9, rolling=True)
        AOT40_123['Jergkara'] = compute_aot(data_jergkara, time_start=1, time_end=23, month_start=5, month_end=9, rolling=True)
        AOT40_123['Pallas'] = compute_aot(data['Pallas'], time_start=1, time_end=23, month_start=5, month_end=9, rolling=True)
        AOT40_123['Esrange'] = compute_aot(data['Esrange'], time_start=1, time_end=23, month_start=5, month_end=9, rolling=True)
        AOT40_123['Svanvik'] = compute_aot(data['Svanvik'], time_start=1, time_end=23, month_start=5, month_end=9, rolling=True)
        #compute_aot(data_svanvik_OzoNorClim, time_start=1, time_end=23, month_start=5, month_end=9, rolling=True)
    
        AOT30_820 = {}
        AOT30_820['Prestebakke'] = compute_aot(data['Prestebakke'], month_start=5, month_end=9, level=30, rolling=True)
        AOT30_820['Jergkara'] = compute_aot(data_jergkara, month_start=5, month_end=9, level=30, rolling=True)
        AOT30_820['Pallas'] = compute_aot(data['Pallas'], month_start=5, month_end=9, level=30, rolling=True)
        AOT30_820['Esrange'] = compute_aot(data['Esrange'], month_start=5, month_end=9, level=30, rolling=True)
        AOT30_820['Svanvik'] = compute_aot(data['Svanvik'], month_start=5, month_end=9, level=30, rolling=True)
        #compute_aot(data_svanvik_OzoNorClim, month_start=5, month_end=9, rolling=True)

        AOT30_123 = {}
        AOT30_123['Prestebakke'] = compute_aot(data['Prestebakke'], time_start=1, time_end=23, month_start=5, month_end=9, level=30, rolling=True)
        AOT30_123['Jergkara'] = compute_aot(data_jergkara, time_start=1, time_end=23, month_start=5, month_end=9, level=30, rolling=True)
        AOT30_123['Pallas'] = compute_aot(data['Pallas'], time_start=1, time_end=23, month_start=5, month_end=9, level=30, rolling=True)
        AOT30_123['Esrange'] = compute_aot(data['Esrange'], time_start=1, time_end=23, month_start=5, month_end=9, level=30, rolling=True)
        AOT30_123['Svanvik'] = compute_aot(data['Svanvik'], time_start=1, time_end=23, month_start=5, month_end=9, level=30, rolling=True)
    #compute_aot(data_svanvik_OzoNorClim, time_start=1, time_end=23, month_start=5, month_end=9, rolling=True)
    #(AOT30_820['Prestebakke']['2005']).plot(ax=ax121, label='Prestebakke (NOR)', ls='-', color='red')
    (AOT30_820['Jergkara']['2005']).plot(ax=ax121,label='Jergul/Karasjok (NOR)', ls='-', color='orange')
    (AOT30_820['Pallas']['2005']).plot(ax=ax121,label='Pallas (FIN)', ls='-', color='black')
    (AOT30_820['Esrange']['2005']).plot(ax=ax121,label='Esrange (SWE)', ls='-', color='blue')
    
    #(AOT30_123['Prestebakke']['2005']).plot(ax=ax121, ls='-.', color='red', label='')
    (AOT30_123['Jergkara']['2005']).plot(ax=ax121,ls='-.', color='orange', label='')
    (AOT30_123['Pallas']['2005']).plot(ax=ax121,ls='-.', color='black', label='')
    (AOT30_123['Esrange']['2005']).plot(ax=ax121,ls='-.', color='blue', label='')
    
    #(AOT40_820['Prestebakke']['2005']).plot(ax=ax121, label='', ls='--', color='red')
    (AOT40_820['Jergkara']['2005']).plot(ax=ax121,label='', ls='--', color='orange')
    (AOT40_820['Pallas']['2005']).plot(ax=ax121,label='', ls='--', color='black')
    (AOT40_820['Esrange']['2005']).plot(ax=ax121,label='', ls='--', color='blue')

    #(AOT40_123['Prestebakke']['2005']).plot(ax=ax121, ls=':', color='red', label='')
    (AOT40_123['Jergkara']['2005']).plot(ax=ax121,ls=':', color='orange', label='')
    (AOT40_123['Pallas']['2005']).plot(ax=ax121,ls=':', color='black', label='')
    (AOT40_123['Esrange']['2005']).plot(ax=ax121,ls=':', color='blue', label='')

    
    #(AOT30_820['Prestebakke']['2006']).plot(ax=ax122, label='Prestebakke (NOR)', ls='-', color='red')
    (AOT30_820['Jergkara']['2006']).plot(ax=ax122,label='Jergul/Karasjok (NOR)', ls='-', color='orange')
    (AOT30_820['Pallas']['2006']).plot(ax=ax122,label='Pallas (FIN)', ls='-', color='black')
    (AOT30_820['Esrange']['2006']).plot(ax=ax122,label='Esrange (SWE)', ls='-', color='blue')

    #(AOT30_123['Prestebakke']['2006']).plot(ax=ax122, ls='-.', color='red', label='')
    (AOT30_123['Jergkara']['2006']).plot(ax=ax122,ls='-.', color='orange', label='')
    (AOT30_123['Pallas']['2006']).plot(ax=ax122,ls='-.', color='black', label='')
    (AOT30_123['Esrange']['2006']).plot(ax=ax122,ls='-.', color='blue', label='')
    
    #(AOT40_820['Prestebakke']['2006']).plot(ax=ax122, label='', ls='--', color='red')
    (AOT40_820['Jergkara']['2006']).plot(ax=ax122,label='', ls='--', color='orange')
    (AOT40_820['Pallas']['2006']).plot(ax=ax122,label='', ls='--', color='black')
    (AOT40_820['Esrange']['2006']).plot(ax=ax122,label='', ls='--', color='blue')

    #(AOT40_123['Prestebakke']['2006']).plot(ax=ax122, ls=':', color='red', label='')
    (AOT40_123['Jergkara']['2006']).plot(ax=ax122,ls=':', color='orange', label='')
    (AOT40_123['Pallas']['2006']).plot(ax=ax122,ls=':', color='black', label='')
    (AOT40_123['Esrange']['2006']).plot(ax=ax122,ls=':', color='blue', label='')

    #ax121.plot(dt.date(2005, 5, 1), 21000, ls='-', color='grey', label='AOT30_8-20')

    ax122.set_xlabel("Time (months)")


    # Add legend
    #lines = ax121.get_lines()
    #legend1 = plt.legend([lines[i] for i in [0,1,2]], ["algo1", "algo2", "algo3"], loc=1)
    #legend2 = plt.legend([lines[i] for i in [0,3,6]], ["balgo1", "balgo2", "balgo3"], loc=4)
    #ax121.add_artist(legend1)
    #ax121.add_artist(legend2)

    for ax, year in zip(fig12.axes, (2005, 2006)):
        ax.set_title("%d" % year)
        ax.set_xlim(dt.date(year, 5, 1), dt.date(year, 7, 10))
        ax.set_ylim(0,20000)
        text_y = ax.get_ylim()[1]
        ax.axvline(dt.date(day=11, month=5, year=year), color='black', ls=':')
        ax.axvspan(dt.date(day=10, month=5, year=year), dt.date(day=12, month=5, year=year), color='black', alpha=0.25)
        ax.text(dt.date(day=10, month=5, year=year), text_y-1000, "Extrapol. SGS 2100")
        ax.axhline(3000, ls=':', color='black')
        ax.axvline(dt.date(day=1, month=7, year=year), color='orange', ls=':')
        ax.text(dt.date(day=1, month=7, year=year), text_y-1000, "SGS present", color='orange')
        ax.set_ylabel("SUM$x_i$ (ppb h)")

    ax121.legend(ncol=3)
    plt.legend(('-','--','-.',':'),('AOT30_8-20','AOT40_8-20','AOT30_1-23','AOT40_1-23'),bbox_to_anchor=(1.05, 1), loc='lower right', borderaxespad=0.)


    fig13 =  plt.figure(13, figsize=(16,14))
    fig13.canvas.set_window_title("sum40_int3month_change_greeningseason_2003-2012_esrange")
    start_year_esrange = 2003#(np.array(AOT30_820['Esrange'].keys())).astype(int).min()
    for i in range(10):
        ax = plt.subplot(3,4,i+1)
        year = start_year_esrange+i
        ax.set_title('%d' % (year))
        (AOT30_820['Esrange']['%d' % (year)]).plot(ax=ax, ls='-', color='blue', label='AOT30_8-20')
        (AOT40_820['Esrange']['%d' % (year)]).plot(ax=ax, ls='--', color='blue', label='AOT40_8-20')
        (AOT30_123['Esrange']['%d' % (year)]).plot(ax=ax, ls='-.', color='blue', label='AOT30_1-23')
        (AOT40_123['Esrange']['%d' % (year)]).plot(ax=ax, ls=':', color='blue', label='AOT40_1-23')

        ax.set_xlim(dt.date(year, 5, 1), dt.date(year, 7, 10))
        ax.set_ylim(0,20000)
        text_y = ax.get_ylim()[1]
        ax.axvline(dt.date(day=11, month=5, year=year), color='black', ls=':')
        ax.axvspan(dt.date(day=10, month=5, year=year), dt.date(day=12, month=5, year=year), color='black', alpha=0.25)
        ax.text(dt.date(day=10, month=5, year=year), text_y-1000, "Extrapol. SGS 2100")
        ax.axhline(3000, ls=':', color='black')
        ax.axvline(dt.date(day=1, month=7, year=year), color='orange', ls=':')
        ax.text(dt.date(day=1, month=7, year=year), text_y-1000, "SGS present", color='orange')
        ax.set_ylabel("SUM$x_i$ (ppb h)", y=-1)
        ax.legend(bbox_to_anchor=(3.15, -2.05), loc='lower right', borderaxespad=0.)

        if i==0:
            print("%s \t %s \t %s \t %s" % ("AOT30_820", "AOT40_820", "AOT30_123", "AOT40_123" ))
        
        print("%3.2f \t %3.2f \t %3.2f \t %3.2f" % (
            ((AOT30_820['Esrange']['%d' % (year)]['%d%s' % (year,"-05-10")])/(AOT30_820['Esrange']['%d' % (year)]['%d%s' % (year,"-07-01")])-1)*100,
            ((AOT40_820['Esrange']['%d' % (year)]['%d%s' % (year,"-05-10")])/(AOT40_820['Esrange']['%d' % (year)]['%d%s' % (year,"-07-01")])-1)*100,
            ((AOT30_123['Esrange']['%d' % (year)]['%d%s' % (year,"-05-10")])/(AOT30_123['Esrange']['%d' % (year)]['%d%s' % (year,"-07-01")])-1)*100,
            ((AOT40_123['Esrange']['%d' % (year)]['%d%s' % (year,"-05-10")])/(AOT40_123['Esrange']['%d' % (year)]['%d%s' % (year,"-07-01")])-1)*100
    ))
        #date_fmt = '%m'
        #formatter = dates.DateFormatter(date_fmt)
        #ax.xaxis.set_major_formatter(formatter)

        plt.gcf().autofmt_xdate()
    

    for ax in fig13.axes[1:]:
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.get_legend().remove()

    for ax in [ fig13.axes[i] for i in (1, 2, 3, 5, 6, 7, 9) ]:
        ax.set_yticklabels('')

if plot_ttest:
    fig14 = plt.figure(14, figsize=(10,12))
    fig14.canvas.set_window_title("ozone_climatology_fenoscandic_obs_residuals-Svanvik")
    ax141 = plt.subplot(211)
    ax142 = plt.subplot(212)

    hist141 = ax141.hist(((svanvik_daily.dropna()-sample)/svanvik_daily_stderr).dropna().values, range=(-10,10), density=True, color='blueviolet', label='Svanvik clim.')
    hist141_2018 = ax141.hist(((svanvik_daily_2018.dropna()-sample_2018)/svanvik_daily_stderr_2018).dropna().values, range=(-10,10), density=True, histtype='step', color='violet', label='Svanvik 2018')
    hist141_2019 = ax141.hist(((svanvik_daily_2019.dropna()-sample_2019)/svanvik_daily_stderr_2019).dropna().values, range=(-10,10), density=True, histtype='step', color='purple', label='Svanvik 2019')

    hist142 = ax142.hist(((svanvik_daily_2018.dropna()-sample_2018_svanvik)/svanvik_daily_stderr_2018).dropna().values, range=(-10,10), density=True, color='violet', label='Svanvik 2018')
    hist142 = ax142.hist(((svanvik_daily_2019.dropna()-sample_2019_svanvik)/svanvik_daily_stderr_2019).dropna().values, range=(-10,10), density=True, histtype='step', color='purple', label='Svanvik 2019')

    ax141.set_xlabel("$t_{score}$")
    ax142.set_xlabel("$t_{score}^{clim_{Svanvik}}$")
    ax141.set_ylabel("Probability density", y=-0.25)

    ax141.legend()
    ax142.legend()

# Show it
plt.show(block=False)





