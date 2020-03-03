# Fit the distributions
x_sample, pdf, fit, stat = fit_skew_normal((svanvik_daily.dropna()-sample).values)
x_sample_2018, pdf_2018, fit_2018, stat_2018 = fit_skew_normal((svanvik_daily_2018.dropna()-sample_2018).values)
x_sample_2019, pdf_2019, fit_2019, stat_2019 = fit_skew_normal((svanvik_daily_2019.dropna()-sample_2019).values)
x_sample_svanvik, pdf_svanvik, fit_svanvik, stat_svanvik = fit_skew_normal((svanvik_daily_2018.dropna()-sample_2018_svanvik).values)
x_sample_svanvik_2019, pdf_svanvik_2019, fit_svanvik_2019, stat_svanvik_2019 = fit_skew_normal((svanvik_daily_2019.dropna()-sample_2019_svanvik).values)

# Fit the distributions for fennoscandia and Prestebakke
x_sample_esrange, pdf_esrange, fit_esrange, stat_esrange = fit_skew_normal((esrange_daily_2018.dropna()-sample_2018_esrange).values)
x_sample_pallas, pdf_pallas, fit_pallas, stat_pallas = fit_skew_normal((pallas_daily_2018.dropna()-sample_2018_pallas).values)
x_sample_prestebakke, pdf_prestebakke, fit_prestebakke, stat_prestebakke = fit_skew_normal((prestebakke_daily_2018.dropna()-sample_2018_prestebakke).values)


# Fit the test distribution
fitable_range = lambda x : x[np.where((x >= -20) & (x <= 20))[0]]
x_test, pdf_test, fit_test, stat_test = fit_skew_normal(fitable_range(score_svanvik))
x_test_2018, pdf_test_2018, fit_test_2018, stat_test_2018 = fit_skew_normal(fitable_range(score_svanvik_2018))
x_test_2019, pdf_test_2019, fit_test_2019, stat_test_2019 = fit_skew_normal(fitable_range(score_svanvik_2019))
x_test_svanvik, pdf_test_svanvik, fit_test_svanvik, stat_test_svanvik = fit_skew_normal(fitable_range(score_svanvik_clim_2018))
x_test_svanvik_2019, pdf_test_svanvik_2019, fit_test_svanvik_2019, stat_test_svanvik_2019 = fit_skew_normal(fitable_range(score_svanvik_clim_2019))

# Fit the test distributions for fennoscandia and Prestebakke
x_test_esrange_2018, pdf_test_esrange_2018, fit_test_esrange_2018, stat_test_esrange_2018 = fit_skew_normal(fitable_range(score_esrange_2018))
x_test_pallas_2018, pdf_test_pallas_2018, fit_test_pallas_2018, stat_test_pallas_2018 = fit_skew_normal(fitable_range(score_pallas_2018))
x_test_prestebakke_2018, pdf_test_prestebakke_2018, fit_test_prestebakke_2018, stat_test_prestebakke_2018 = fit_skew_normal(fitable_range(score_prestebakke_2018))

# Fit the hourly residuals
x_test_pallas_hourly_2018, pdf_test_pallas_hourly_2018, fit_test_pallas_hourly_2018, stat_test_pallas_hourly_2018 = fit_skew_normal(fitable_range((data['Pallas'].loc['2018-07'].values-clim_hourly.loc[7].values.flatten())/clim_hourly_err.loc[7].values.flatten()))
x_test_esrange_hourly_2018, pdf_test_esrange_hourly_2018, fit_test_esrange_hourly_2018, stat_test_esrange_hourly_2018 = fit_skew_normal(fitable_range((data['Esrange'].loc['2018-07'].values-clim_hourly.loc[7].values.flatten())/clim_hourly_err.loc[7].values.flatten()))
x_test_svanvik_hourly_2018, pdf_test_svanvik_hourly_2018, fit_test_svanvik_hourly_2018, stat_test_svanvik_hourly_2018 = fit_skew_normal(fitable_range((data_svanvik_OzoNorClim.loc['2018-07']-sample_clim_hourly_svanvik[0][data_svanvik_OzoNorClim.loc['2018-07'].index])/sample_clim_hourly_err_svanvik[0][data_svanvik_OzoNorClim.loc['2018-07'].index]))

data_svanvik_rra = (data_rra.sel(lat=station_location['Svanvik'].lat, lon=station_location['Svanvik'].lon, method='nearest', time='2018-07')['O3']*0.5)

x_test_svanvik_rra_hourly_2018, pdf_test_svanvik_rra_hourly_2018, fit_test_svanvik_rra_hourly_2018, stat_test_svanvik_rra_hourly_2018 = fit_skew_normal(fitable_range((data_svanvik_rra-sample_clim_hourly_svanvik[0][data_svanvik_rra.time.values])/sample_clim_hourly_err_svanvik[0][data_svanvik_rra.time.values]))
