# Resample
svanvik_daily = data['Svanvik'].resample('1d').apply(np.nanmean)
svanvik_daily_2018 = data_svanvik_OzoNorClim['2018'].resample('1d').apply(np.nanmean)
svanvik_daily_2019 = data_svanvik_OzoNorClim['2019'].resample('1d').apply(np.nanmean)

svanvik_daily_std = data['Svanvik'].resample('1d').apply(np.nanstd)
svanvik_daily_std_2018 = data_svanvik_OzoNorClim['2018'].resample('1d').apply(np.nanstd)
svanvik_daily_std_2019 = data_svanvik_OzoNorClim['2019'].resample('1d').apply(np.nanstd)

svanvik_daily_stderr = data['Svanvik'].resample('1d').apply(lambda x: np.nanstd(x)/np.sqrt(x.count()))
svanvik_daily_stderr_2018 = data_svanvik_OzoNorClim['2018'].resample('1d').apply(lambda x: np.nanstd(x)/np.sqrt(x.count()))
svanvik_daily_stderr_2019 = data_svanvik_OzoNorClim['2019'].resample('1d').apply(lambda x: np.nanstd(x)/np.sqrt(x.count()))

# Draw sample from climatology of Jergul/Karsjok, Esrange, Pallas -> fig9
sample = fitSpl_dmean(svanvik_daily.dropna().index.dayofyear)
sample_2018 = fitSpl_dmean(svanvik_daily_2018.dropna().index.dayofyear)
sample_2019 = fitSpl_dmean(svanvik_daily_2019.dropna().index.dayofyear)
# Draw sample from Svanvik climatology
sample_2018_svanvik = fitSpl_dmean_svanvik(svanvik_daily_2018.dropna().index.dayofyear)
sample_2019_svanvik = fitSpl_dmean_svanvik(svanvik_daily_2019.dropna().index.dayofyear)


# Select 2018 data -> fig10
esrange_daily_2018 = data['Esrange']['2018'].resample('1d').apply(np.nanmean)
pallas_daily_2018 = data['Pallas']['2018'].resample('1d').apply(np.nanmean)
prestebakke_daily_2018 = data['Prestebakke']['2018'].resample('1d').apply(np.nanmean)

esrange_daily_2018_std = data['Esrange']['2018'].resample('1d').apply(np.nanstd)
pallas_daily_2018_std = data['Pallas']['2018'].resample('1d').apply(np.nanstd)
prestebakke_daily_2018_std = data['Prestebakke']['2018'].resample('1d').apply(np.nanstd)

esrange_daily_2018_stderr = data['Esrange']['2018'].resample('1d').apply(lambda x: np.nanstd(x)/np.sqrt(x.count()))
pallas_daily_2018_stderr = data['Pallas']['2018'].resample('1d').apply(lambda x: np.nanstd(x)/np.sqrt(x.count()))
prestebakke_daily_2018_stderr = data['Prestebakke']['2018'].resample('1d').apply(lambda x: np.nanstd(x)/np.sqrt(x.count()))

# Sample accordingly from climatology
sample_2018_esrange = fitSpl_dmean(esrange_daily_2018.dropna().index.dayofyear)
sample_2018_pallas = fitSpl_dmean(pallas_daily_2018.dropna().index.dayofyear)
sample_2018_prestebakke = fitSpl_dmean_prestebakke(prestebakke_daily_2018.dropna().index.dayofyear)

# Sample only from June-September
esrange_jja = esrange_daily_2018.where((esrange_daily_2018.index.month>=6) & (esrange_daily_2018.index.month<9)).dropna()
pallas_jja = pallas_daily_2018.where((pallas_daily_2018.index.month>=6) & (pallas_daily_2018.index.month<9)).dropna()
prestebakke_jja = prestebakke_daily_2018.where((prestebakke_daily_2018.index.month>=6) & (prestebakke_daily_2018.index.month<9)).dropna()
svanvik_jja = svanvik_daily_2018.where((svanvik_daily_2018.index.month>=6) & (svanvik_daily_2018.index.month<9)).dropna()

sample_jja_esrange = fitSpl_dmean(esrange_jja.index.dayofyear)
sample_jja_pallas = fitSpl_dmean(pallas_jja.index.dayofyear)
sample_jja_prestebakke = fitSpl_dmean_prestebakke(prestebakke_jja.index.dayofyear)
sample_jja_svanvik = fitSpl_dmean_svanvik(svanvik_jja.index.dayofyear)


