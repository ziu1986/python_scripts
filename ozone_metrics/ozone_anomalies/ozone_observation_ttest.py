score_svanvik = ((svanvik_daily.dropna()-sample)/svanvik_daily_std).dropna().values
score_svanvik_2018 = ((svanvik_daily_2018.dropna()-sample_2018)/svanvik_daily_std_2018).dropna().values
score_svanvik_2019 = ((svanvik_daily_2019.dropna()-sample_2019)/svanvik_daily_std_2019).dropna().values

score_svanvik_clim_2018 = ((svanvik_daily_2018.dropna()-sample_2018_svanvik)/svanvik_daily_std_2018).dropna().values
score_svanvik_clim_2019 = ((svanvik_daily_2019.dropna()-sample_2019_svanvik)/svanvik_daily_std_2019).dropna().values



score_esrange_2018 = ((esrange_daily_2018.dropna()-sample_2018_esrange)/esrange_daily_2018_std).dropna().values
score_pallas_2018 = ((pallas_daily_2018.dropna()-sample_2018_pallas)/pallas_daily_2018_std).dropna().values
score_prestebakke_2018 = ((prestebakke_daily_2018.dropna()-sample_2018_prestebakke)/prestebakke_daily_2018_std).dropna().values
