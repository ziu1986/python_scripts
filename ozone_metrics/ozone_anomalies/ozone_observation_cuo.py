aot0_svanvik = compute_aot(data['Svanvik'], time_start=1, time_end=23, month_start=6, month_end=7,level=0) + compute_aot(data['Svanvik'], time_start=6, time_end=17, month_start=8, month_end=8,level=0)

aot0_svanvik_OzoNorClim = compute_aot(data_svanvik_OzoNorClim, time_start=1, time_end=23, month_start=6, month_end=7,level=0) + compute_aot(data_svanvik_OzoNorClim, time_start=6, time_end=17, month_start=8, month_end=8,level=0)

aot0_esrange = compute_aot(data['Esrange'], time_start=1, time_end=23, month_start=6, month_end=7,level=0) + compute_aot(data['Esrange'], time_start=6, time_end=17, month_start=8, month_end=8,level=0)

aot0_pallas = compute_aot(data['Pallas'], time_start=1, time_end=23, month_start=6, month_end=7,level=0) + compute_aot(data['Pallas'], time_start=6, time_end=17, month_start=8, month_end=8,level=0)

aot0_jerkara = compute_aot(data_jergkara, time_start=1, time_end=23, month_start=6, month_end=7,level=0) + compute_aot(data_jergkara, time_start=6, time_end=17, month_start=8, month_end=8,level=0)

aot0_prestebakke = compute_aot(data['Prestebakke'], time_start=1, time_end=23, month_start=6, month_end=7,level=0) + compute_aot(data['Prestebakke'], time_start=6, time_end=17, month_start=8, month_end=8,level=0)

aot0_climatology = compute_aot(sample_clim_hourly[0], time_start=1, time_end=23, month_start=6, month_end=7,level=0) + compute_aot(sample_clim_hourly[0], time_start=6, time_end=17, month_start=8, month_end=8,level=0)

aot0_climatology_prestebakke = compute_aot(sample_clim_hourly_prestebakke[0], time_start=1, time_end=23, month_start=6, month_end=7,level=0) + compute_aot(sample_clim_hourly_prestebakke[0], time_start=6, time_end=17, month_start=8, month_end=8,level=0)

aot0_climatology_svanvik = compute_aot(sample_clim_hourly_svanvik[0], time_start=1, time_end=23, month_start=6, month_end=7,level=0) + compute_aot(sample_clim_hourly_svanvik[0], time_start=6, time_end=17, month_start=8, month_end=8,level=0)


# Correction for missing data in 2018 in Svanvik
scaling_max_pallas = ((data['Pallas'].loc['2018-07']-sample_clim_hourly[0].loc['2018-07'])/sample_clim_hourly[0].loc['2018-07'])
scaling_max_esrange =  ((data['Esrange'].loc['2018-07']-sample_clim_hourly[0].loc['2018-07'])/sample_clim_hourly[0].loc['2018-07'])
scaling_max_svanvik = ((data_svanvik_OzoNorClim.loc['2018-07']-sample_clim_hourly_svanvik[0][data_svanvik_OzoNorClim.loc['2018-07'].index])/sample_clim_hourly_svanvik[0][data_svanvik_OzoNorClim.loc['2018-07'].index])

aot0_svanvik_corr = compute_aot((scaling_max_pallas*sample_clim_hourly_svanvik[0]['2018-07']+sample_clim_hourly_svanvik[0]['2018-07'].loc['2018-07-09':'2018-07-23']), time_start=1, time_end=23, month_start=7, month_end=7,level=0)

aot0_svanvik_OzoNorClim = aot0_svanvik_OzoNorClim.add(aot0_svanvik_corr, fill_value=0)

# [aot0] = ppm h or nmol/mol h 
aot0 = pd.DataFrame({"Esrange":aot0_esrange, "Pallas":aot0_pallas, "Jergul/Karasjok":aot0_jerkara, "Prestebakke":aot0_prestebakke, "Svanvik":aot0_svanvik_OzoNorClim}) #, "Svanvik (2018/19)":aot0_svanvik_OzoNorClim

delta_aot0 = pd.DataFrame({"Esrange":aot0_esrange-aot0_climatology.values[0], "Pallas":aot0_pallas-aot0_climatology.values[0], "Jergul/Karasjok":aot0_jerkara-aot0_climatology.values[0], "Prestebakke":aot0_prestebakke-aot0_climatology_prestebakke.values[0], "Svanvik":aot0_svanvik_OzoNorClim-aot0_climatology_svanvik.values[0]})

rel_aot0 = pd.DataFrame({"Esrange":(aot0_esrange-aot0_climatology.values[0])/aot0_climatology.values[0], "Pallas":(aot0_pallas-aot0_climatology.values[0])/aot0_climatology.values[0], "Jergul/Karasjok":(aot0_jerkara-aot0_climatology.values[0])/aot0_climatology.values[0], "Prestebakke":(aot0_prestebakke-aot0_climatology_prestebakke.values[0])/aot0_climatology_prestebakke.values[0], "Svanvik":(aot0_svanvik_OzoNorClim-aot0_climatology_svanvik.values[0])/aot0_climatology_svanvik.values[0]})


# Exclude years where more then 3 percent of data are mising!
for each in ['Esrange', 'Pallas', 'Prestebakke']:
    years = np.arange(2012,2019)#np.unique(data[each].index.year)
    
    for iyear in years:
        opt_obs = len(pd.date_range("%d-06" % iyear,"%d-08-31 23:0" % iyear, freq='H'))
        missing_data = (len(data[each].loc["%d-06" % iyear:"%d-08" % iyear].dropna())-opt_obs)/float(opt_obs)*100
        if missing_data < -3:
            #print(each, iyear, missing_data)
            #print(delta_aot0[each][iyear])
            aot0[each][iyear] = np.nan
            delta_aot0[each][iyear] = np.nan


