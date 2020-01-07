# Calculate climatology until 2017 based on Esrange, Pallas and Jergul/Karasjok data -> fig8
ozone_days = []
doys = []
ozone_days_svanvik = []
doys_svanvik = []
# Bin hourly data into daily for each doy -> max number of points per bin num_years x 24hours
for idoy in np.arange(1,367):
    # Get the data for all relevant stations and concatenate them
    ozone_days.append(data['Esrange'][:'2017'].where(data['Esrange'][:'2017'].index.dayofyear==idoy).dropna().values)
    doys.append(np.full(ozone_days[-1].size,idoy))
    ozone_days.append(data['Pallas'][:'2017'].where(data['Pallas'][:'2017'].index.dayofyear==idoy).dropna().values)
    doys.append(np.full(ozone_days[-1].size,idoy))
    ozone_days.append(data_jergkara.where(data_jergkara.index.dayofyear==idoy).dropna().values)
    doys.append(np.full(ozone_days[-1].size,idoy))
    # Get the date for Svanvik
    ozone_days_svanvik.append(data['Svanvik'].where(data['Svanvik'].index.dayofyear==idoy).dropna().values)
    doys_svanvik.append(np.full(ozone_days_svanvik[-1].size,idoy))

doys = np.concatenate(np.array(doys))
ozone_days = np.concatenate(np.array(ozone_days))
doys_svanvik = np.concatenate(np.array(doys_svanvik))
ozone_days_svanvik = np.concatenate(np.array(ozone_days_svanvik))


