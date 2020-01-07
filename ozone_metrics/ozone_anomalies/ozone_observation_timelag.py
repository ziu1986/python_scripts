# Time lags -> fig6
time_lag = range(-32,33)
lag_jergkara_esrange = []
lag_jergkara_pallas = []
lag_svanvik_esrange = []
lag_svanvik_pallas = []
lag_svanvik_jergkara = []
lag_svanvik_janiskoski = []
lag_jergkara_janiskoski = []

lag_label = ("jergkara_esrange","jergkara_pallas","svanvik_esrange","svanvik_pallas","svanvik_jergkara","lag_svanvik_janiskoski","lag_jergkara_janiskoski")
for i in time_lag:
    #print("%d %1.2f" % (i, time_lagged_corr(data_jergkara, data['Esrange'], lag=i, pandas=True)))
    lag_jergkara_esrange.append(time_lagged_corr(data_jergkara, data['Esrange'], lag=i, pandas=True))
    lag_jergkara_pallas.append(time_lagged_corr(data_jergkara, data['Pallas'], lag=i, pandas=True))
    lag_svanvik_esrange.append(time_lagged_corr(data['Svanvik'], data['Esrange'], lag=i, pandas=True))
    lag_svanvik_pallas.append(time_lagged_corr(data['Svanvik'], data['Pallas'], lag=i, pandas=True))
    lag_svanvik_jergkara.append(time_lagged_corr(data['Svanvik'], data_jergkara, lag=i, pandas=True))
    lag_svanvik_janiskoski.append(time_lagged_corr(data['Svanvik'], data['Janiskoski'], lag=i, pandas=True))
    lag_jergkara_janiskoski.append(time_lagged_corr(data_jergkara, data['Janiskoski'], lag=i, pandas=True))
# Print maximum in lag
print("Lag correlation")
for i,lag in zip(lag_label,(lag_jergkara_esrange, lag_jergkara_pallas, lag_svanvik_esrange, lag_svanvik_pallas, lag_svanvik_jergkara, lag_svanvik_janiskoski, lag_jergkara_janiskoski)):
    print("%s max at %d h" % (i, np.array(time_lag)[np.where(np.array(lag)==np.array(lag).max())[0]]))


