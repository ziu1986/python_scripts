# Climatology from Esrange, Pallas, Jergul/Karasjok data
climatology = pd.concat((data['Esrange'][:'2012'], data['Pallas'][:'2012'], data_jergkara[:'2012']))


# Daily mean climatology from Esrange, Pallas, Jergul/Karasjok data
yozone, yerr, yerr_mean = compute_climatology(climatology)
yozone_max, yerr_max, yerr_mean_max = compute_climatology(climatology, mode='max')
yozone_min, yerr_min, yerr_mean_min = compute_climatology(climatology, mode='min')

# Svanvik climatology
yozone_svanvik, yerr_svanvik, yerr_mean_svanvik = compute_climatology(data['Svanvik'])
yozone_max_svanvik, yerr_max_svanvik, yerr_mean_max_svanvik = compute_climatology(data['Svanvik'], mode='max')
yozone_min_svanvik, yerr_min_svanvik, yerr_mean_min_svanvik = compute_climatology(data['Svanvik'], mode='min')

# Prestebakke
yozone_prestebakke, yerr_prestebakke, yerr_mean_prestebakke = compute_climatology(data['Prestebakke'])

# Hourly climatology
clim_hourly, clim_hourly_err, clim_hourly_err_mean = compute_climatology(climatology, mode='hourly')
clim_hourly_svanvik, clim_hourly_err_svanvik, clim_hourly_err_mean_svanvik = compute_climatology(data['Svanvik'], mode='hourly')
clim_hourly_prestebakke, clim_hourly_err_prestebakke, clim_hourly_err_mean_prestebakke = compute_climatology(data['Prestebakke'][:'2012'], mode='hourly')

# Compute spline fits
from scipy.interpolate import UnivariateSpline
# Fennoscandic climatology
w = 1/yerr_mean
fitSpl_dmean = UnivariateSpline(np.unique(doys), climatology.groupby(climatology.index.dayofyear).apply(np.nanmean), w=w)
dmax = climatology.resample('1d').apply(np.nanmax)
fitSpl_dmax = UnivariateSpline(np.unique(doys), dmax.groupby(dmax.index.dayofyear).apply(np.nanmean))
# SVanvik
w_svanvik = 1/yerr_mean_svanvik
fitSpl_dmean_svanvik = UnivariateSpline(np.unique(doys_svanvik), data['Svanvik'].groupby(data['Svanvik'].index.dayofyear).apply(np.nanmean), w=w_svanvik)
dmax_svanvik = data['Svanvik'].resample('1d').apply(np.nanmax)
fitSpl_dmax_svanvik = UnivariateSpline(np.unique(doys_svanvik), dmax_svanvik.groupby(dmax_svanvik.index.dayofyear).apply(np.nanmean))
# Prestebakke
w_prestebakke = 1/yerr_mean_prestebakke
doys_prestebakke = data['Prestebakke'][:'2012'].index.dayofyear
fitSpl_dmean_prestebakke = UnivariateSpline(np.unique(doys_prestebakke), data['Prestebakke'][:'2012'].groupby(doys_prestebakke).apply(np.nanmean), w=w_prestebakke)
dmax_prestebakke = data['Prestebakke'].resample('1d').apply(np.nanmax)

# Compute Savgol filtered data
#from scipy.signal import savgol_filter
#ytest = savgol_filter(climatology.groupby(climatology.index.dayofyear).apply(np.nanmean),31,3)

# Pickle splines for comparison with other data
import pickle
with open('obs_climatologies.pkl','wb') as output:
    pickle.dump(fitSpl_dmean, output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(fitSpl_dmean_svanvik, output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(fitSpl_dmean_prestebakke, output, pickle.HIGHEST_PROTOCOL)
    
    pickle.dump(yerr_mean, output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(yerr_mean_svanvik, output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(yerr_mean_prestebakke, output, pickle.HIGHEST_PROTOCOL)
