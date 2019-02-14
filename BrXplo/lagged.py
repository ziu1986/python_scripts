import os, glob # Access environment variables
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates # handling of dates in plotting
import datetime as dt  # Python standard library datetime module
import scipy.stats
from mytools.met_tools import *

def read_file(infile, **kargs):
    verbose = kargs.pop('v', False)
    file = open(infile, 'r')
    lines = file.readlines()
    raw_data = []
    for line in lines:
        cols = line.strip()
        raw_data.append(float(cols))
    return raw_data

def time_lagg_corr(data_1, data_truth, **kargs):
    lag = kargs.pop('lag',0)
    verbose = kargs.pop('v', False)
    if lag >= 0:
        corr_coef = np.ma.corrcoef(np.roll(data_1, lag)[lag:], data_truth[lag:])
    else:
        corr_coef = np.ma.corrcoef(np.roll(data_1, lag)[:lag], data_truth[:lag])
    corr_sign = corr_coef[0,1]*np.sqrt((len(data_1)-np.fabs(lag)-2)/(1-corr_coef[0,1]**2))
    if verbose:
        print "%d %d %1.2f -> %2.2f" % (lag, len(data_1), corr_coef[0,1], corr_sign)
    return corr_coef


src = os.environ['DATA']
src_dir = '/BrXplo/'
file_emac_ref = ['emac_ref_alert.txt','emac_ref_barrow.txt']
file_emac = ['emac_alert.txt','emac_barrow.txt']
file_emac_rs = ['emac_rs_alert.txt','emac_rs_barrow.txt']
file_obs = ['obs_alert.txt', 'obs_barrow.txt']

data_obs = np.ma.masked_values([read_file(src+src_dir+infile) for infile in file_obs],(-9999999.9,))
data_obs = np.ma.masked_where(data_obs>100, data_obs)
data_emac = np.array([read_file(src+src_dir+infile) for infile in file_emac])
data_emac_ref = np.array([read_file(src+src_dir+infile) for infile in file_emac_ref])
data_emac_rs = np.array([read_file(src+src_dir+infile) for infile in file_emac_rs])

lag_span = np.arange(-120,121)
data_lag_alert = [time_lagg_corr(data_emac[0],data_obs[0],lag=lag)[0,1] for lag in lag_span]
corr_max_alert = (lag_span[np.argmax(data_lag_alert)], np.max(data_lag_alert))
data_lag_barrow = [time_lagg_corr(data_emac[1],data_obs[1],lag=lag)[0,1] for lag in lag_span]
corr_max_barrow = (lag_span[np.argmax(data_lag_barrow)], np.max(data_lag_barrow))

data_lag_alert_ref = [time_lagg_corr(data_emac_ref[0],data_obs[0],lag=lag)[0,1] for lag in lag_span]
corr_max_alert_ref = (lag_span[np.argmax(data_lag_alert_ref)], np.max(data_lag_alert_ref))
data_lag_barrow_ref = [time_lagg_corr(data_emac_ref[1],data_obs[1],lag=lag)[0,1] for lag in lag_span]
corr_max_barrow_ref = (lag_span[np.argmax(data_lag_barrow_ref)], np.max(data_lag_barrow_ref))

data_lag_alert_rs = [time_lagg_corr(data_emac_rs[0],data_obs[0],lag=lag)[0,1] for lag in lag_span]
corr_max_alert_rs = (lag_span[np.argmax(data_lag_alert_rs)], np.max(data_lag_alert_rs))
data_lag_barrow_rs = [time_lagg_corr(data_emac_rs[1],data_obs[1],lag=lag)[0,1] for lag in lag_span]
corr_max_barrow_rs = (lag_span[np.argmax(data_lag_barrow_rs)], np.max(data_lag_barrow_rs))
# Plot it
plt.close('all')

fig1 = plt.figure(1,figsize=(12,8))
fig1.canvas.set_window_title("lagged_correlation_all")
ax11 = plt.subplot()
ax11.plot(lag_span, data_lag_alert_ref, color='black', label="Alert_ref", ls='-.')
ax11.vlines(corr_max_alert_ref[0], 0, corr_max_alert_ref[1], color='black', linestyle=':')
ax11.text(corr_max_alert_ref[0]-4, corr_max_alert_ref[1]+0.001, '%d h' % corr_max_alert_ref[0], color='black', size='large')
ax11.plot(lag_span, data_lag_alert, color='blue', label="Alert_mysic", ls='-.')
ax11.vlines(corr_max_alert[0], 0, corr_max_alert[1], color='blue', linestyle=':')
ax11.text(corr_max_alert[0]-4, corr_max_alert[1]+0.001, '%d h' % corr_max_alert[0], color='blue', size='large')
ax11.plot(lag_span, data_lag_alert_rs, label="Alert_rs", color='c', ls='-.')
ax11.vlines(corr_max_alert_rs[0], 0, corr_max_alert_rs[1], color='c', linestyle=':')
ax11.text(corr_max_alert_rs[0]-4, corr_max_alert_rs[1]+0.001, '%d h' % corr_max_alert_rs[0], color='c', size='large')

ax11.plot(lag_span, data_lag_barrow_ref, color='black', label='Utqia$\mathrm{\dot{g}}$vik_ref')
ax11.vlines(corr_max_barrow_ref[0], 0, corr_max_barrow_ref[1], color='black', linestyle=':')
ax11.text(corr_max_barrow_ref[0]-2, corr_max_barrow_ref[1]+0.001, '%d h' % corr_max_barrow_ref[0], color='black', size='large')
ax11.plot(lag_span, data_lag_barrow, color='blue', label='Utqia$\mathrm{\dot{g}}$vik_mysic')
ax11.vlines(corr_max_barrow[0], 0, corr_max_barrow[1], color='blue', linestyle=':')
ax11.text(corr_max_barrow[0]-2, corr_max_barrow[1]+0.001, '%d h' % corr_max_barrow[0], color='blue', size='large')
ax11.plot(lag_span, data_lag_barrow_rs, color='c', label='Utqia$\mathrm{\dot{g}}$vik_rs')
ax11.vlines(corr_max_barrow_rs[0], 0, corr_max_barrow_rs[1], color='c', linestyle=':')
ax11.text(corr_max_barrow_rs[0]-2, corr_max_barrow_rs[1]+0.001, '%d h' % corr_max_barrow_rs[0], color='c', size='large')

ax11.set_xlabel("Lag (hours)")
ax11.set_xticks(lag_span[::24])
ax11.set_xlim(-120,120)
ax11.set_ylabel("Correlation coefficient")
ax11.legend(loc=(0,0), ncol=2)
ax11.axhline(0, ls='--', color='black')

plt.show(block=False)





