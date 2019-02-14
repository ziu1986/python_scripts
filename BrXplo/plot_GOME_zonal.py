# Clean up
plt.close('all')
BrO_data_zonalmean = BrO_data_mean.mean(dim='lon')
fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("GOME_BrO_tot_zonal_mean")
ax21 = plt.subplot(221)
ax21.set_title("JFM")
for imonth in (1,2,3):
    ax21.plot(BrO_data_zonalmean.lat, BrO_data_zonalmean.sel(month=imonth).BrO, label=month_name[imonth-1])
ax22 = plt.subplot(222)
ax22.set_title("AMJ")
for imonth in (4,5,6):
    ax22.plot(BrO_data_zonalmean.lat, BrO_data_zonalmean.sel(month=imonth).BrO, label=month_name[imonth-1])
ax23 = plt.subplot(223)
ax23.set_title("JAS")
for imonth in (7,8,9):
    ax23.plot(BrO_data_zonalmean.lat, BrO_data_zonalmean.sel(month=imonth).BrO, label=month_name[imonth-1])
ax24 = plt.subplot(224)
ax24.set_title("OND")
for imonth in (10,11,12):
    ax24.plot(BrO_data_zonalmean.lat, BrO_data_zonalmean.sel(month=imonth).BrO, label=month_name[imonth-1])

for ax in fig2.axes:
    ax.set_xticks(np.arange(-90, 91, 15))
    ax.set_xlim(-90,90)
    ax.set_ylim(0,7)
    ax.axvline(-45, color='grey', ls='--')
    ax.axvline(45, color='grey', ls='--')
    ax.legend(frameon=False,loc='lower right')
ax24.set_xlabel("Latitude (deg)", x=-0.1)
ax23.set_ylabel("%s (%s %s)" % ('Vertical column BrO', '1e+13', 'molecules/cm2'), y=1.1)
plt.show(block=False)
