import os, sys, glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all, get_month_name, seconds_in_month

# Clean up
plt.close('all')

# Source
src = os.environ['DATA'] + "/astra_data/ECMWF/ERA5/*.nc"

for file in(glob.glob(src)):
    data = xr.open_dataset(file)

selection = data.sel(latitude=69.5, longitude=30.0, method="nearest")
# Plot it
#fig1 = plt.figure(1)
#ax11 = plt.subplot()

#for iyear, icolor in zip(('2018','2019'), ('violet','purple')):
#    selection.sel(time=iyear)['tcc'].plot.hist(ax=ax11, color=icolor, label=iyear)


#plt.legend()

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("fraction_total_cloud_cover")
ax21 = plt.subplot()



pd_sel = []
for iyear in ('2018','2019'):
    hours_month = np.array([seconds_in_month(i, int(iyear)) for i in range(1,13)])/60**2
    pd_sel.append(selection.sel(time=iyear)['tcc'].groupby(selection.sel(time=iyear).time.dt.month).sum().to_dataframe()['tcc'].values/hours_month)

pd_sel = pd.DataFrame({'2018':pd_sel[0], '2019':pd_sel[1]})

pd_sel.plot.bar(ax=ax21, color=('violet','purple'))

plt.xticks(range(12), [get_month_name(i, length=3) for i in range(1,13)], rotation='horizontal')
ax21.set_xlabel("Time (months)")
ax21.set_ylabel("Faction of max. monthly cloud cover")

# Show it
plt.show(block=False)



