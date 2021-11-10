
import os, sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

src = os.environ['MODELS'] + "/icon/build_lndonly/experiments/jsbalone_R2B4_sfa/jsbalone_R2B4_sfa_lnd_basic_ml_20000201T000000Z.nc"
data = xr.open_dataset(src)

x = data['clon'].data
y = data['clat'].data
lai = data['pheno_lai_veg'].isel(time=0).data

f, ax = plt.subplots(1,2, sharex=True, sharey=True)
ax[0].tripcolor(x,y,lai)
ax[1].tricontourf(x,y,lai, 20) # choose 20 contour levels, just to show how good its interpolation is

plt.show(block=False)