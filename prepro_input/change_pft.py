import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from mytools.met_tools import *
from mytools.netcdf_tools import *

pft_data_dir = os.environ['DATA']+"/astra_data/processed_data/pftlandusedyn.0.5x0.5.gcp2015.simyr1860-2015.c160715.nc"

try:
    pft_data
except NameError:
    # Read plant function types
    pft_data = read_data(pft_data_dir)


pft_new = pft_data.copy()
pft_zonal = pft_data.mean(dim='lon')
pft_zonal_2005 = pft_zonal.sel(time=2005)

colors = ("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99", "#d3d303", "#fbceb6", "#b15928", "#6cf7c5", "#00b373", "#999999", "#000000")
labels = ("Barren land",
          "Needle. evergr. temp.",
          "Needle. evergr. boreal", 
          "Needle. decid. boreal",
          "Broad. evergr. trop.",
          "Broad. evergr. temp.",
          "Broad. decid. trop.",
          "Broad. decid. temp.",
          "Broad. decid. boreal",
          "Broad. evergr. temp. sh.",
          "Broad. decid. temp. sh.",
          "Broad. decid. boreal sh.",
          "Arctic C3 grass (cold)",
          "C3 grass (cool)",
          "C4 grass (warm)",
          "Crop1",
          "Crop2")
# Clean up
plt.close('all')
# Plot zonals
import matplotlib.gridspec as gridspec
fig1 = plt.figure(1, figsize=(16,9))
gs0 = gridspec.GridSpec(1, 5, fig1)
ax11 = plt.subplot2grid((1,5),(0,0),colspan=4)

stack = np.zeros_like(pft_zonal_2005.isel(pft=0)['PFT_PCT'].values)
pbar = []

for i in np.arange(pft_zonal_2005['pft'].size):
    p1 = plt.bar(np.linspace(-90,90,pft_zonal_2005.isel(pft=i).lat.size), pft_zonal_2005.isel(pft=i)['PFT_PCT'].values, bottom=stack, width=2, color=colors[i])
    pbar.append(p1[0])
    stack +=  np.nan_to_num(pft_zonal_2005.isel(pft=i)['PFT_PCT'].values)

ax11.set_ylim(-5,105)
ax11.set_xlim(-90,90)
ax11.set_xlabel("Latitude (deg)")
ax11.set_ylabel("PFT (%)")
plt.legend(pbar,labels,bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
# Show it
plt.show(block=False)
