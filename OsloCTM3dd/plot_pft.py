import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
import datetime as dt
from mytools.met_tools import *
from mytools.netcdf_tools import *
from collections import OrderedDict

pft_data_dir = "./pft_landusedyn.3x3.nc"

try:
    pft_data
except NameError:
    # Read plant function types
    pft_data = read_data(pft_data_dir)
# Get PFT selection for use as masks
try:
    pft_sel
except NameError:
    pft_sel = []
    for each in pft_data['PFT_PCT'].sel(time=2005.00):
        print("Changing coordinates...")
        each.coords['x'] = each.lat.lon[0].data
        each.coords['y'] = each.lat[:,0].data
        each = each.drop(('lon','lat'))
        each = each.rename({'x':'lon','y':'lat'})
        pft_sel.append(each)

#"0: Barren land";
#"1: Needleaf evergreen temperate tree";
#"2: Needleaf evergreen boreal tree";
#"3: Needleaf deciduous boreal tree";
#"4: Broadleaf evergreen tropical tree";
#"5: Broadleaf evergreen temperate tree";
#"6: Broadleaf deciduous tropical tree";
#"7: Broadleaf deciduous temperate tree";
#"8: Broadleaf deciduous boreal tree";
#"9: Broadleaf evergreen temperate shrub";
#"10: Broadleaf deciduous temperate shrub";
#"11: Broadleaf deciduous boreal shrub";
#"12: Arctic C3 grass (cold)";
#"13: C3 grass (cool)";
#"14: C4 grass (warm)";
#"15: Crop1";
#"16: Crop2";

pft_de = pft_sel[0].where((pft_sel[0].lat >=-60) & (pft_sel[0].lat <= 70))
pft_si = pft_sel[0].where(~((pft_sel[0].lat >=-60) & (pft_sel[0].lat <= 70)))
pft_oc = (xr.concat(pft_sel,dim='pft').sum(dim='pft')-100)*-1
pft_cf = xr.concat(pft_sel[1:4], dim='pft').sum(dim='pft') #xr.concat(pft_sel[10], dim='pft')
pft_df = xr.concat(pft_sel[7:11], dim='pft').sum(dim='pft') #xr.concat(pft_sel[11:13], dim='pft')
pft_tf = xr.concat(pft_sel[4:7], dim='pft').sum(dim='pft')
pft_ac = pft_sel[15]
pft_gr = xr.concat(pft_sel[13:15], dim='pft').sum(dim='pft')
pft_tu = xr.concat(pft_sel[11:13], dim='pft').sum(dim='pft')

pft_de_zonal = pft_de.mean(dim='lon')
pft_si_zonal = pft_si.mean(dim='lon')
pft_oc_zonal = pft_oc.mean(dim='lon')
pft_cf_zonal = pft_cf.mean(dim='lon')
pft_df_zonal = pft_df.mean(dim='lon')
pft_tf_zonal = pft_tf.mean(dim='lon')
pft_ac_zonal = pft_ac.mean(dim='lon')
pft_gr_zonal = pft_gr.mean(dim='lon')
pft_tu_zonal = pft_tu.mean(dim='lon')

pft_zonal = (pft_cf_zonal, pft_df_zonal, pft_tf_zonal, pft_ac_zonal, pft_gr_zonal, pft_tu_zonal, pft_si_zonal, pft_oc_zonal, pft_de_zonal)

labels = ("Coniferous forest", "Deciduous forest", "Tropical forest", "Cropland", "Grassland", "Tundra", "Ice/Snow", "Ocean", "Desert")
colors = ("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999")
#colors = ("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6")


# Clean up
plt.close('all')

# Plot it
# Pfts
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("pft_categories")

ax11 = plt.subplot(331)
ax12 = plt.subplot(332)
ax13 = plt.subplot(333)
ax14 = plt.subplot(334)
ax15 = plt.subplot(335)
ax16 = plt.subplot(336)
ax17 = plt.subplot(337)
ax18 = plt.subplot(338)
ax19 = plt.subplot(339)

(pft_cf/100).plot(ax=ax11, levels=np.arange(0,1.01,0.1), cmap=plt.cm.nipy_spectral_r, extend='max')
(pft_df/100).plot(ax=ax12, levels=np.arange(0,1.01,0.1), cmap=plt.cm.nipy_spectral_r, extend='max')
(pft_tf/100).plot(ax=ax13, levels=np.arange(0,1.01,0.1), cmap=plt.cm.nipy_spectral_r, extend='max')
(pft_ac/100).plot(ax=ax14, levels=np.arange(0,1.01,0.1), cmap=plt.cm.nipy_spectral_r, extend='max')
(pft_gr/100).plot(ax=ax15, levels=np.arange(0,1.01,0.1), cmap=plt.cm.nipy_spectral_r, extend='max')
(pft_tu/100).plot(ax=ax16, levels=np.arange(0,1.01,0.1), cmap=plt.cm.nipy_spectral_r, extend='max')
(pft_si/100).plot(ax=ax17, levels=np.arange(0,1.01,0.1), cmap=plt.cm.nipy_spectral_r, extend='max')
(pft_oc/100).plot(ax=ax18, levels=np.arange(0,1.01,0.1), cmap=plt.cm.nipy_spectral_r, extend='max')
(pft_de/100).plot(ax=ax19, levels=np.arange(0,1.01,0.1), cmap=plt.cm.nipy_spectral_r, extend='max')

ax11.set_title("Coniferous forest", y=0.97)
ax12.set_title("Deciduous forest", y=0.97)
ax13.set_title("Tropical forest", y=0.97)
ax14.set_title("Cropland", y=0.97)
ax15.set_title("Grassland", y=0.97)
ax16.set_title("Tundra", y=0.97)
ax17.set_title("Ice/Snow", y=0.97)
ax18.set_title("Ocean", y=0.97)
ax19.set_title("Desert", y=0.97)

for ax in fig1.axes:
    ax.set_xlabel("")
    ax.set_xticks(np.arange(0,361,60))
    ax.set_ylabel("")
    ax.set_yticks(np.arange(-90,91,45))
ax14.set_ylabel("Latitude (deg)")
ax18.set_xlabel("Longitude (deg)")

# Zonal plot
import matplotlib.gridspec as gridspec
fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("pft_categories_zonal")

gs0 = gridspec.GridSpec(1, 5, fig2)
ax21 = plt.subplot2grid((1,5),(0,0),colspan=4)
#ax21 = plt.subplot(121)
bottom = np.zeros_like(pft_de_zonal.data)
pbar = []
for i in np.arange(len(pft_zonal)):
    #pft_zonal[i].plot(label=labels[i])
    p1 = plt.bar(pft_zonal[i].lat, pft_zonal[i].data, bottom=bottom, width=2, color=colors[i])
    pbar.append(p1[0])
    bottom +=  np.nan_to_num(pft_zonal[i].data)
    
#ax21 = plt.gca()
ax21.set_ylim(-5,105)
ax21.set_xlabel("Latitude (deg)")
ax21.set_ylabel("PFT (%)")
plt.legend(pbar,labels,bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Show it
plt.show(block=False)

