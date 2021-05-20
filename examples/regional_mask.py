import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import regionmask

lon = np.arange(-179.5, 180)
lat = np.arange(-89.5, 90)

land = regionmask.defined_regions.natural_earth.land_110

land.plot()

mask = land.mask(np.arange(14,33), np.arange(65,71.5))

f, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
ax.coastlines()

mask.plot(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False);

plt.show(block=False)
