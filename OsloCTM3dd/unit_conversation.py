import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from mytools.met_tools import *

plt.close('all')

t_range = np.arange(-30,51)+273.15
p_range = np.arange(913,1114)

t_mesh, p_mesh = np.meshgrid(t_range, p_range)
converter = R*t_mesh/p_mesh*1/hecto

# Plot it
fig1 = plt.figure(1, figsize=(14,9))
fig1.canvas.set_window_title("Molar_density_conversion")
ax11 = plt.subplot(121)
im11 = ax11.imshow(1/converter)
cbar11 = plt.colorbar(im11)

ax12 = plt.subplot(122)
im12 = ax12.imshow((1/converter-41)/41*100, vmin=-100., vmax=100., cmap=plt.cm.seismic)
cbar12 = plt.colorbar(im12)

for ax in fig1.axes[::2]:
    ax.set_xlabel("Temperature (K)")
    ax.set_xticks(np.arange(0,len(t_range),10))
    ax.set_xticklabels((t_range[::10]).astype(int))
    ax.set_ylabel("Pressure (hPa)")
    ax.set_yticks(np.arange(0,len(p_range),25))
    ax.set_yticklabels((p_range[::25]).astype(int)/hecto)
    ax.axhline(100,color='black', ls='--')
    ax.axvline(50,color='black', ls='--')

fig1.axes[1].set_ylabel("Molar Densisty (mol m$^{-3}$)")
fig1.axes[-1].set_ylabel("$\Delta$Molar Densisty/41 (%)")
plt.show(block=False)
