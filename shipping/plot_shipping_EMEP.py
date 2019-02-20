import os, glob, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt   # Plotting data
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.line_styles import *

# Clean-up
plt.close('all')

# Data source
src = "EMEP_shipping_emission_data.txt"
# Open and read file
data = pd.read_csv(src,';', header=1)
# Select data
ocean_accros = data['# Format: ISO2'].unique()
ocean_names = ('Baltic Sea', 'Black Sea', 'Mediterranean Sea', 'North Sea', 'Atlantic Ocean')
oceans = dict(zip(ocean_accros,ocean_names))

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("Total_shipping_emissions_EMEP")
ax11 = plt.subplot()
plt.subplots_adjust(right=0.75)

# Tweek around with parasite axises
ax11sox = ax11.twinx()
ax11nox = ax11.twinx()
ax11nmvoc = ax11.twinx()

ax11sox.set_ylabel("$SO_x$ (Gg)", color='orange')
ax11sox.set_ylim(0,800)
ax11sox.tick_params(axis='y', colors='orange')
ax11sox.spines["right"].set_edgecolor('orange')

ax11nox.set_ylabel("$NO_x$ (Gg)", color='blue')
ax11nox.set_ylim(0,1600)
ax11nox.tick_params(axis='y', colors='blue')
ax11nox.spines['right'].set_position(("axes", 1.1))
ax11nox.spines["right"].set_edgecolor('blue')

ax11nmvoc.set_ylabel("$NMVOC$ (Gg)", color='red')
ax11nmvoc.set_ylim(0,16)
ax11nmvoc.tick_params(axis='y', colors='red')
ax11nmvoc.spines['right'].set_position(("axes", 1.2))
ax11nmvoc.spines["right"].set_edgecolor('red')

# Change the color of the grid lines
#for each in ax11sox.get_ygridlines():
#    each.set_color('orange')
#for each in ax11nox.get_ygridlines():
#    each.set_color('blue')
#for each in ax11nmvoc.get_ygridlines():
#    each.set_color('red')  

line_cycle = ('-','--','-.', custom_linestyle['s-'], custom_linestyle['-..'])
color_cycle = ('black', 'orange', 'blue', 'red')

# Plot the data
for ocean,i in zip(oceans,range(len(line_cycle))):
    data_ocean_sel = data.where(data['# Format: ISO2']==ocean).dropna()
    if i < 3:
        data_ocean_sel.where(data_ocean_sel['POLLUTANT']=='CO').dropna().plot('YEAR', 'NUMBER/FLAG', color=color_cycle[0], ls=line_cycle[i], ax=ax11, label=oceans[ocean])
        data_ocean_sel.where(data_ocean_sel['POLLUTANT']=='SOx').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11sox, color=color_cycle[1], ls=line_cycle[i], label='_')
        data_ocean_sel.where(data_ocean_sel['POLLUTANT']=='NOx').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11nox, color=color_cycle[2], ls=line_cycle[i], label='_')
        data_ocean_sel.where(data_ocean_sel['POLLUTANT']=='NMVOC').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11nmvoc, color=color_cycle[3], ls=line_cycle[i], label='_')
    else:
        data_ocean_sel.where(data_ocean_sel['POLLUTANT']=='CO').dropna().plot('YEAR', 'NUMBER/FLAG', color=color_cycle[0], dashes=line_cycle[i], ax=ax11, label=oceans[ocean])
        data_ocean_sel.where(data_ocean_sel['POLLUTANT']=='SOx').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11sox, color=color_cycle[1], dashes=line_cycle[i], label='_')
        data_ocean_sel.where(data_ocean_sel['POLLUTANT']=='NOx').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11nox, color=color_cycle[2], dashes=line_cycle[i], label='_')
        data_ocean_sel.where(data_ocean_sel['POLLUTANT']=='NMVOC').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11nmvoc, color=color_cycle[3], dashes=line_cycle[i], label='_')

ax11.set_xlabel("Time (year)")
ax11.set_ylabel("$CO$ (Gg)")
ax11.set_ylim(0,160)

# Show it
plt.show(block=False)



