import os, glob, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt   # Plotting data
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *

# Clean-up
plt.close('all')

# Data source
src = "EMEP_shipping_emission_data.txt"
# Open and read file
data = pd.read_csv(src,';', header=1)
# Select data
oceans = data['# Format: ISO2'].unique()
atlantic = data.where(data['# Format: ISO2']=='ATL').dropna()
medi = data.where(data['# Format: ISO2']=='MED').dropna()

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot()
plt.subplots_adjust(right=0.75)

ax11sox = ax11.twinx()
ax11nox = ax11.twinx()
ax11nmvoc = ax11.twinx()

ax11sox.set_ylabel("$SO_x$ (Gg)", color='orange')
#ax11sox.set_ylim(0,500)
ax11sox.tick_params(axis='y', colors='orange')
ax11sox.spines["right"].set_edgecolor('orange')

ax11nox.set_ylabel("$NO_x$ (Gg)", color='blue')
#ax11nox.set_ylim(500,1000)
ax11nox.tick_params(axis='y', colors='blue')
ax11nox.spines['right'].set_position(("axes", 1.1))
ax11nox.spines["right"].set_edgecolor('blue')

ax11nmvoc.set_ylabel("$NMVOC$ (Gg)", color='red')
#ax11nmvoc.set_ylim(5,10)
ax11nmvoc.tick_params(axis='y', colors='red')
ax11nmvoc.spines['right'].set_position(("axes", 1.2))
ax11nmvoc.spines["right"].set_edgecolor('red')

atlantic.where(atlantic['POLLUTANT']=='CO').dropna().plot('YEAR', 'NUMBER/FLAG', color='black', ax=ax11, label='Atlantic')
atlantic.where(atlantic['POLLUTANT']=='SOx').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11sox, color='orange', label='_')
atlantic.where(atlantic['POLLUTANT']=='NOx').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11nox, color='blue', label='_')
atlantic.where(atlantic['POLLUTANT']=='NMVOC').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11nmvoc, color='red', label='_')

medi.where(medi['POLLUTANT']=='CO').dropna().plot('YEAR', 'NUMBER/FLAG', ls='--', color='black', ax=ax11, label='Mediterranean')
medi.where(medi['POLLUTANT']=='SOx').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11sox, ls='--', color='orange', label='_')
medi.where(medi['POLLUTANT']=='NOx').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11nox, ls='--', color='blue', label='_')
medi.where(medi['POLLUTANT']=='NMVOC').dropna().plot('YEAR', 'NUMBER/FLAG', ax=ax11nmvoc, ls='--', color='red', label='_')

ax11.set_xlabel("Time (year)")
ax11.set_ylabel("$CO$ (Gg)")
#ax11.set_ylim(0,100)

# Show it
plt.show(block=False)



