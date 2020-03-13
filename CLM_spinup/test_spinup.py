import os, sys, glob, calendar
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all

def max_years_spinup(data, **karg):
    '''
    Get an estimate for the maximum of spin-up years.
    Parameters
    ----------
    data : xarray dataset
    Keyword arguments
    -----------------
    ncycle : int
       The number of forcing data years chosen for spin-up
    bound : float
       The threashold for the similarity to zero. 
       Standard value: 1e-3
    Returns
    -------
    index : int
       Index at which the first derivative of the spin-up date is roughly zero.
    '''
    ncycle = karg.pop('ncycle', 11)
    bound = karg.pop('bound', 1e-3)
    y = data.values.flatten()
    test = np.abs((np.roll(y,ncycle)-y)/y)
    try:
        stop_spinup = np.where(test < bound)[0][0]
    except IndexError:
        return(0)
    return(stop_spinup)

def equi_state(data, **karg):
    '''
    Get the asymptode - the equilibrium value - of the spin-up.
    Parameters
    ----------
    data: xarray dataset
    Keyword arguments
    -----------------
    ncycle : int
       The number of forcing data years chosen for spin-up
    bound : float
       The threashold for the similarity to zero. 
       Standard value: 1e-3
    Returns
    -------
    mean : float
       Equilibrium value of spin-up.
    '''
    index = max_years_spinup(data, **karg)
    return(data[index:].mean().values)
    
    
# Clean up
plt.close('all')
#src = os.environ['CESM_RUN']+"/work/test_brazil_spin-up/run/*.clm2.h0.*"
#src = os.environ['CESM_RUN']+"/archive/test_brazil_spin-up/lnd/hist/*.clm2.h0.*"
#nyears = 11
src = os.environ['CESM_RUN']+"/work/test_2000_brazil_spin-up_ozone/run/*.clm2.h0.*"
#src = os.environ['CESM_RUN']+"/archive/test_2000_brazil_spin-up_ozone/lnd/hist/*.clm2.h0.*"
#src = os.environ['CESM_RUN']+"/work/test_2000_brazil_spin-up/run/*.clm2.h0.*"
#src = os.environ['CESM_RUN']+"/archive/test_2000_brazil_spin-up/lnd/hist/*.clm2.h0.*"
nyears = 20
data_list = []
try:
    data
except NameError:
    for file in sorted(glob.glob(src))[:-1]:
        print(file)
        data = xr.open_dataset(file)
        data_list.append(data[['NPP','GPP','TLAI','TOTECOSYSC','TOTECOSYSN',"TOTSOMC", "TOTSOMN", "TOTVEGC", "TOTVEGN"]])
    data = xr.concat(data_list, dim='time')

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("test_spinup_%s" % (file[file.rfind('/')+1:file.find('h0')+2]))

ax11 = plt.subplot(231)
ax12 = plt.subplot(232)
#ax13 = plt.subplot(233)
ax14 = plt.subplot(234)
ax15 = plt.subplot(235)
ax16 = plt.subplot(236)

data['GPP'].plot(ax=ax11, label='GPP')
data['NPP'].plot(ax=ax11, label='NPP')
data['TLAI'].plot(ax=ax12)
data['TOTECOSYSC'].plot(ax=ax14, label='TOTECOSYSC')
data['TOTSOMC'].plot(ax=ax14, label='TOTSOMC')
data['TOTVEGC'].plot(ax=ax14, label='TOTVEGC')
data['TOTECOSYSN'].plot(ax=ax15, label='TOTECOSYSN')
data['TOTSOMN'].plot(ax=ax15, label='TOTSOMN')
data['TOTVEGN'].plot(ax=ax15, label='TOTVEGN')

ax11.legend(loc='upper right')
ax14.legend(loc='upper right')
ax15.legend(loc='upper right')
ax16.set_xticklabels("")
ax16.set_yticklabels("")
ax16.set_xticks([])
ax16.set_yticks([])
ax16.set_frame_on(False)


ypos = 1
for each in ('GPP','NPP','TLAI','TOTECOSYSC','TOTECOSYSN',"TOTSOMC", "TOTSOMN", "TOTVEGC", "TOTVEGN"):
    spinup_index = max_years_spinup(data[each],ncycle=nyears)
    ypos = ypos-0.1
    if spinup_index > 0:
        max_date = data.time.isel(time=spinup_index).values
        equi_value = equi_state(data[each],ncycle=nyears)
        print("%s \t %s \t %s " % (each, str(max_date)[:4], equi_value))
        ax16.text(0.2,ypos,"%s - %s - %s" % (each, str(max_date)[:4], equi_value))
        if ((each == 'GPP') or (each == 'NPP')):
            ax11.axhline(equi_value, ls=':', color='grey')
            ax11.text(0, equi_value, "%2.2f" % equi_value)
        elif each == 'TLAI':
            ax12.axhline(equi_value, ls=':', color='grey')
            ax12.text(0, equi_value, "%2.2f" % equi_value)
        elif ((each == 'TOTECOSYSC') or (each == 'TOTSOMC') or (each == 'TOTVEGC')):
            ax14.axhline(equi_value, ls=':', color='grey')
            ax14.text(0, equi_value, "%2.2f" % equi_value)
        elif ((each == 'TOTECOSYSN') or (each == 'TOTSOMN') or (each == 'TOTVEGN')):
            ax15.axhline(equi_value, ls=':', color='grey')
            ax15.text(0, equi_value, "%2.2f" % equi_value)
    else:
        print("%s \t Not in euqilibrium yet!" % (each))
        ax16.text(0.2,ypos,"%s - Not in euqilibrium yet!" % (each))
  
    
        
# Show it
plt.show(block=False)
        

