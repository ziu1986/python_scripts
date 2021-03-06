import os, sys, glob, calendar
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mytools.plot_tools import print_all
from mytools.clm_tools import *

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
        print("WARNING: No equilibrium found.")
        print(test)
        return(-1)
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
    final : bool
       Is this the final spin-up (post AD)?
       Standard value: False
    start_year : str
       Start of the final spin-up
       Standart value: '0000'
    stop_year : str
       End of the final spin-up.
       Standart value: '1000'
    Returns
    -------
    mean : float
       Equilibrium value of spin-up.
    '''
    final = karg.pop('final', False)
    start_year = karg.pop('start', '0001')
    stop_year = karg.pop('stop', '1000')
    if(final):
        print("Start year: %s Stop year: %s" % (start_year, stop_year))
        data = data.sel(time=slice(start_year, stop_year))
    index = max_years_spinup(data, **karg)
    if index >= 0:
        #print("Index:" % index)
        max_date = data.time.isel(time=index).values
        mean_state = data[index:].mean().values
        return(max_date, mean_state)
    else:
        print("Index: %d" % index)
        return(0,0)
    
    
# Clean up
plt.close('all')
basedir = os.environ['CESM_RUN']
case = ('test_brazil_spin-up', 
        'test_2000_brazil_spin-up_ozone', 
        'test_2000_brazil_spin-up_ozone2', 
        'test_2000_brazil_spin-up_ozone_luna', 
        'spin-up_brazil_2000', 
        'spin-up_brazil_2000_ozone',
        'spin-up_brazil_2000_ozone_luna_0',
        'spin-up_brazil_2000_ozone_luna_100',
        'spin-up_brazil_2000_ozone_luna_100_pwu', 
        'spin-up_brazil_2000_5.0.34', 
        'spin-up_brazil_2000_5.0.34_ozone', 
        'spin-up_brazil_2000_5.0.34_ozone_luna_100', 
        'spin-up_brazil_2000_5.0.34_wohydr', 
        'spin-up_brazil_2000_5.0.34_wohydr_ozone', 
        'spin-up_brazil_2000_5.0.34_wohydr_ozone_luna_100',
        'spin-up_brazil_2000_5.0.34_ozone_luna_100_sha',
        'test_merge_master_issue1224',
        'test_merge_master_issue1224_ozone_luna')
subdir1 = ('work', 'archive')
subdir2 = {'work':'run/', 'archive':'lnd/hist/'}
filename = "*.clm2.h0.*"
nyears = 20 #(20, 11)

start = '0001' #('0411','0111')
postAD = False
icase = -1
src = basedir + '/' + subdir1[1] + '/' + case[-2] + '/' + subdir2[subdir1[1]] + filename
src2 = basedir + '/' + subdir1[1] + '/' + case[-1] + '/' + subdir2[subdir1[1]] + filename

ozone_luna = False

data = load_data(src, var=['NPP','GPP','TLAI','TOTECOSYSC','TOTECOSYSN',"TOTSOMC", "TOTSOMN", "TOTVEGC", "TOTVEGN"])
data2 = load_data(src2, var=['NPP','GPP','TLAI','TOTECOSYSC','TOTECOSYSN',"TOTSOMC", "TOTSOMN", "TOTVEGC", "TOTVEGN"])
    
# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("test_spinup_%s" % (case[icase]))

ax11 = plt.subplot(231)
ax12 = plt.subplot(232)
ax13 = plt.subplot(233)
ax14 = plt.subplot(234)
ax15 = plt.subplot(235)
ax16 = plt.subplot(236)

for date, label in zip((data, data2),('ref', 'ozone_luna')):
    date['GPP'].plot(ax=ax11, label='GPP_%s' % label)
    date['NPP'].plot(ax=ax11, label='NPP_%s' % label)
    date['TLAI'].plot(ax=ax12)
    date['TOTECOSYSC'].plot(ax=ax14, label='TOTECOSYSC_%s' % label)
    date['TOTSOMC'].plot(ax=ax14, label='TOTSOMC_%s' % label)
    date['TOTVEGC'].plot(ax=ax14, label='TOTVEGC_%s' % label)
    date['TOTECOSYSN'].plot(ax=ax15, label='TOTECOSYSN_%s' % label)
    date['TOTSOMN'].plot(ax=ax15, label='TOTSOMN_%s' % label)
    date['TOTVEGN'].plot(ax=ax15, label='TOTVEGN_%s' % label)

ax11.legend()
ax14.legend()
ax15.legend()
#ax11.set_ylim()
for ax in (ax13, ax16):
    ax.set_xticklabels("")
    ax.set_yticklabels("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)


for date, ax_stats, icase in zip((data, data2), (ax13,ax16), (-2,-1)):
    ypos = 1
    ax_stats.text(0.,ypos,"%s" % (case[icase]))
    outfile = open(case[icase]+".dat", 'w')
    for each in ('GPP','NPP','TLAI','TOTECOSYSC','TOTECOSYSN',"TOTSOMC", "TOTSOMN", "TOTVEGC", "TOTVEGN"):
    
        ypos = ypos-0.1
        max_date, equi_value = equi_state(date[each],ncycle=nyears, final=postAD, start=start)
        if max_date:
            print("%s \t %s \t %s " % (each, str(max_date)[:4], equi_value), file=outfile)

            ax_stats.text(0.,ypos,"%s - %s - %s" % (each, str(max_date)[:4], equi_value))
        
            if ((each == 'GPP') or (each == 'NPP')):
                ax11.axhline(equi_value, ls=':', color='grey')
            elif each == 'TLAI':
                ax12.axhline(equi_value, ls=':', color='grey')
            elif ((each == 'TOTECOSYSC') or (each == 'TOTSOMC') or (each == 'TOTVEGC')):
                ax14.axhline(equi_value, ls=':', color='grey')
            elif ((each == 'TOTECOSYSN') or (each == 'TOTSOMN') or (each == 'TOTVEGN')):
                ax15.axhline(equi_value, ls=':', color='grey')
        else:
            print("%s \t Not in euqilibrium yet!" % (each))
            ax_stats.text(0.2,ypos,"%s - Not in euqilibrium yet!" % (each))
 
""" 
if ozone_luna:
    fig2 = plt.figure(2, figsize=(16,9))
    fig2.canvas.set_window_title("luna_test_spinup_%s" % (file[file.rfind('/')+1:file.find('h0')+2]))

    ax21 = plt.subplot(221)
    ax23 = plt.subplot(223)
    ax24 = plt.subplot(224)

    data_ozone['O3COEFJMAXSHA'].plot(ax=ax21, label='O3COEFJMAXSHA')
    data_ozone['O3COEFJMAXSUN'].plot(ax=ax21, label='O3COEFJMAXSUN')

    data_ozone['O3COEFVCMAXSHA'].plot(ax=ax23, label='O3COEFVCMAXSHA')
    data_ozone['O3COEFVCMAXSUN'].plot(ax=ax23, label='O3COEFVCMAXSUN')

    ax21.legend()
    ax23.legend()

    ax24.set_xticklabels("")
    ax24.set_yticklabels("")
    ax24.set_xticks([])
    ax24.set_yticks([])
    ax24.set_frame_on(False)

"""    
outfile.close()
        
# Show it
plt.show(block=False)
        

