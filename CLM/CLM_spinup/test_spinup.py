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
#basedir = os.environ['CESM_RUN']
basedir = os.environ['DATA'] + '/saga/'
case = ('acc_spinup_ozone_luna_off_brazil',
        'acc_spinup_ozone_luna_lombardozzi_brazil',
        'acc_spinup_ozone_luna_falk_brazil')
subdir1 = ('work', 'archive')
subdir2 = {'work':'run/', 'archive':'lnd/hist/'}
filename = "*.clm2.h0.*"
nyears = 20 #(20, 11)

start = '0001' #('0411','0111')
postAD = False
sel_case = (0,-1)
name_case = ('ref', 'lombardozzi', 'falk')

src = []
for icase in sel_case:
    src.append(os.path.join(basedir, subdir1[1], case[sel_case[icase]], subdir2[subdir1[1]], filename))

ozone_luna = False
variables = ['NPP','GPP','TLAI','TOTECOSYSC','TOTECOSYSN',"TOTSOMC", "TOTSOMN", "TOTVEGC", "TOTVEGN"]
data = []
for each in src:
    data.append(load_data(each, var=variables))
    
# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("test_spinup_%s" % (case[sel_case[1]]))

ax11 = plt.subplot(231)
ax12 = plt.subplot(232)
ax13 = plt.subplot(233)
ax14 = plt.subplot(234)
ax15 = plt.subplot(235)
ax16 = plt.subplot(236)

for date, label in zip(data, [name_case[i] for i in sel_case]):
    ax11.plot(date['GPP'].data, label='GPP_%s' % label)
    ax11.plot(date['NPP'].data, label='NPP_%s' % label)
    ax12.plot(date['TLAI'].data)
    ax14.plot(date['TOTECOSYSC'].data, label='TOTECOSYSC_%s' % label)
    ax14.plot(date['TOTSOMC'].data, label='TOTSOMC_%s' % label)
    ax14.plot(date['TOTVEGC'].data, label='TOTVEGC_%s' % label)
    ax15.plot(date['TOTECOSYSN'].data, label='TOTECOSYSN_%s' % label)
    ax15.plot(date['TOTSOMN'].data, label='TOTSOMN_%s' % label)
    ax15.plot(date['TOTVEGN'], label='TOTVEGN_%s' % label)

    ax11.set_ylabel("($%s$)" % (date['GPP'].units))
    ax12.set_ylabel("($%s$)" % (date['TLAI'].units))
    ax14.set_ylabel("($%s$)" % (date['TOTECOSYSC'].units))
    ax15.set_ylabel("($%s$)" % (date['TOTVEGC'].units))
    ax15.set_xlabel("Time since start (years)", x=-0.25)
    
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


for date, ax_stats, icase in zip(data, (ax13,ax16), sel_case):
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
plt.tight_layout()
plt.show(block=False)
        

