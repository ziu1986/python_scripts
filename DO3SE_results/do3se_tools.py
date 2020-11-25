import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_data(src):
    # Read as Excel object
    data = pd.ExcelFile(src)
    return(data)

def uncumsum(date, var):
    test = date[var].values
    # Undo cum sum and return new pandas series
    return(pd.Series(np.append(test[0], (test[2:]-test[1:-1])), index=date.index[:-1]))

def biomass(uptake, **karg):
    species = karg.pop('species', 'Norway spruce')
    b_above_gr = karg.pop('above_ground', False)

    biomass_tot = {'Birch': (100.2, 0.93), 'Norway spruce': (99.8, 0.22), "Perennial grass": (94.7, 0.62)}
    biomass_above_gr = {'Perennial grass': (93.9, 0.99)}


    if (b_above_gr) & (species == 'Perennial grass'):
        return(np.ones_like(uptake)*biomass_above_gr[species][0]-biomass_above_gr[species][1]*uptake)
    else:
        return(np.ones_like(uptake)*biomass[species][0]-biomass[species][1]*uptake)

def correlations(date, var, **karg):
    '''
    Compute correlation for each veriable in DO3SE model.
    
    '''
    b_uncum = karg.pop('uncum', False)
    b_daily = karg.pop('daily', False)
    b_verbose = karg.pop('verbose', False)
    keys = karg.pop('keys', date.keys())
    output = {}
    if b_uncum:
        date_var = uncumsum(date, var)
    else:
        date_var = date[var]

    for each in keys:
        if b_daily:
            tmp = date_var.groupby(date['Day']).mean().corr(date.groupby('Day').mean()[each])
        else:
            tmp = date_var.corr(date[each])
            
        if ~np.isnan(tmp):
            output.update({each: tmp})
        if b_verbose:
            print(each, tmp)
    return(output)

def char_range(c1, c2):
    """Generates the characters from `c1` to `c2`, inclusive."""
    for c in xrange(ord(c1), ord(c2)+1):
        yield chr(c)

def plot_pody_gsto_o3(ax, date, **karg):
    o3_color = karg.pop('o3color', 'blueviolet')
    start_day = karg.pop('start_day', 90)
    stop_day = karg.pop('stop_day', 300)
    # Define additional axis
    ax2 = ax.twinx()
    ax3 = ax.twinx()

    date['Gsto_l (mmol/m^2/s)'].plot(ax=ax, ls='None', marker='o', label='$G_{sto}^{leaf}$')
    date['O3_zR (ppb)'].plot(ax=ax2, ls='None', marker='x', color=o3_color, alpha=0.5, label='$[O_3]$')
    date['PODY (mmol/m^2 PLA)'].plot(ax=ax3, ls='None', marker='s', color='darkgreen', label='$POD_y$')
    
    # Adjust axis
    ax.set_ylim(0,250)
    ax.set_ylabel('$G_{sto}^{leaf}$ ($mmol\,m^{-2}s^{-1})$')
    
    ax2.set_ylim(0,250)
    ax2.set_yticklabels("")
    ax2.set_yticks(())
    
    ax3.set_ylim(0,16)
    ax3.set_ylabel('$POD_y$ ($mmol\,m^{-2}$ PLA)')
    ax3.grid(b=False)

    ax.set_xlabel("Time (doy)")
    ax.set_xlim(start_day*24,stop_day*24)
    ax.set_xticks(np.arange(start_day*24,stop_day*24,30*24))
    ax.set_xticklabels(np.arange(start_day*24,stop_day*24,30*24)/24)

    # Legend
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    lines3, labels3 = ax3.get_legend_handles_labels()
    ax.legend(lines + lines2 + lines3, labels + labels2 + labels3, loc="upper left")
