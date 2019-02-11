import os, glob, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # Plotting data
import cartopy as cp        # Globe projections
import cartopy.util as ccrs_util  # Add cyclic
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.constants import *     # Get physics constants
from mytools.met_tools import *
from mytools.netcdf_tools import *
from astropy.io import ascii

celsius = 273.15
# Data source
nc_src = os.environ['DATA']+'/abel/C3RUN_nitrate_corr_2005/avgsav/avgsav_*.nc'
src = os.environ['DATA']+"/EMEP_parameters.dat"
src_2 = os.environ['DATA']+"/EMEP_parameters_2.dat"
data_emep = ascii.read(src)
data_emep2 = ascii.read(src_2)
print(data_emep2)
# Read temperature data
try:
    data
except NameError:
    data_list = []
    # Open dataset
    for file in sorted(glob.glob(nc_src)):
        print(file)
        data = xr.open_dataset(file)
        # Define new time coordinates and drop the old one
        #data['time'].reset_coords(drop=True)
        data.coords['time'] = dt.datetime(data['START_TIME'][0],data['START_TIME'][1],data['START_TIME'][2])
        data_list.append(data['temperature'])
    # Concatinating the list
    data = xr.concat(data_list, dim='time')

#temperature_range = np.linspace(0,50)
temperature_range = data.isel(lev=0).sel(lat=50, lon=24, method='nearest')-celsius
doy = (1,32,60,91,121,152,182,213,244,274,305,335)
temperature_range = np.interp(np.arange(1,366),doy,temperature_range)
doy = np.arange(1,366)

def fT(T2m, Tmin, Topt, Tmax):
    if np.all((T2m > Tmin) & (T2m < Tmax)):
        fTemp = (T2m-Tmin)/(Topt-Tmin)*((Tmax-T2m)/(Tmax-Topt))**((Tmax-Topt)/(Topt-Tmin))
    else:
        fTemp = 0
    return(fTemp)
def fD(D, fmin, Dmax, Dmin):
    fDemp = fmin+(1-fmin)*(Dmin-D)/(Dmin-Dmax)
    return(fDemp)
def fSW():
    return(0)
        
def fPhen(ptype, temperature, phia, phib, phic, phid, phie, phif, phiAS, phiAE):
    if not np.any((ptype == 'GR') | (ptype == 'SNL') | (ptype == 'CF')):
        if temperature >= 5:
            if not fPhen.sgs:
                if fPhen.dcounter==5:
                    fPhen.sgs = True
                    fPhen.dsgs += 1
                    res_pheno = phia
                else:
                    fPhen.dcounter += 1
                    res_pheno = 0
            else:
                if fPhen.dsgs <= phiAS:
                    res_pheno = phia
                elif fPhen.dsgs <= phiAS+phie:
                    res_pheno = (phic-phia)/phie * fPhen.dsgs + phia
                else:
                    res_pheno = phic
                fPhen.dsgs += 1
        else:        
            if not fPhen.sgs:
                res_pheno = 0
            else:
                if fPhen.dcounter <= -phif:
                    fPhen.sgs = False
                    fPhen.dsg = 0
                    res_pheno = 0
                else:
                    if phif>0:
                        res_pheno = (phid-phic)/(-phif)*(fPhen.dcounter-5)+phic
                        if res_pheno <=0:
                            fPhen.sgs = False
                            fPhen.dsg = 0
                            res_pheno = 0
                        print res_pheno, (phid-phic)/(-phif), fPhen.dcounter-5, phic
                    else:
                        res_pheno = phid
                    fPhen.dcounter -= 1
                    fPhen.dsgs += 1
    else:
        res_pheno = phic
    return(res_pheno)

fPhen.sgs = False
fPhen.dsgs = 0
fPhen.dcounter = 0

def gsto(temperature, params, **kwargs):
    verbose = kwargs.pop('v', False)
    bphen = kwargs.pop('phen',True)
    gmax = params[1]
    if bphen:
        Phen = fPhen(params[0], temperature,float(params[3]), float(params[4]),float(params[5]),float(params[6]),
                     float(params[7]), float(params[8]), float(params[9]),float(params[10]))
    else:
        Phen = 1
    light = params[11]
    fmin = params[2]
    T = fT(temperature,float(params[12]),float(params[13]),float(params[14]))
    D = 0
    SW = fSW()
    if verbose:
        print temperature, gmax, Phen, light, fmin, T, D, SW, np.max((fmin,T,D,SW)), gmax*Phen*light*np.max((fmin,T,D,SW))
    return(gmax*Phen*light*np.max((fmin,T,D,SW)))

#station_gsto = [gsto(itemp, data_emep[2], v=True) for itemp in temperature_range]    
# Plot it
plt.close('all')
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot(121)
for itype in range(10):
    pheno = [fPhen(data_emep[itype][0],itemp,float(data_emep[itype][3]),float(data_emep[itype][4]),float(data_emep[itype][5]),
                   float(data_emep[itype][6]),float(data_emep[itype][7]), float(data_emep[itype][8]),
                   float(data_emep[itype][9]),float(data_emep[itype][10])) for itemp in temperature_range]  
    ax11.plot(doy, pheno, label=data_emep[itype][0])
ax11.set_xlabel("Time (Days)")
ax11.set_ylabel("fPhen")
ax11.legend()
ax12 = plt.subplot(122)
ax12.plot(doy, temperature_range)
ax12.set_xlabel("Time (Days)")
ax12.set_ylabel("Temperature ($^\circ$C)")
fig2 = plt.figure(2, figsize=(16,9))
ax21 = plt.subplot()
for itype in (0,1,9):
    station_gsto = [gsto(itemp, data_emep[itype]) for itemp in temperature_range]    
    ax21.plot(doy, station_gsto, label=data_emep[itype][0])
    station_gsto = [gsto(itemp, data_emep[itype], phen=False) for itemp in temperature_range]    
    ax21.plot(doy, station_gsto, label=data_emep[itype][0], ls='--', color=ax21.lines[-1].get_color())
ax21.set_xlabel("Time (Days)")
ax21.set_ylabel("g$_{sto}$ (mmol O$_3$ m$^{-2}$ s$^{-1}$)")
ax21.legend()
# Show it
plt.show(block=False)
