import numpy as np
import pandas as pd
import xarray as xr


def growing_season(temperature,**kwargs):
    '''
    Classical definition of growing season:
    5 consecutive days above 5 degree Celsius.
    '''
    s_shift = kwargs.pop('s_shift',365/2)
    degree_days_crit = kwargs.pop('ddc', 5)
    temperature_crit = kwargs.pop('tc', 5)
    sh = kwargs.pop('sh', False)
    verbose = kwargs.pop('verbose', False)
    count_gdd = 0
    count_days = 0
    start_gs = (0,False)
    end_gs = (0,False)
    if sh:
        temp = temperature.roll(time=s_shift)
    else:
        temp = temperature
    for itemp in temp:
        count_days += 1
        if not start_gs[1]:
            if itemp > temperature_crit:
                count_gdd += 1
            else:
                count_gdd = 0
            if count_gdd == degree_days_crit:
                if sh:
                    start_gs = (count_days+s_shift, True)
                else:
                    start_gs = (count_days, True)
                count_gdd = 0
        elif start_gs[1] and not end_gs[1] and count_days>s_shift:
            if itemp <= temperature_crit:
                count_gdd += 1
                if verbose:
                    print(count_days, count_gdd)
            else:
                count_gdd = 0
            if count_gdd == degree_days_crit:
                if sh:
                    end_gs = (count_days-s_shift, True)
                else:
                    end_gs = (count_days, True)
                count_gdd = 0
       
    return({'sgs':start_gs,'egs':end_gs})

def growing_season_sigmoid(temperature, **kwargs):
    degree_days_crit = kwargs.pop('ddc', 5)
    temperature_crit = kwargs.pop('tc', 5)
    precision = kwargs.pop('precision',1e-5)
    verbose = kwargs.pop('verbose', False)
    count_days = 0
    sgs = False
    growing = []
    for itemp in temperature:
        if itemp <= 0 and not sgs:
            count_days = max(0, count_days-1)
        if itemp > temperature_crit:
            count_days += 1
       
        growing.append(1/(1+np.exp(-2*(count_days-degree_days_crit))))
        if ((not sgs) and (np.isclose(growing[-1],[1],rtol=precision))):
            sgs = True
            if verbose:
                print(len(growing))
    return np.array(growing)

def falling_season_sigmoid(temperature, **kwargs):
    degree_days_crit = kwargs.pop('ddc', 5)
    temperature_crit = kwargs.pop('tc', 5)
    s_shift = kwargs.pop('s_shift',182)
    sgs = kwargs.pop('sgs',True)
    count_days = 0
    falling = []
    if not sgs:
        return(0)
    for itemp in temperature:
        if itemp <= temperature_crit and len(falling)>s_shift:
            count_days += 1
        falling.append(np.exp(-2*(count_days-degree_days_crit))/(1+np.exp(-2*(count_days-degree_days_crit))))    
    return np.array(falling)

def start_growing_season_fixed(lat, **kwargs):
    bLeap = kwargs.pop('leap', False)
    if (bLeap):
       MAY31 = 152
       NOV1  = 306
       DEC31 = 366
    else:
       MAY31 = 151
       NOV1  = 305
       DEC31 = 365
    
    if lat < -65:
        return((0,))
    elif lat < -23:
        return((NOV1,1))
    elif lat <= 23:
        return((1,))
    elif lat < 65:
        return(((lat-23)*4.5,))
    else:
        return((0,))
    
def end_growing_season_fixed(lat, **kwargs):
    bLeap = kwargs.pop('leap', False)
    if (bLeap):
       MAY31 = 152
       NOV1  = 306
       DEC31 = 366
    else:
       MAY31 = 151
       NOV1  = 305
       DEC31 = 365
       
    if lat < -65:
        return((0,))
    elif lat < -23:
        return((DEC31,MAY31))
    elif lat <= 23:
        return((DEC31,))
    elif lat < 65:
        return((DEC31-(lat-23)*3.3,))
    else:
        return((0,))

def growing_season_moving(temperature, **kwargs):
    degree_days_crit = kwargs.pop('ddc', 5)
    temperature_crit = kwargs.pop('tc', 5)
    verbose = kwargs.pop('verbose', False)
    msum_temp = (temperature-temperature_crit).rolling(time=degree_days_crit).sum()
    index = np.where(msum_temp>degree_days_crit*temperature_crit)
    sgs = index[0][0]
    egs = index[0][max(np.where(index[0][1:]-index[0][:-1]>1)[0])]
    if verbose:
        print(np.where(index[0][1:]-index[0][:-1]>1)[0])
    #return((sgs, egs))
    return((temperature).where(temperature>temperature_crit).cumsum(dim='time'))

def growing_season_stadyn(xr_temp):
    '''
    Preprocess the greening season for OsloCTM3.
    '''
    # Standard prescribed greeening
    bLeap = False
    if xr_temp.time.size==366:
       bLeap = True
       
    if (xr_temp.lat > 45 and xr_temp.lat < 85):
        gs = growing_season(xr_temp)
        sgs = (gs['sgs'][0],)
        egs = (gs['egs'][0],)
        if ( not gs['sgs'][1] or not gs['egs'][1] or sgs[0]>=egs[0]):
            # Fall back to fixed
            sgs = start_growing_season_fixed(xr_temp.lat.data, leap=bLeap)
            egs = end_growing_season_fixed(xr_temp.lat.data, leap=bLeap)
    elif (xr_temp.lat < -35 and xr_temp.lat >= -65):
        gs = growing_season(xr_temp, sh=True)
        sgs = (gs['sgs'][0],1)
        egs = (xr_temp.size, gs['egs'][0])
        if (not gs['sgs'][1] or not gs['egs'][1] or sgs[0]>=egs[0]):
            sgs = start_growing_season_fixed(xr_temp.lat.data, leap=bLeap)
            egs = end_growing_season_fixed(xr_temp.lat.data, leap=bLeap)
    else:
        sgs = start_growing_season_fixed(xr_temp.lat.data, leap=bLeap)
        egs = end_growing_season_fixed(xr_temp.lat.data, leap=bLeap)   
    if (sgs[0]==egs[0]==0):
        gday = np.repeat(0,len(xr_temp.time))
        glen = 0
    else:
        if len(sgs)==1:
            gday = np.concatenate((np.repeat(0,int(sgs[0])-1),
                                   np.arange(1,int(egs[0])-int(sgs[0])+1),
                                   np.repeat(0,len(xr_temp.time)-int(egs[0])+1)))
            glen = int(egs[0])-int(sgs[0])+1
        else:
            gday = np.concatenate((np.arange(1,int(egs[1])+1)+len(xr_temp.time)-int(sgs[0])+1,
                                   np.repeat(0,int(sgs[0])-int(egs[1])),
                                   np.arange(1,len(xr_temp.time)-int(sgs[0])+1)))
            glen = len(xr_temp.time)-(int(sgs[0])-int(egs[1]))+1
    if not(len(gday)==len(xr_temp.time)):
        print(xr_temp, sgs, egs)

    return(xr.DataArray(gday.astype(int),[('time', xr_temp.time.data)]), (int)(glen))
