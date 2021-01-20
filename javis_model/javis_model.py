def f_temp(temperature, **karg):
    '''
    Compute Javis f_temp function.
    Parameters
    ----------
    temperature : float
        The 2m temperature in deg C
    Keyword arguments
    -----------------
    Tmin : float
        Minimum temperature for stomatal opening in deg C
    Tmax : float
        Maximum temperature for stomatal opening in deg C
    Topt : float
        Optimal temperature for stomatal opening in deg C
    Returns
    -------
    f_temp : float, range 0,1
    '''
    import numpy as np
       
    TMIN = karg.pop("Tmin")
    TMAX = karg.pop("Tmax")
    TOPT = karg.pop("Topt")

    beta = ((TMAX-TOPT)/(TOPT-TMIN))
    f_temp = (temperature-TMIN)/(TOPT-TMIN)*((TMAX-temperature)/(TMAX-TOPT))**beta

    # Catch negative values
    f_temp[np.where(f_temp<0)] = 0
    # Catch values larger 1
    f_temp[np.where(f_temp>1)] = 1
    
    return(f_temp)

def f_vpd(VPD, **karg):
    '''
    Compute Javis f_vpd function.
    Parameters
    ----------
    VPD : float
        Water vapor pressure deficit
    Keyword arguments
    -----------------
    fmin : float
        Minimum value for stomatal opening
    Dmin : float
        Maximum for stomatal opening
    Dmax : float
        Maximum for stomatal opening
    Returns
    -------
    f_vpd : float, range 0,1
    
    '''
    import numpy as np
    
    FMIN = karg.pop("fmin")
    DMIN = karg.pop("Dmin")
    DMAX = karg.pop("Dmax")

    f_vpd = FMIN + (1-FMIN) *(DMIN-VPD)/(DMIN-DMAX)

    # Catch negative values
    f_vpd[np.where(f_vpd<0)] = 0
    # Catch values larger 1
    f_vpd[np.where(f_vpd>1)] = 1
    
    return(f_vpd)


def f_light(ppfd, **karg):
    '''
    Compute Javis f_light function.
    Parameters
    ----------
    ppfd : float
        Photoactive photonflaux density in W/m^2/s
    Keyword arguments
    -----------------
    alpha : float
        slope parameter in m^2*s/W
    Returns
    -------
    f_light : float, range 0,1
    
    '''

    import numpy as np
    
    ALPHA = karg.pop("alpha")

    f_light = 1 - np.exp(-ALPHA*ppfd)

    # Catch negative values
    f_light[np.where(f_light<0)] = 0
    # Catch values larger 1
    f_light[np.where(f_light>1)] = 1
        
    return(f_light)


def stomatal_conductance(f_phen, f_light, f_min, f_temp, f_vpd, f_sw, gmax):
    '''
    Compute Javis stomatal conductance.
    Parameters
    ----------
    f_phen : float

    f_light : float

    f_min : float

    f_temp : float

    f_vpd : float

    f_sw : float

    Arguments
    ---------
    gmax : float
        Maximum stomatal conductance
    Returns
    -------
    gsto : float
        Stomatal conductance

    '''
    gsto = gmax * f_phen * f_light * np.max((f_min, f_temp, f_vpd, f_sw))

    return(gsto)
    




    
