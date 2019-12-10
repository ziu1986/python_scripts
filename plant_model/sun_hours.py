import numpy as np
import pandas as pd

def daylength(dayOfYear, lat):
    """Computes the length of the day (the time between sunrise and
    sunset) given the day of the year and latitude of the location.
    Function uses the Brock model for the computations.
    For more information see, for example,
    Forsythe et al., "A model comparison for daylength as a
    function of latitude and day of year", Ecological Modelling,
    1995.
    Parameters
    ----------
    dayOfYear : int
        The day of the year. 1 corresponds to 1st of January
        and 365 to 31st December (on a non-leap year).
    lat : float
        Latitude of the location in degrees. Positive values
        for north and negative for south.
    Returns
    -------
    d : float
        Daylength in hours.
    """
    latInRad = np.deg2rad(lat)
    declinationOfEarth = 23.45*np.sin(np.deg2rad(360.0*(283.0+dayOfYear)/365.0))
    if -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) <= -1.0:
        return 24.0
    elif -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) >= 1.0:
        return 0.0
    else:
        hourAngle = np.rad2deg(np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))))
        return 2.0*hourAngle/15.0

def get_avg_daylenght(start_date, end_date, lat):
    doy_start = pd.DatetimeIndex((np.datetime64(start_date),)).dayofyear[0]
    doy_end = pd.DatetimeIndex((np.datetime64(end_date),)).dayofyear[0]
    print(doy_start, doy_end)
    light_hours = [daylength(each, lat) for each in np.arange(0,365)]
    mean_daylight = np.mean(light_hours[doy_start-1:doy_end-1])
    std_daylight = np.std(light_hours[doy_start-1:doy_end-1])
    return(mean_daylight, std_daylight)
