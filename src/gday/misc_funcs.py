""" various functions used throughout """

import datetime
import calendar
import math

from utilities import float_lt, float_gt

def day_length(date, latitude):
    """ Figure out number of sunlight hours, (hours day-1)
    
    Routine from sdgvm. date is a python object, see datetime library for 
    more info 
    
    Parameters:
    -----------
    date : date format string
        date object, yr/month/day
    latitude : float    
        latitude [degrees]
        
    Returns:
    --------
    dayl : float 
        daylength [hrs]
    
    """
    conv = math.pi / 180.0
    
    # day of year 1-365/366
    doy = int(date.strftime('%j'))
    
    # Total number of days in year
    if calendar.isleap(date.year):
        yr_days = 366.
    else:
        yr_days = 365.
    
    solar_declin = -23.4 * math.cos(conv * yr_days * (doy + 10.0) / yr_days)
    temx = -math.tan(latitude * conv) * math.tan(solar_declin * conv) 
    
    if float_lt(math.fabs(temx), 1.0):
        has = math.acos(temx) / conv
        dayl = 2.0 * has / 15.0
    elif float_gt(temx, 0.0):
        dayl = 0.0
    else:
        dayl = 24.0
    
    return dayl

def mate_day_length(date, latitude):
    """ Calculate number of sunlight hours (units = h d-1)
    
    Routine comes from MATE, though not sure how right this is, are the 
    hemispheres inverted? Check
    
    Parameters:
    -----------
    date : date format string
        date object, yr/month/day
    latitude : float    
        latitude [degrees]
        
    Returns:
    --------
    dayl : float 
        daylength [hrs]
    
    """
    conv = math.pi / 180.0
    
    # day of year 1-365/366
    doy = int(date.strftime('%j'))
    
    # Total number of days in year
    if calendar.isleap(date.year):
        yr_days = 366.
    else:
        yr_days = 365.
    
    solar_dec = (23.4 * math.pi / 180.0 * math.cos(2.0 * math.pi / 
                    yr_days * (doy + 10.0)))
    
    if float_lt(latitude, 0.0):
        solar_dec *= -1.0
	 
    dayl = (math.acos(-math.tan(latitude * conv) * math.tan(solar_dec)) * 
	            24.0 / math.pi)

    return dayl
    
def clip(value, min=None, max=None):
    """ clip value btw defined range """
    if float_lt(value, min):
        value = min
    elif float_gt(value, max):
        value = max
    return value
    

 
    
    