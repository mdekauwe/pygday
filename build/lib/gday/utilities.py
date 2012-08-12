""" Various misc funcs """

from math import fabs, exp, sqrt, sin, pi, cos, tan, acos
import sys

__author__  = "Martin De Kauwe"
__version__ = "1.0 (09.03.2011)"
__email__   = "mdekauwe@gmail.com"

def float_eq(arg1, arg2, tol=1E-14):
    """arg1 == arg2"""
    return fabs(arg1 - arg2) < tol + tol * fabs(arg2)
    
def float_ne(arg1, arg2, tol=1E-14):
    """arg1 != arg2"""
    return not float_eq(arg1, arg2)

def float_lt(arg1, arg2, tol=1E-14):
    """arg1 < arg2"""
    return arg2 - arg1 > fabs(arg1) * tol
    
def float_le(arg1, arg2, tol=1E-14):
    """arg1 <= arg2"""
    return float_lt(arg1, arg2)

def float_gt(arg1, arg2, tol=1E-14):
    """arg1 > arg2"""
    return arg1 - arg2 > fabs(arg1) * tol

def float_ge(arg1, arg2, tol=1E-14):
    """arg1 >= arg2"""
    return float_gt(arg1, arg2)

def day_length(doy, yr_days, latitude):
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
    conv = pi / 180.0
    solar_declin = -23.4 * cos(conv * yr_days * (doy + 10.0) / yr_days)
    temx = -tan(latitude * conv) * tan(solar_declin * conv)

    if float_lt(fabs(temx), 1.0):
        has = acos(temx) / conv
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
    conv = pi / 180.0

    # day of year 1-365/366
    doy = int(date.strftime('%j'))

    # Total number of days in year
    if calendar.isleap(date.year):
        yr_days = 366.
    else:
        yr_days = 365.

    solar_dec = (23.4 * pi / 180.0 * cos(2.0 * pi /
                    yr_days * (doy + 10.0)))

    if float_lt(latitude, 0.0):
        solar_dec *= -1.0

    dayl = acos(-tan(latitude * conv) * tan(solar_dec)) * 24.0 / pi

    return dayl

def clip(value, min=None, max=None):
    """clip(value [, min [, max]]) => value

    Return value clipped to the range [min, max] inclusive. If either
    min or max is None, no clipping is performed on that side.
    """
    if min is not None and value < min:
        value = min
    if max is not None and value > max:
        value = max
    return value

def uniq(inlist): 
    # order preserving
    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques

def calculate_daylength(yr_days, latitude):
    daylen = []
    for d in xrange(yr_days):
        daylen.append(day_length(d+1, yr_days, latitude))   
    return daylen

if __name__ == '__main__':

    print float_eq(0.0, 0.0)
    
    print float_ge(2.1000001, 2.1000001)