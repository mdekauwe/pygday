""" Various misc funcs """

import datetime
import calendar
import math
import sys

__author__  = "Martin De Kauwe"
__version__ = "1.0 (09.03.2011)"
__email__   = "mdekauwe@gmail.com"


def get_attrs(obj):
    """ Get user class attributes and exclude builitin attributes
    Returns a list

    Parameters:
    ----------
    obj : object
        clas object

    Returns:
    -------
    attr list : list
        list of attributes in a class object

    """
    return [i for i in obj.__dict__.keys()
                if not i.startswith('__') and not i.endswith('__')]

class Bunch(object):
    """ group a few variables together

    advantage is that it avoids the need for dictionary syntax, taken from the
    python cookbook; although nothing to guard against python reserved words

    >>> point = Bunch(datum=y, squared=y*y, coord=x)
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def float_eq(arg1, arg2, tol=1E-14):
    """arg1 == arg2"""
    return math.fabs(arg1 - arg2) < tol + tol * math.fabs(arg2)
    
def float_ne(arg1, arg2, tol=1E-14):
    """arg1 != arg2"""
    return not float_eq(arg1, arg2)

def float_lt(arg1, arg2, tol=1E-14):
    """arg1 < arg2"""
    return arg2 - arg1 > math.fabs(arg1) * tol
    
def float_le(arg1, arg2, tol=1E-14):
    """arg1 <= arg2"""
    return float_lt(arg1, arg2)

def float_gt(arg1, arg2, tol=1E-14):
    """arg1 > arg2"""
    return arg1 - arg2 > math.fabs(arg1) * tol

def float_ge(arg1, arg2, tol=1E-14):
    """arg1 >= arg2"""
    return float_gt(arg1, arg2)

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
    """clip(value [, min [, max]]) => value

    Return value clipped to the range [min, max] inclusive. If either
    min or max is None, no clipping is performed on that side.
    """
    if min is not None and value < min:
        value = min
    if max is not None and value > max:
        value = max
    return value


if __name__ == '__main__':

    print float_eq(0.0, 0.0)
    
    print float_ge(2.1000001, 2.1000001)