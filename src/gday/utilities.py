""" Logic tests for floating point numbers """

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


if __name__ == '__main__':

    print float_eq(0.0, 0.0)
    
    print float_ge(2.1000001, 2.1000001)