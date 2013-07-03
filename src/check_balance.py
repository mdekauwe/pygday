#!/usr/bin/env python
""" Check model C, N and water balances """

import sys
from math import fabs, exp, sqrt, sin, pi, log
import constants as const
from utilities import float_eq, float_lt, float_le, float_gt, day_length

__author__  = "Martin De Kauwe"
__version__ = "1.0 (01.07.2013)"
__email__   = "mdekauwe@gmail.com"


class CheckBalance(object):
    """ Check the model is balancing C, N and water
    
    - Haven't decided if it makes sense to run this all the time or as a check
      when changes are made to the code base.
    """
    def __init__(self, control, params, state, fluxes, met_data):
        """
        Parameters
        ----------
        control : integers, object
            model control flags
        params: floats, object
            model parameters
        state: floats, object
            model state
        fluxes : floats, object
            model fluxes
        met_data : floats, dictionary
            meteorological forcing data

        """
        self.params = params
        self.fluxes = fluxes
        self.control = control
        self.state = state
        self.met_data = met_data
        
    def check_water_balance(self, project_day, tolerance=1E-4):
        
        sources = self.met_data['rain'][project_day]
        sinks = (self.fluxes.runoff + self.fluxes.transpiration + 
                 self.fluxes.soil_evap + self.fluxes.interception)
        stores = self.state.delta_sw_store
        balance = sources - sinks - stores
        if fabs(balance) > tolerance:
            raise ValueError("Water balance check error on project day: %d" % project_day)
        