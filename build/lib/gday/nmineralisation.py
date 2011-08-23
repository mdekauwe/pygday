#!/usr/bin/env python
""" Nitrogen mineralisation rate is given by the excess of N out over N in """

import constants as const
from utilities import float_eq, float_lt, float_gt, Bunch


__author__  = "Martin De Kauwe"
__version__ = "1.0 (25.02.2011)"
__email__   = "mdekauwe@gmail.com"


class Mineralisation(): 
    """ Inputs from plant litter partitioned btw metabolic and structural"""
    
    def __init__(self, control, params, state, fluxes): 
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
    
        """
        self.params = params
        self.fluxes = fluxes
        self.control = control
        self.state = state
    
    def calculate_mineralisation(self):
        """Mineralisation calculations.
        Plant N uptake, soil N loss, N gross mineralisation, N immobilisation
        """ 
        # gross n mineralisation (t/ha/yr) 
        self.fluxes.ngross = (sum(self.fluxes.nstruct) + 
                                sum(self.fluxes.nmetab) +
                                sum(self.fluxes.nactive) + 
                                sum(self.fluxes.nslow) + self.fluxes.npassive) 
        
        # N:C new SOM - slope is an object containing active, slow and passive
        slope = self.calc_new_nc_slope_params()
        
        # Mineral N pool
        (numer1, numer2, denom) = self.calc_mineral_npool(slope)
        
        # evaluate n immobilisation in new SOM
        self.fluxes.nimmob = numer1 + denom * self.state.inorgn 
        if float_gt(self.fluxes.nimmob, numer2):
            self.fluxes.nimmob = numer2
            
    def calc_new_nc_slope_params(self):
        """ General eqn for new soil n:c ratio vs nmin 
        
        Expressed as linear eqn passing thru point (nmin0, actnc0)
        values can be nmin0=0, actnc0= actncmin or 
        can be eqm values from a previous run (see notes from bordeaux 0999)
        nb n:c ratio of new passive som can change even if assume passiveconst
        
        Returns:
        --------
        slope : object, float
            object containting N:C ratio for active, slow and passive soils.
        
        """
        actnc = self.calc_ratio(self.params.actncmax, self.params.actnc0)
        slownc = self.calc_ratio(self.params.slowncmax, self.params.slownc0)
        passnc = self.calc_ratio(self.params.passncmax, self.params.passnc0)
        
        # group the return variables
        slope = Bunch(actnc=actnc, slownc=slownc, passnc=passnc)
        
        return slope
    
    
    def calc_ratio(self, argmax, argmin):
        """
        Parameters:
        -----------
        argmax : float
        
        argmin : float
            nmin of soil pool
        Returns:
        --------
            NC ratio of soil pool vs nmnin
        """
        
        numerator = argmax - argmin
        denominator = self.params.nmincrit - self.params.nmin0
        units_conv = const.M2_AS_HA / const.G_AS_TONNES

        return numerator / denominator * units_conv
        
    def calc_mineral_npool(self, slope):
        """ Need to look this up, but just numerator and demominator terms for
        the mineral N pool calculation 
        
        Parameter:
        ----------
        slope : float
            ?
        
        Returns:
        --------
        numer1 : float 
            first numerator term for mineral N pool calculation
        numer2 : float 
            second numerator term for mineral N pool calculation
        denom : float
             denominator term for mineral N pool calculation               
        """
        
        arg1 = ((self.fluxes.cactive[1] + self.fluxes.cslow[1]) * 
                    (self.params.passnc0 - slope.passnc * self.params.nmin0 / 
                    const.M2_AS_HA * const.G_AS_TONNES))
        arg2 = ((self.fluxes.cstruct[0] + self.fluxes.cstruct[2] + 
                    self.fluxes.cactive[0]) * 
                    (self.params.slownc0 - slope.slownc * 
                    self.params.nmin0 / const.M2_AS_HA * const.G_AS_TONNES))
        arg3 = ((self.fluxes.cstruct[1] + self.fluxes.cstruct[3] + 
                    sum(self.fluxes.cmetab) + self.fluxes.cslow[0] + 
                    self.fluxes.passive) * (self.params.actnc0 - slope.actnc * 
                    self.params.nmin0 / const.M2_AS_HA * const.G_AS_TONNES))
        numer1 = arg1 + arg2 + arg3 
                         
        
        arg1 = ((self.fluxes.cactive[1] + self.fluxes.cslow[1]) * 
                    self.params.passncmax)
        arg2 = ((self.fluxes.cstruct[0] + self.fluxes.cstruct[2] + 
                    self.fluxes.cactive[0]) * self.params.slowncmax)
        arg3 = ((self.fluxes.cstruct[1] + self.fluxes.cstruct[3] + 
                    sum(self.fluxes.cmetab) + self.fluxes.cslow[0] + 
                    self.fluxes.passive) * self.params.actncmax)
        numer2 = arg1 + arg2 + arg3 
                    
        arg1 = (self.fluxes.cactive[1] + self.fluxes.cslow[1]) * slope.passnc
        arg2 = ((self.fluxes.cstruct[0] + self.fluxes.cstruct[2] + 
                    self.fluxes.cactive[0]) * slope.slownc)
        arg3 = ((self.fluxes.cstruct[1] + self.fluxes.cstruct[3] + 
                    sum(self.fluxes.cmetab) + self.fluxes.cslow[0] + 
                    self.fluxes.passive) * slope.actnc)
        denom = arg1 + arg2 + arg3
        
        return numer1, numer2, denom