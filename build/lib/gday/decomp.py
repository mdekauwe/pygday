""" Soil and water decomposition rates """

import math

from water_balance import WaterBalance
from utilities import float_eq, float_lt, float_le, float_gt, float_ge, clip

__author__  = "Martin De Kauwe"
__version__ = "1.0 (25.02.2011)"
__email__   = "mdekauwe@gmail.com"


class DecompFactors(object):
    """ Calculate C and N litter production rates """
    def __init__(self, control, params, state, fluxes, met_data):
        """
        Parameters
        ----------
        control : integers, structure
            model control flags
        params: floats, structure
            model parameters
        state: floats, structure
            model state
        fluxes : floats, structure
            model fluxes
        met_data : floats, dictionary
            meteorological forcing data

        """
        self.params = params
        self.fluxes = fluxes
        self.control = control
        self.state = state
        self.met_data = met_data

        self.wb = WaterBalance(self.control, self.params, self.state,
                                self.fluxes, self.met_data)

    def decay_rates(self, project_day):
        """ Model decay rates - temperature dependency (i.e. increase with temp)
        [See section A8 in Comins and McMurtrie 1993].

        Parameters:
        -----------
        project_day : int
            current simulation day (index)

        """
        # temperature and water factors for decomposition
        tempact = self.soil_temp_factor(project_day)
        wtfac = self.wb.calculate_soil_water_fac(topsoil=True)

        # decay rate of surface structural pool
        self.params.decayrate[0] = (self.params.kdec1 *
                                        math.exp(-3. * self.params.ligshoot) *
                                        tempact * wtfac)

        # decay rate of surface metabolic pool
        self.params.decayrate[1] = self.params.kdec2 * tempact * wtfac


        # decay rate of soil structural pool
        self.params.decayrate[2] = (self.params.kdec3 *
                                        math.exp(-3. * self.params.ligroot) *
                                        tempact * wtfac)

        # decay rate of soil metabolic pool
        self.params.decayrate[3] = self.params.kdec4 * tempact * wtfac


        # decay rate of active pool
        self.params.decayrate[4] = (self.params.kdec5 *
                                        (1.0 - 0.75 * self.params.finesoil) *
                                        tempact * wtfac)

        # decay rate of slow pool
        self.params.decayrate[5] = self.params.kdec6 * tempact * wtfac

        # decay rate of passive pool
        self.params.decayrate[6] = self.params.kdec7 * tempact * wtfac

    def soil_temp_factor(self, project_day):
        """Soil-temperature activity factor (A9).

        Parameters:
        -----------
        project_day : int
            current simulation day (index)

        Returns:
        --------
        tfac : float
            soil temperature factor [degC]

        """
        tsoil = self.met_data['tsoil'][project_day]

        if float_gt(tsoil, 0.0):
            tfac = (0.0326 + 0.00351 * tsoil**1.652 - (tsoil / 41.748)**7.19)
            if float_lt(tfac, 0.0):
                tfac = 0.0
        else:
            # negative number cannot be raised to a fractional power
            # number would need to be complex
            tfac = 0.0

        return tfac
