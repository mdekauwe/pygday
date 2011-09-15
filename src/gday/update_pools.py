""" Update model C and N pools """

import constants as const
from utilities import float_eq, float_ne, float_lt, float_le, float_gt, float_ge

__author__  = "Martin De Kauwe"
__version__ = "1.0 (25.02.2011)"
__email__   = "mdekauwe@gmail.com"


def ncflux(cflux, nflux, nc_ratio):
    """Returns the amount of N fixed

    Release N to Inorgn or fix N from Inorgn, in order to normalise
    the N: C ratio of a net flux.
    """
    return cflux * nc_ratio - nflux

class CarbonPools(object):
    """ Update model C pools"""

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

    def calculate_cpools(self):
        """Calculate new soil carbon pools.

        Returns:
        --------
        cact : float
            C source flux from the active pool
        cslo : float
            C source flux from the slow pool
        cpas : float
            C source flux from the passive pool
        """

        # net source fluxes
        cstsu = self.fluxes.cresid[0] # s surf
        cstsl = self.fluxes.cresid[1] # s soil
        cmtsu = self.fluxes.cresid[2] # m surf
        cmtsl = self.fluxes.cresid[3] # m soil
        cact = ((self.fluxes.cstruct[1] + self.fluxes.cstruct[3] +
                    sum(self.fluxes.cmetab) + self.fluxes.cslow[0] +
                     self.fluxes.passive))
        cslo = (self.fluxes.cstruct[0] + self.fluxes.cstruct[2] +
                    self.fluxes.cactive[0])
        cpas = (self.fluxes.cactive[1] + self.fluxes.cslow[1] )


        # update pools
        self.state.structsurf += (cstsu - (self.fluxes.cstruct[0] +
                                    self.fluxes.cstruct[1] +
                                    self.fluxes.co2_to_air[0]))
        self.state.structsoil += (cstsl - (self.fluxes.cstruct[2] +
                                    self.fluxes.cstruct[3] +
                                    self.fluxes.co2_to_air[1]))
        self.state.metabsurf += (cmtsu - (self.fluxes.cmetab[0] +
                                    self.fluxes.co2_to_air[2]))
        self.state.metabsoil += (cmtsl - (self.fluxes.cmetab[1] +
                                    self.fluxes.co2_to_air[3]))
        self.state.activesoil += (cact - (self.fluxes.cactive[0] +
                                    self.fluxes.cactive[1] +
                                    self.fluxes.co2_to_air[4]))
        self.state.slowsoil += (cslo - (self.fluxes.cslow[0] +
                                    self.fluxes.cslow[1] +
                                    self.fluxes.co2_to_air[5]))
        self.state.passivesoil += (cpas - (self.fluxes.passive +
                                    self.fluxes.co2_to_air[6]))
        self.state.carbon_loss += self.fluxes.hetero_resp

        return cact, cslo, cpas

class NitrogenPools(object):
    """ Update model N pools"""
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


    def calculate_npools(self, cact, cslo, cpas, day):
        """Calculate new soil N pools.

        Parameters:
        -----------
        cact : float
            C source flux from the active pool
        cslo : float
            C source flux from the slow pool
        cpas : float
            C source flux from the passive pool
        day : integer
            project day of simulation
        """
        # net source fluxes.
        nstsu = self.fluxes.nresid[0]  # s surf
        nstsl = self.fluxes.nresid[1]  # s soil
        nmtsu = self.fluxes.nresid[2]  # m surf
        nmtsl = self.fluxes.nresid[3]  # m soil
        nact = (self.fluxes.nstruct[1] + self.fluxes.nstruct[3] +
                    self.fluxes.nmetab[0] + self.fluxes.nmetab[1] +
                    self.fluxes.nslow[0] + self.fluxes.npassive)
        nslo = (self.fluxes.nstruct[0] + self.fluxes.nstruct[2] +
                    self.fluxes.nactive[0])
        npas = self.fluxes.nactive[1] + self.fluxes.nslow[1]

        # net effluxes.
        lstsu = (self.fluxes.nstruct[0] + self.fluxes.nstruct[1])   # s surf
        lstsl = (self.fluxes.nstruct[2] + self.fluxes.nstruct[3])   # s soil
        lmtsu = self.fluxes.nmetab[0]                               # m surf
        lmtsl = self.fluxes.nmetab[1]                               # m soil
        lact = (self.fluxes.nactive[0] + self.fluxes.nactive[1])
        lslo = (self.fluxes.nslow[0] + self.fluxes.nslow[1])
        lpas = self.fluxes.npassive

        # net n release implied by separation of litter into structural
        # & metabolic
        # n:c of struct litter = 1/150
        # n:c of metab litter = 1/25 to 1/10
        # update pools.
        #
        # the following pools only fix or release n at their limiting n:c
        # values. Conformity with the limiting range is achieved in a single
        # step.

        self.fluxes.nlittrelease = 0.0
        self.state.structsurfn += nstsu - lstsu
        if not self.control.strfloat:
            self.state.structsurfn += self.nclimit(self.state.structsurf,
                                                    self.state.structsurfn,
                                                    1.0/self.params.structcn,
                                                    1.0/self.params.structcn)
        self.state.structsoiln += nstsl - lstsl
        if not self.control.strfloat:
            self.state.structsoiln += self.nclimit(self.state.structsoil,
                                                    self.state.structsoiln,
                                                    1.0/self.params.structcn,
                                                    1.0/self.params.structcn)
        self.state.metabsurfn += nmtsu - lmtsu
        self.state.metabsurfn += self.nclimit(self.state.metabsurf,
                                                self.state.metabsurfn,
                                                1.0/25.0, 1.0/10.0)

        self.state.metabsoiln += nmtsl - lmtsl
        self.state.metabsoiln += self.nclimit(self.state.metabsoil,
                                                self.state.metabsoiln,
                                                1.0/25.0, 1.0/10.0)
        
        # N:C of the SOM pools increases linearly btw prescribed min and max 
        # values as the Nconc of the soil increases.
        arg = (self.state.inorgn - self.params.nmin0 / const.M2_AS_HA * 
                const.G_AS_TONNES)
        # active
        actnc = self.params.actnc0 + self.state.actncslope * arg
        if float_gt(actnc, self.params.actncmax):
            actnc = self.params.actncmax

        fixn = ncflux(cact, nact, actnc)
        self.state.activesoiln += nact + fixn - lact

        # slow
        slownc = self.params.slownc0 + self.state.slowncslope * arg
        if float_gt(slownc, self.params.slowncmax):
            slownc = self.params.slowncmax

        fixn = ncflux(cslo, nslo, slownc)
        self.state.slowsoiln += nslo + fixn - lslo

        # passive
        passnc = self.params.passnc0 + self.state.passncslope * arg
        if float_gt(passnc, self.params.passncmax):
            passnc = self.params.passncmax

        fixn = ncflux(cpas, npas, passnc)

        #update passive pool only if passiveconst=0
        self.state.passivesoiln += npas + fixn - lpas

        # Daily increment of soil inorganic N pool, diff btw in and effluxes
        # (grazer urine n goes directly into inorganic pool) nb inorgn may be
        # unstable if rateuptake is large
        self.state.inorgn += ((self.fluxes.ngrossmin + self.fluxes.ninflow +
                                self.fluxes.nurine - self.fluxes.nimmob -
                                self.fluxes.nloss - self.fluxes.nuptake) +
                                self.fluxes.nlittrelease)


    def nclimit(self, cpool, npool, ncmin, ncmax):
        """ Release N to 'Inorgn' pool or fix N from 'Inorgn', in order to keep
        the  N:C ratio of a litter pool within the range 'ncmin' to 'ncmax'.

        Parameters:
        -----------
        cpool : float
            various C pool (state)
        npool : float
            various N pool (state)
        ncmin : float
            maximum N:C ratio
        ncmax : float
            minimum N:C ratio

        Returns:
        --------
        fix/rel : float
            amount of N to be added/released from the inorganic pool

        """
        nmax = cpool * ncmax
        nmin = cpool * ncmin

        if float_gt(npool, nmax):  #release
            rel = npool - nmax
            self.fluxes.nlittrelease += rel
            return -rel
        elif float_lt(npool, nmin):   #fix
            fix = nmin - npool
            self.fluxes.nlittrelease -= fix
            return fix
        else:
            return 0.0
