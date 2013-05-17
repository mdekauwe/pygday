""" Soil C and N flows into 4 litter pools (structural and metabolic, both
above and belowground) and 3 SOM pools (Active, slow and passive). In
essence the CENTURY model.

Active pool -> soil microbes and microbial products, turnover time of mths-yrs.
Slow pool -> resistant plant material, turnover time of 20-50 yrs.
Passive pool -> very resistant to decomp, turnover time of > 400 yrs.
"""

from math import exp

import constants as const
from utilities import float_eq, float_lt, float_le, float_gt, float_ge

__author__  = "Martin De Kauwe"
__version__ = "1.0 (22.04.2013)"
__email__   = "mdekauwe@gmail.com"


class CarbonSoilFlows(object):
    """ Plant litter C production is divided btw metabolic and structural """
    
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

    def calculate_csoil_flows(self, project_day):
        """ C from decomposing litter -> active, slow and passive SOM pools.

        Parameters:
        -----------
        project_day : integer
            simulation day

        """
        # calculate model decay rates
        self.calculate_decay_rates(project_day)

        # plant litter inputs to the metabolic and structural pools determined 
        # by ratio of lignin/N ratio 
        (lnleaf, lnroot) = self.ligin_nratio()
        self.params.fmleaf = self.metafract(lnleaf)
        self.params.fmroot = self.metafract(lnroot)
       
        self.flux_from_grazers() # input from faeces
        self.cflux_from_plants()
        self.cfluxes_from_struct_pool()
        self.cfluxes_from_metabolic_pool()
        self.cfluxes_from_passive_pool()
        
        # update the C pools
        self.calculate_cpools()
        
    def calculate_decay_rates(self, project_day):
        """ Model decay rates - decomposition rates have a strong temperature 
        and moisture dependency. Note same temperature is assumed for all 3 
        SOM pools, found by Knorr et al (2005) to be untrue. N mineralisation
        depends on top soil moisture (most variable) (Connell et al. 1995)
        
        References:
        -----------
        Knorr et al. (2005) Nature, 433, 298-301.
        Connell et al. (1995) Biol. Fert. Soils, 20, 213-220.
    
        Parameters:
        -----------
        project_day : int
            current simulation day (index)

        """
        # temperature factor for decomposition
        tempact = self.soil_temp_factor(project_day)
        rate_scalar = self.state.wtfac_tsoil * tempact
        
        # decay rate of surface structural pool
        self.params.decayrate[0] = (self.params.kdec1 *
                                    exp(-3. * self.params.ligshoot) *
                                    rate_scalar)
        
        # decay rate of surface metabolic pool
        self.params.decayrate[1] = self.params.kdec2 * rate_scalar

        # decay rate of soil structural pool
        self.params.decayrate[2] = (self.params.kdec3 * 
                                    exp(-3.0 * self.params.ligroot) * 
                                    rate_scalar)

        # decay rate of soil metabolic pool
        self.params.decayrate[3] = self.params.kdec4 * rate_scalar

        # decay rate of active pool
        self.params.decayrate[4] = (self.params.kdec5 *
                                    (1.0 - 0.75 * self.params.finesoil) *
                                    rate_scalar)
                                        
        # decay rate of slow pool
        self.params.decayrate[5] = self.params.kdec6 * rate_scalar

        # decay rate of passive pool
        self.params.decayrate[6] = self.params.kdec7 * rate_scalar

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
    
    def flux_from_grazers(self):
        """ Input from faeces """
        if self.control.grazing:
            arg = (self.params.ligfaeces * self.params.faecescn /
                   self.params.cfracts)
            self.params.fmfaeces = self.metafract(arg)
            self.fluxes.faecesc = self.fluxes.ceaten * self.params.fracfaeces
        else:
            self.params.fmfaeces = 0.0
            self.fluxes.faecesc = 0.0

    def ligin_nratio(self):
        """ Estimate Lignin/N ratio, as this dictates the how plant litter is 
        seperated between metabolic and structural pools.

        Returns:
        --------
        lnleaf : float
            lignin:N ratio of leaf
        lnroot : float
            lignin:N ratio of fine root
        """
        nceleaf = self.ratio_of_litternc_to_live_leafnc()
        nceroot = self.ratio_of_litternc_to_live_rootnc()
        
        if float_eq(nceleaf, 0.0):
            #lnleaf = 1E20 # This is in the code, but why, this seems a mental thing to do???
            lnleaf = 0.0 
            
        else:
            lnleaf = self.params.ligshoot / self.params.cfracts / nceleaf
            
        if float_eq(nceroot, 0.0):
            #lnroot = 1E20 # This is in the code, but why, this seems a mental thing to do???
            lnroot = 0.0
        else:
            lnroot = self.params.ligroot / self.params.cfracts / nceroot

        return (lnleaf, lnroot)

    def ratio_of_litternc_to_live_leafnc(self):
        """ratio of litter N:C to live leaf N:C

        Returns:
        --------
        nceleaf : float
            N:C ratio of litter to foliage

        """
        if float_eq(self.fluxes.deadleaves, 0.0):
            ncleaf = 0.0
        else:
            ncleaf = self.fluxes.deadleafn / self.fluxes.deadleaves
        
        if self.control.use_eff_nc:
            nceleaf = self.params.liteffnc  * (1. - self.params.fretrans)
        else:
            nceleaf = ncleaf

        return nceleaf

    def ratio_of_litternc_to_live_rootnc(self):
        """ratio of litter N:C to live root N:C

        Returns:
        --------
        nceroot : float
            N:C ratio of litter to live root

        """
        if float_eq(self.fluxes.deadroots, 0.0):
            ncroot = 0.0
        else:
            ncroot = self.fluxes.deadrootn / self.fluxes.deadroots

        if self.control.use_eff_nc:
            nceroot = (self.params.liteffnc * self.params.ncrfac *
                      (1.0 - self.params.rretrans))
        else:
            nceroot = ncroot

        return nceroot

    def metafract(self, lig2n):
        """Partition fraction to metabolic pool.
        As a function of lignin to nitrogen ratio. First two eqn in section A7

        Parameters:
        -----------
        lig2n : float
            lignin to N ratio

        Returns:
        --------
        frac : float
            partitioned fraction to metabolic pool

        """
        frac = self.params.metfrac0 + self.params.metfrac1 * lig2n
        if float_gt(frac, 0.0):
            return frac
        else:
            return 0.0

    def cflux_from_plants(self):
        """ C flux from plants into soil/surface struct and metabolic pools """
        # -> into surface structural
        self.fluxes.cresid[0] = (self.fluxes.deadleaves *
                                 (1.0 - self.params.fmleaf) +
                                 self.fluxes.deadbranch * self.params.brabove +
                                 self.fluxes.deadstems + self.fluxes.faecesc *
                                 (1.0 - self.params.fmfaeces))

        # -> into soil structural
        self.fluxes.cresid[1] = (self.fluxes.deadroots *
                                 (1.0 - self.params.fmroot) +
                                 self.fluxes.deadbranch *
                                 (1.0 - self.params.brabove))

        # -> into metabolic surface
        self.fluxes.cresid[2] = (self.fluxes.deadleaves * self.params.fmleaf +
                                 self.fluxes.faecesc * self.params.fmfaeces)

        # -> into metabolic soil
        self.fluxes.cresid[3] = (self.fluxes.deadroots * self.params.fmroot + 
                                 self.fluxes.cprootexudate)

    def cfluxes_from_struct_pool(self):
        """C fluxes from structural pools """

        structout1 = self.state.structsurf * self.params.decayrate[0]
        structout2 = self.state.structsoil * self.params.decayrate[2]
        
        # surface -> slow
        self.fluxes.cstruct[0] = structout1 * self.params.ligshoot * 0.7

        # surface -> active
        self.fluxes.cstruct[1] = structout1 * (1. - self.params.ligshoot) * 0.55
        self.fluxes.co2_to_air[0] = (structout1 * (self.params.ligshoot * 0.3 +
                                    (1. - self.params.ligshoot) * 0.45))

        # soil -> slow
        self.fluxes.cstruct[2] = structout2 * self.params.ligroot * 0.7

        # soil -> active
        self.fluxes.cstruct[3] = structout2 * (1. - self.params.ligroot) * 0.45
        self.fluxes.co2_to_air[1] = (structout2 * (self.params.ligroot * 0.3 +
                                    (1. - self.params.ligroot) * 0.55))

    def cfluxes_from_metabolic_pool(self):
        """C fluxes from metabolic pools """

        # surface metabolic pool -> active soil pool
        self.fluxes.cmetab[0] = (self.state.metabsurf *
                                 self.params.decayrate[1] * 0.45)
        self.fluxes.co2_to_air[2] = (self.state.metabsurf *
                                     self.params.decayrate[1] * 0.55)

        # soil -> act
        self.fluxes.cmetab[1] = (self.state.metabsoil *
                                 self.params.decayrate[3] * 0.45)
        self.fluxes.co2_to_air[3] = (self.state.metabsoil *
                                     self.params.decayrate[3] * 0.55)

        # c fluxes from active pool
        self.fluxes.activelossf = 0.85 - 0.68 * self.params.finesoil
        activeout = self.state.activesoil * self.params.decayrate[4]

        # -> slow
        self.fluxes.cactive[0] = activeout * (0.996 - self.fluxes.activelossf)

        # -> passive
        self.fluxes.cactive[1] = activeout * 0.004
        self.fluxes.co2_to_air[4] = activeout * self.fluxes.activelossf

        # c fluxes from slow pool
        slowout = self.state.slowsoil * self.params.decayrate[5]

        # -> active
        self.fluxes.cslow[0] = slowout * 0.42

        # -> passive
        self.fluxes.cslow[1] = slowout * 0.03
        self.fluxes.co2_to_air[5] = slowout * 0.55


    def cfluxes_from_passive_pool(self):
        """ C fluxes from passive pool """

        # -> act
        self.fluxes.passive = (self.state.passivesoil *
                               self.params.decayrate[6] * 0.45)
        self.fluxes.co2_to_air[6] = (self.state.passivesoil *
                                     self.params.decayrate[6] * 0.55)

        # total co2 production
        self.fluxes.hetero_resp = (sum(self.fluxes.co2_to_air) + 
                                   self.fluxes.microbial_resp)


        # insert following line so value of resp obeys c conservn if fix
        # passive pool
        if self.control.passiveconst == True:
            self.fluxes.hetero_resp = (self.fluxes.hetero_resp +
                                       self.fluxes.cactive[1] +
                                       self.fluxes.cslow[1] -
                                       self.state.passivesoil *
                                       self.params.decayrate[6])
   
    def calculate_cpools(self):
        """Calculate new soil carbon pools. """
        
        # net source fluxes
        cstsu = self.fluxes.cresid[0] # s surf
        cstsl = self.fluxes.cresid[1] # s soil
        cmtsu = self.fluxes.cresid[2] # m surf
        cmtsl = self.fluxes.cresid[3] # m soil
        
        # store the C SOM fluxes for Nitrogen calculations
        self.fluxes.cact = (self.fluxes.cstruct[1] + self.fluxes.cstruct[3] +
                            sum(self.fluxes.cmetab) + self.fluxes.cslow[0] +
                            self.fluxes.passive)
        self.fluxes.cslo = (self.fluxes.cstruct[0] + self.fluxes.cstruct[2] +
                            self.fluxes.cactive[0])
        self.fluxes.cpas = (self.fluxes.cactive[1] + self.fluxes.cslow[1])
        
        # update pools
        self.state.structsurf += (cstsu - (self.fluxes.cstruct[0] +
                                  self.fluxes.cstruct[1] +
                                  self.fluxes.co2_to_air[0]))
        self.state.structsoil += (cstsl - (self.fluxes.cstruct[2] +
                                  self.fluxes.cstruct[3] +
                                  self.fluxes.co2_to_air[1]))
        
        # When nothing is being added to the metabolic pools, there is the 
        # potential scenario with the way the model works for tiny bits to be
        # removed with each timestep. Effectively with time this value which is
        # zero can end up becoming zero but to a silly decimal place
        self.state.metabsurf += (cmtsu - (self.fluxes.cmetab[0] +
                                 self.fluxes.co2_to_air[2]))
        
        self.state.metabsoil += (cmtsl - (self.fluxes.cmetab[1] +
                                 self.fluxes.co2_to_air[3]))
        self.state.activesoil += (self.fluxes.cact - (self.fluxes.cactive[0] +
                                  self.fluxes.cactive[1] +
                                  self.fluxes.co2_to_air[4]))
        self.state.slowsoil += (self.fluxes.cslo - (self.fluxes.cslow[0] +
                                self.fluxes.cslow[1] +
                                self.fluxes.co2_to_air[5]))
        self.state.passivesoil += (self.fluxes.cpas - (self.fluxes.passive +
                                   self.fluxes.co2_to_air[6]))
        self.state.carbon_loss += self.fluxes.hetero_resp
        
        
                    
class NitrogenSoilFlows(object):
    """ Calculate daily nitrogen fluxes"""
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

    def calculate_nsoil_flows(self):

        self.grazer_inputs()
        (nsurf, nsoil) = self.inputs_from_plant_litter()
        self.inputs_from_structrual_pool(nsurf, nsoil)

        # remaining N goes to metabolic pools
        self.fluxes.nresid[2] = nsurf - self.fluxes.nresid[0]
        self.fluxes.nresid[3] = nsoil - self.fluxes.nresid[1]

        # SOM nitrogen effluxes.  These are assumed to have the source n:c
        # ratio prior to the increase of N:C due to co2 evolution.
        self.nfluxes_from_structural_pools()
        self.nfluxes_from_metabolic_pool()
        self.nfluxes_from_active_pool()
        self.nfluxes_from_slow_pool()
        self.nfluxes_from_passive_pool()
        
        # gross N mineralisation 
        self.fluxes.ngross = self.calculate_nmineralisation()

        # calculate N immobilisation
        self.fluxes.nimmob = self.calculate_nimmobilisation()
        
        # Update model soil N pools
        self.calculate_npools()
        
    def grazer_inputs(self):
        """ Grazer inputs from faeces and urine, flux detd by faeces c:n """
        if self.control.grazing:
            self.params.faecesn = self.fluxes.faecesc / self.params.faecescn
        else:
            self.params.faecesn = 0.0

        # make sure faecesn <= total n input to soil from grazing
        arg = self.fluxes.neaten * self.params.fractosoil
        if float_gt(self.params.faecesn, arg):
            self.params.faecesn = self.fluxes.neaten * self.params.fractosoil

        # urine=total-faeces
        if self.control.grazing:
            self.fluxes.nurine = (self.fluxes.neaten * self.params.fractosoil -
                                  self.params.faecesn)
        else:
            self.fluxes.nurine = 0.0

        if float_lt(self.fluxes.nurine, 0.0):
            self.fluxes.nurine = 0.0

    def inputs_from_plant_litter(self):
        """ inputs from plant litter.

        surface and soil pools are independent. Structural input flux n:c can
        be either constant or a fixed fraction of metabolic input flux.

        Returns:
        --------
        nsurf : float
            N input from surface pool
        nsoil : float
            N input from soil pool

        """
        # surface and soil inputs (faeces n goes to abovgrd litter pools)
        nsurf = (self.fluxes.deadleafn + self.fluxes.deadbranchn *
                 self.params.brabove + self.fluxes.deadstemn +
                 self.params.faecesn)
        
        nsoil = (self.fluxes.deadrootn + self.fluxes.deadbranchn *
                 (1. - self.params.brabove))

        return nsurf, nsoil

    def inputs_from_structrual_pool(self, nsurf, nsoil):
        """structural pool input fluxes

        Parameters:
        -----------
        nsurf : float
            N input from surface pool
        nsoil : float
            N input from soil pool
        """
        # constant structural input n:c as per century
        
        if not self.control.strfloat:
            # dead plant -> structural

            # surface
            self.fluxes.nresid[0] = self.fluxes.cresid[0] / self.params.structcn

            # soil
            self.fluxes.nresid[1] = self.fluxes.cresid[1] / self.params.structcn

            # if not enough N for structural, all available N goes to structural
            if float_gt(self.fluxes.nresid[0], nsurf):
                self.fluxes.nresid[0] = nsurf
            if float_gt(self.fluxes.nresid[1], nsoil):
                self.fluxes.nresid[1] = nsoil
        else:
            # structural input n:c is a fraction of metabolic
            cwgtsu = (self.fluxes.cresid[0] * self.params.structrat +
                      self.fluxes.cresid[2])
            if float_eq(cwgtsu, 0.0):
                self.fluxes.nresid[0] = 0.0
            else:
                self.fluxes.nresid[0] = (nsurf * self.fluxes.cresid[0] *
                                         self.params.structrat / cwgtsu)

            cwgtsl = (self.fluxes.cresid[1] * self.params.structrat +
                      self.fluxes.cresid[3])
            if float_eq(cwgtsl, 0.0):
                self.fluxes.nresid[1] = 0.
            else:
                self.fluxes.nresid[1] = (nsurf * self.fluxes.cresid[1] *
                                         self.params.structrat / cwgtsl)

    def nfluxes_from_structural_pools(self):
        """ from structural pool """
        structout1 = self.state.structsurfn * self.params.decayrate[0]
        structout2 = self.state.structsoiln * self.params.decayrate[2]
        
        sigwt = (structout1 / (self.params.ligshoot * 0.7 +
                (1. - self.params.ligshoot) * 0.55))

        # surface -> slow
        self.fluxes.nstruct[0] = sigwt * self.params.ligshoot * 0.7

        # surface -> active
        self.fluxes.nstruct[1] = sigwt * (1. - self.params.ligshoot) * 0.55

        sigwt = (structout2 / (self.params.ligroot * 0.7 +
                (1. - self.params.ligroot) * 0.45))

        # soil -> slow
        self.fluxes.nstruct[2] = sigwt * self.params.ligroot * 0.7

        # soil -> active
        self.fluxes.nstruct[3] = sigwt * (1. - self.params.ligroot) * 0.45


    def nfluxes_from_metabolic_pool(self):
        """ N fluxes from metabolic pool"""

        # surf -> active
        self.fluxes.nmetab[0] = (self.state.metabsurfn *
                                 self.params.decayrate[1])

        # soil -> active
        self.fluxes.nmetab[1] = (self.state.metabsoiln *
                                 self.params.decayrate[3])

    def nfluxes_from_active_pool(self):
        """ from active pool """
        activeout = self.state.activesoiln * self.params.decayrate[4]
        sigwt = activeout / (1. - self.fluxes.activelossf)

        # -> slow
        self.fluxes.nactive[0] = sigwt * (1. - self.fluxes.activelossf - 0.004)

        # -> passive
        self.fluxes.nactive[1] = sigwt * 0.004


    def nfluxes_from_slow_pool(self):
        """ from slow pool """
        slowout = self.state.slowsoiln * self.params.decayrate[5]
        sigwt = slowout / 0.45

        # -> active
        self.fluxes.nslow[0] = sigwt * 0.42

        # -> passive
        self.fluxes.nslow[1] = sigwt * 0.03

    def nfluxes_from_passive_pool(self):
        """ from passive pool """

        # -> active
        self.fluxes.npassive = (self.state.passivesoiln *
                                self.params.decayrate[6])

    def calculate_nmineralisation(self):
        """ N gross mineralisation rate is given by the excess of N outflows 
        over inflows
        
        Returns:
        --------
        value : float
            Gross N mineralisation 
        """
        return  (sum(self.fluxes.nstruct) + sum(self.fluxes.nmetab) +
                 sum(self.fluxes.nactive) + sum(self.fluxes.nslow) + 
                 self.fluxes.npassive)
    
    def calculate_nimmobilisation(self):
        """ Calculated N immobilised in new soil organic matter
        
         General equation for new soil N:C ratio vs Nmin, expressed as linear 
         equation passing through point Nmin0, actnc0 (etc). Values can be 
         Nmin0=0, Actnc0=Actncmin 
         
         if Nmin < Nmincrit:
            New soil N:C = soil N:C (when Nmin=0) + slope * Nmin
         
         if Nmin > Nmincrit
            New soil N:C = max soil N:C       
        
        NB N:C ratio of new passive SOM can change even if assume Passiveconst
        
        Returns:
        --------
        nimob : float
            N immobilsed
        """
        conv = const.M2_AS_HA / const.G_AS_TONNES
        
        # N:C new SOM - active, slow and passive
        self.state.actncslope = self.calculate_nc_slope(self.params.actncmax, 
                                                        self.params.actnc0)
        self.state.slowncslope = self.calculate_nc_slope(self.params.slowncmax, 
                                                        self.params.slownc0)
        self.state.passncslope = self.calculate_nc_slope(self.params.passncmax, 
                                                        self.params.passnc0) 
        
        arg1 = ((self.fluxes.cactive[1] + self.fluxes.cslow[1]) *
                (self.params.passnc0 - self.state.passncslope * 
                self.params.nmin0 / conv))
        arg2 = ((self.fluxes.cstruct[0] + self.fluxes.cstruct[2] +
                self.fluxes.cactive[0]) *
                (self.params.slownc0 - self.state.slowncslope *
                self.params.nmin0 / conv))
        arg3 = ((self.fluxes.cstruct[1] + self.fluxes.cstruct[3] +
                sum(self.fluxes.cmetab) + self.fluxes.cslow[0] +
                self.fluxes.passive) * (self.params.actnc0 - 
                self.state.actncslope * self.params.nmin0 / conv))
        numer1 = arg1 + arg2 + arg3


        arg1 = ((self.fluxes.cactive[1] + self.fluxes.cslow[1]) *
                self.params.passncmax)
        arg2 = ((self.fluxes.cstruct[0] + self.fluxes.cstruct[2] +
                self.fluxes.cactive[0]) * self.params.slowncmax)
        arg3 = ((self.fluxes.cstruct[1] + self.fluxes.cstruct[3] +
                sum(self.fluxes.cmetab) + self.fluxes.cslow[0] +
                self.fluxes.passive) * self.params.actncmax)
        numer2 = arg1 + arg2 + arg3

        arg1 = ((self.fluxes.cactive[1] + self.fluxes.cslow[1]) * 
                self.state.passncslope)
        arg2 = ((self.fluxes.cstruct[0] + self.fluxes.cstruct[2] +
                self.fluxes.cactive[0]) * self.state.slowncslope)
        arg3 = ((self.fluxes.cstruct[1] + self.fluxes.cstruct[3] +
                sum(self.fluxes.cmetab) + self.fluxes.cslow[0] +
                self.fluxes.passive) * self.state.actncslope)
        denom = arg1 + arg2 + arg3
        
        # evaluate N immobilisation in new SOM
        nimmob = numer1 + denom * self.state.inorgn
        if float_gt(nimmob, numer2):
            nimmob = numer2
        
        return nimmob
 
    def calculate_nc_slope(self, pool_ncmax, pool_ncmin):
        """ Returns N:C ratio of the mineral pool slope """
        arg1 = (pool_ncmax - pool_ncmin)
        arg2 = (self.params.nmincrit - self.params.nmin0) 
        arg3 = const.M2_AS_HA / const.G_AS_TONNES
        
        return arg1 / arg2 * arg3 #slope
    
    def calculate_npools(self):
        """ Calculate new soil N pools. """

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

        # net N release implied by separation of litter into structural
        # & metabolic. The following pools only fix or release N at their 
        # limiting n:c values. 
        
        # N released or fixed from the N inorganic pool is incremented with
        # each call to nclimit and stored in self.fluxes.nlittrelease
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
        
        # When nothing is being added to the metabolic pools, there is the 
        # potential scenario with the way the model works for tiny bits to be
        # removed with each timestep. Effectively with time this value which is
        # zero can end up becoming zero but to a silly decimal place
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
        fixn = ncflux(self.fluxes.cact, nact, actnc)
        self.state.activesoiln += nact + fixn - lact

        # slow
        slownc = self.params.slownc0 + self.state.slowncslope * arg
        if float_gt(slownc, self.params.slowncmax):
            slownc = self.params.slowncmax
        fixn = ncflux(self.fluxes.cslo, nslo, slownc)
        self.state.slowsoiln += nslo + fixn - lslo

        # passive
        passnc = self.params.passnc0 + self.state.passncslope * arg
        if float_gt(passnc, self.params.passncmax):
            passnc = self.params.passncmax
        fixn = ncflux(self.fluxes.cpas, npas, passnc)
        # update passive pool only if passiveconst=0
        self.state.passivesoiln += npas + fixn - lpas

        # Daily increment of soil inorganic N pool, diff btw in and effluxes
        # (grazer urine n goes directly into inorganic pool) nb inorgn may be
        # unstable if rateuptake is large
        self.state.inorgn += ((self.fluxes.ngross + self.fluxes.ninflow + 
                               self.fluxes.nrootexudate + self.fluxes.nurine - 
                               self.fluxes.nimmob - self.fluxes.nloss - 
                               self.fluxes.nuptake) + self.fluxes.nlittrelease)
        
        
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


def ncflux(cflux, nflux, nc_ratio):
    """Returns the amount of N fixed

    Release N to Inorgn or fix N from Inorgn, in order to normalise
    the N: C ratio of a net flux.
    """
    return cflux * nc_ratio - nflux
