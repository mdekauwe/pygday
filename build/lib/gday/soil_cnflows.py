""" Soil C and N flows into 4 litter pools (structural and metabolic, both
above and belowground) and 3 SOM pools (Active, slow and passive). In
essence the CENTURY model.

Active pool -> soil microbes and microbial products, turnover time of mths-yrs.
Slow pool -> resistant plant material, turnover time of 20-50 yrs.
Passive pool -> very resistant to decomp, turnover time of > 400 yrs.
"""

import math

import constants as const
from utilities import float_eq, float_lt, float_le, float_gt, float_ge

__author__  = "Martin De Kauwe"
__version__ = "1.0 (05.09.2011)"
__email__   = "mdekauwe@gmail.com"


class CarbonFlows(object):
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

    def calculate_cflows(self, project_day):
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
       
        # input from faeces
        self.flux_from_grazers()

        self.cflux_from_plants()
        self.cfluxes_from_struct_pool()
        self.cfluxes_from_metabolic_pool()
        self.cfluxes_from_passive_pool()
    
    
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
        
        # decay rate of surface structural pool
        self.params.decayrate[0] = (self.params.kdec1 *
                                        math.exp(-3. * self.params.ligshoot) *
                                        tempact * self.state.wtfac_tsoil)

        # decay rate of surface metabolic pool
        self.params.decayrate[1] = (self.params.kdec2 * tempact * 
                                        self.state.wtfac_tsoil)


        # decay rate of soil structural pool
        self.params.decayrate[2] = (self.params.kdec3 *
                                        math.exp(-3. * self.params.ligroot) *
                                        tempact * self.state.wtfac_tsoil)

        # decay rate of soil metabolic pool
        self.params.decayrate[3] = (self.params.kdec4 * tempact * 
                                        self.state.wtfac_tsoil)

        # decay rate of active pool
        self.params.decayrate[4] = (self.params.kdec5 *
                                        (1.0 - 0.75 * self.params.finesoil) *
                                        tempact * self.state.wtfac_tsoil)
                                        
        # decay rate of slow pool
        self.params.decayrate[5] = (self.params.kdec6 * tempact * 
                                        self.state.wtfac_tsoil)

        # decay rate of passive pool
        self.params.decayrate[6] = (self.params.kdec7 * tempact * 
                                        self.state.wtfac_tsoil)

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
            lnleaf = 1E20
        else:
            lnleaf = self.params.ligshoot / self.params.cfracts / nceleaf
            #print lnleaf, self.params.ligshoot, self.params.cfracts, nceleaf
        if float_eq(nceroot, 0.0):
            lnroot = 1E20
        else:
            lnroot = self.params.ligroot / self.params.cfracts / nceroot

        return lnleaf, lnroot

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
        #print frac, self.params.metfrac0, self.params.metfrac1, lig2n
        if float_gt(frac, 0.0):
            return frac
        else:
            return 0.0

    def cflux_from_plants(self):
        """ C flux from plants into soil/surface struct and metabolic pools """
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
        self.fluxes.cresid[3] = self.fluxes.deadroots * self.params.fmroot


        # switch off production flows
        if self.control.sel_noprod1:
            self.fluxes.cresid[0] = 0.
        if self.control.sel_noprod2:
            self.fluxes.cresid[1] = 0.
        if self.control.sel_noprod3:
            self.fluxes.cresid[2] = 0.
        if self.control.sel_noprod4:
            self.fluxes.cresid[3] = 0.

    def cfluxes_from_struct_pool(self):
        """C fluxes from structural pools """

        structout = self.state.structsurf * self.params.decayrate[0]

        # surface -> slow
        self.fluxes.cstruct[0] = structout * self.params.ligshoot * 0.7

        # surface -> active
        self.fluxes.cstruct[1] = structout * (1. - self.params.ligshoot) * 0.55
        self.fluxes.co2_to_air[0] = (structout * (self.params.ligshoot * 0.3 +
                                    (1. - self.params.ligshoot) * 0.45))

        structout = self.state.structsoil * self.params.decayrate[2]

        # soil -> slow
        self.fluxes.cstruct[2] = structout * self.params.ligroot * 0.7

        # soil -> active
        self.fluxes.cstruct[3] = structout * (1. - self.params.ligroot) * 0.45
        self.fluxes.co2_to_air[1] = (structout * (self.params.ligroot * 0.3 +
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

        #print self.state.activesoil , self.params.decayrate[4]


    def cfluxes_from_passive_pool(self):
        """ C fluxes from passive pool """

        # -> act
        self.fluxes.passive = (self.state.passivesoil *
                                    self.params.decayrate[6] * 0.45)
        self.fluxes.co2_to_air[6] = (self.state.passivesoil *
                                    self.params.decayrate[6] * 0.55)

        # total co2 production
        self.fluxes.hetero_resp = sum(self.fluxes.co2_to_air)


        # insert following line so value of resp obeys c conservn if fix
        # passive pool
        if self.control.passiveconst != 0:
            self.fluxes.hetero_resp = (self.fluxes.hetero_resp +
                                            self.fluxes.cactive[1] +
                                            self.fluxes.cslow[1] -
                                            self.state.passivesoil *
                                            self.params.decayrate[6])


class NitrogenFlows(object):
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
        self.conv = const.M2_AS_HA / const.G_AS_TONNES

    def calculate_nflows(self):

        self.grazer_inputs()
        (nsurf, nsoil) = self.inputs_from_plant_litter()
        self.inputs_from_structrual_pool(nsurf, nsoil)

        # remaining n goes to metabolic pools
        self.fluxes.nresid[2] = nsurf - self.fluxes.nresid[0]
        self.fluxes.nresid[3] = nsoil - self.fluxes.nresid[1]

        # SOM nitrogen effluxes.  These are assumed to have the source n:c
        # ratio prior to the increase of n:c due to co2 evolution.
        self.nfluxes_from_structural_pools()
        self.nfluxes_from_metabolic_pool()
        self.nfluxes_from_active_pool()
        self.nfluxes_from_slow_pool()
        self.nfluxes_from_passive_pool()
        
        # gross N mineralisation 
        self.fluxes.ngross = self.calculate_nmineralisation()

        # calculate N immobilisation
        self.fluxes.nimmob = self.calculate_nimmobilisation()
        
    def grazer_inputs(self):
        """ Grazer inputs from faeces and urine, flux detd by faeces c:n """
        if self.control.grazing:
            self.params.faecesn = self.fluxes.faecesc / self.params.faecescn
        else:
            self.params.faecesn = 0.0

        #make sure faecesn <= total n input to soil from grazing
        arg = self.fluxes.neaten * self.params.fractosoil
        if float_gt(self.params.faecesn, arg):
            self.params.faecesn = self.fluxes.neaten * self.params.fractosoil

        #urine=total-faeces
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

            # if not enough N for structural, all available goes to structural
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
        structout = self.state.structsurfn * self.params.decayrate[0]
        sigwt = (structout / (self.params.ligshoot * 0.7 +
                    (1. - self.params.ligshoot) * 0.55))

        # surface -> slow
        self.fluxes.nstruct[0] = sigwt * self.params.ligshoot * 0.7

        # surface -> active
        self.fluxes.nstruct[1] = sigwt * (1. - self.params.ligshoot) * 0.55

        structout = self.state.structsoiln * self.params.decayrate[2]
        sigwt = (structout / (self.params.ligroot * 0.7 +
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
        # N:C new SOM - active, slow and passive
        self.calculate_ncratio_slope_of_mineral_pools()
        
        arg1 = ((self.fluxes.cactive[1] + self.fluxes.cslow[1]) *
                    (self.params.passnc0 - self.state.passncslope * 
                    self.params.nmin0 / self.conv))
        arg2 = ((self.fluxes.cstruct[0] + self.fluxes.cstruct[2] +
                    self.fluxes.cactive[0]) *
                    (self.params.slownc0 - self.state.slowncslope *
                    self.params.nmin0 / self.conv))
        arg3 = ((self.fluxes.cstruct[1] + self.fluxes.cstruct[3] +
                    sum(self.fluxes.cmetab) + self.fluxes.cslow[0] +
                    self.fluxes.passive) * (self.params.actnc0 - 
                    self.state.actncslope * self.params.nmin0 / self.conv))
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

    def calculate_ncratio_slope_of_mineral_pools(self):
        """ N:C ratio of the slope 3 minerals (active, slow, passive) pools
        
        General equation for new soil N:C ratio vs Nmin, expressed as linear 
        equation passing through point Nmin0, actnc0 (etc). Values can be 
        Nmin0=0, Actnc0=Actncmin 
         
        if Nmin < Nmincrit:
            New soil N:C = soil N:C (when Nmin=0) + slope * Nmin
         
        if Nmin > Nmincrit
            New soil N:C = max soil N:C       
        
        NB N:C ratio of new passive SOM can change even if assume Passiveconst
        
        """
        # N:C new SOM - active, slow and passive
        arg = (self.params.nmincrit - self.params.nmin0) / self.conv
    
        self.state.actncslope = ((self.params.actncmax - self.params.actnc0) / 
                                    arg) 
        self.state.slowncslope = ((self.params.slowncmax - self.params.slownc0)/ 
                                    arg) 
        self.state.passncslope = ((self.params.passncmax - self.params.passnc0)/  
                                    arg) 
        