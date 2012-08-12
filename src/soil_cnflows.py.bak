""" Soil C and N flows into 4 litter pools (structural and metabolic, both
above and belowground) and 3 SOM pools (Active, slow and passive). In 
essence the CENTURY model.

Active pool -> soil microbes and microbial products, turnover time of mths-yrs.
Slow pool -> resistant plant material, turnover time of 20-50 yrs.
Passive pool -> very resistant to decomp, turnover time of > 400 yrs.
"""

from decomp import DecompFactors
from utilities import float_eq, float_lt, float_le, float_gt, float_ge

__author__  = "Martin De Kauwe"
__version__ = "1.0 (25.02.2011)"
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
        
        # decomposition rates depend on soil moisture and temperature
        self.dc = DecompFactors(self.control, self.params, self.state, 
                                self.fluxes, self.met_data)
        
    def calculate_cflows(self, project_day):        
        """ C from decomposing litter -> active, slow and passive SOM pools.
        
        Parameters:
        -----------
        project_day : integer
            simulation day
        
        """
        
        # calculate model decay rates
        self.dc.decay_rates(project_day)
            
        # plant litter inputs
        lnleaf, lnroot = self.ligin_nratio()
        self.params.fmleaf = self.metafract(lnleaf)
        self.params.fmroot = self.metafract(lnroot)
        #print self.metafract(lnleaf), lnleaf
        # input from faeces
        self.flux_from_grazers()
        
        self.cflux_from_plants()
        self.cfluxes_from_struct_pool()
        self.cfluxes_from_metabolic_pool()
        self.cfluxes_from_passive_pool()
        
     
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
        """ first equation section A7, Comins and McMurtrie, 1993 
        
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
        
        # surf -> act
        
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
        
    def calculate_nflows(self):     
        
        self.grazer_inputs()
        nsurf, nsoil = self.inputs_from_plant_litter()
        self.inputs_from_structrual_pool(nsurf, nsoil)
        
        # remaining n goes to metabolic pools
        self.fluxes.nresid[2] = nsurf - self.fluxes.nresid[0]
        self.fluxes.nresid[3] = nsoil - self.fluxes.nresid[1]
        
        # SOM nitrogen effluxes.  these are assumed to have the source n:c 
        # ratio prior to the increase of n:c due to co2 evolution.
        self.nfluxes_from_structural_pools()
        self.nfluxes_from_metabolic_pool()
        self.nfluxes_from_active_pool()
        self.nfluxes_from_slow_pool()
        self.nfluxes_from_passive_pool()
    
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
        
            # if not enough n for structural, all available goes to structural
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

     