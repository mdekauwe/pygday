""" Carbon production module, call photosynthesis model """

import math 

import constants as const
from utilities import float_eq, float_lt, float_gt, Bunch
from bewdy import Bewdy
from plant_production_mcmurtrie import PlantProdModel
from water_balance import WaterBalance, WaterLimitedNPP
from mate import Mate
from misc_funcs import day_length

__author__  = "Martin De Kauwe"
__version__ = "1.0 (23.02.2011)"
__email__   = "mdekauwe@gmail.com"


class PlantGrowth(object): 
    """ G'DAY plant growth module. 
    
    Calls photosynthesis model, water balance and evolve plant state.
    Pools recieve C through allocation of accumulated photosynthate and N
    from both soil uptake and retranslocation within the plant.
    
    Key feedback through soil N mineralisation and plant N uptake
    
    * Note met_forcing is an object with radiation, temp and precip data
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
        self.bw = Bewdy(self.control, self.params, self.state, self.fluxes, 
                            self.met_data)
        self.wb = WaterBalance(self.control, self.params, self.state, 
                                self.fluxes, self.met_data)
        self.pp = PlantProdModel(self.control, self.params, self.state, 
                                    self.fluxes, self.met_data)
        self.wl = WaterLimitedNPP(self.control, self.params, self.state, 
                                    self.fluxes)
        
        self.mt = Mate(self.control, self.params, self.state, self.fluxes, 
                            self.met_data)
        
    def calculate_net_c_prodn(self, day, date, fdecay, rdecay):
        """Evolve plant state"
        
        Parameters:
        -----------
        day : intefer
            simulation day
        date : date string object
            date object string (yr/mth/day)
        fdecay : float
            foliage decay rate
        rdecay : float
            fine root decay rate
        """
        daylen = day_length(date, self.params.latitude)
        
        # calculate NPP
        self.carbon_production(date, day, daylen) 
        
        # calculate water balance and adjust C production for any water stress.
        # If we are using the MATE model then water stress is applied directly
        # through the Ci:Ca reln, so do not apply any scalar to production.
        if self.control.water_model == 1:
            self.wb.calculate_water_balance(day, daylen)
            # adjust carbon production for water limitations, all models except 
            # MATE!
            if self.control.model_number != 7: 
                self.wl.adjust_cproduction(self.control.water_model)
        
        # leaf N:C as a fraction of Ncmaxyoung, i.e. the max N:C ratio of 
        # foliage in young stand
        nitfac = min(1.0, self.state.shootnc / self.params.ncmaxfyoung)	
        
        # figure out allocation fractions for C
        allocfrac = self.allocate_carbon(nitfac)
        
        # Distribute new C and N through the system
        wood_nc = self.calculate_ncwood_ratios(nitfac)
        self.nitrogen_distribution(wood_nc, fdecay, rdecay, allocfrac)
        self.carbon_distribution(allocfrac, nitfac)
        self.update_plant_state(fdecay, rdecay)
    
    def calculate_ncwood_ratios(self, nitfac):
        """ Estimate the N:C ratio in the branch and stem. Option to vary 
        the N:C ratio of the stem following Jeffreys (1999) or keep it a fixed
        fraction
        
        Parameters:
        -----------
        nitfac : float
            leaf N:C as a fraction of the max N:C ratio of foliage in young 
            stand
        
        Returns:
        --------
        wood_nc : object
            object containing N:C ratio of branch, immobile and mobile stem
            components
        
        References:
        ----------
        * Jeffreys, M. P. (1999) Dynamics of stemwood nitrogen in Pinus radiata
          with modelled implications for forest productivity under elevated
          atmospheric carbon dioxide. PhD.
        """
        # n:c ratio of new branch wood
        ncbnew = (self.params.ncbnew + nitfac * 
                    (self.params.ncbnew_crit - self.params.ncbnew))
        
        # fixed N:C in the stemwood
        if self.control.fixed_stem_nc == 1:
            # n:c ratio of stemwood - immobile pool and new ring
            ncwimm = (self.params.ncwimm + nitfac * 
                        (self.params.ncwimm_crit - self.params.ncwimm))
            
            # New stem ring N:C at critical leaf N:C (mobile)
            ncwnew = (self.params.ncwnew + nitfac * 
                        (self.params.ncwnew_crit - self.params.ncwnew))
        
        # vary stem N:C based on reln with foliage, see Jeffreys.
        else:
            ncwimm = (0.0282 * self.state.shootnc + 0.000234) * self.params.fhw
        
            # New stem ring N:C at critical leaf N:C (mobile)
            ncwnew = 0.162 * self.state.shootnc - 0.00143
        
        # group stuff just to reduce func args passed through code
        return Bunch(ncbnew=ncbnew, ncwimm=ncwimm, ncwnew=ncwnew)
    
    def carbon_production(self, date, day, daylen):
        """ Calculate GPP, NPP and plant respiration 
        
        Parameters:
        -----------
        day : intefer
            simulation day
        date : date string object
            date object string (yr/mth/day)   
        daylen : float
            daytime length (hrs)
            
        References:
        -----------
        * Jackson, J. E. and Palmer, J. W. (1981) Annals of Botany, 47, 561-565.
        """
        
        # leaf nitrogen content
        self.state.ncontent = (self.state.shootnc * self.params.cfracts / 
                                self.state.sla * const.KG_AS_G)
        
        # fractional ground cover.
        if float_lt(self.state.lai, self.params.lai_cover):
            frac_gcover = self.state.lai / self.params.lai_cover
        else:
            frac_gcover = 1.0
        
        # Radiance intercepted by the canopy, accounting for partial closure
        # Jackson and Palmer (1981), derived from beer's law 
        self.state.light_interception = ((1.0 - math.exp(-self.params.kext * 
                                            self.state.lai / frac_gcover)) * 
                                            frac_gcover)
        
        # Calculate the soil moisture availability factors [0,1] in the topsoil
        # and the entire root zone
        self.state.wtfac_root = self.wb.calculate_soil_water_fac(topsoil=False)
        
        
        # Estimate photosynthesis using an empirical model
        if self.control.model_number >=0 and self.control.model_number <= 4:
            self.pp.calculate_photosynthesis(day)
        # Estimate photosynthesis using the mechanistic BEWDY model
        elif self.control.model_number >=5 and self.control.model_number <= 6:
            # calculate plant C uptake using bewdy
            self.bw.calculate_photosynthesis(frac_gcover, date, day, daylen)    
        # Estimate photosynthesis using the mechanistic MATE model. Also need to
        # calculate a water availability scalar to determine Ci:Ca reln.
        elif self.control.model_number ==7:
            
            self.mt.calculate_photosynthesis(day, daylen)
        else:
            raise AttributeError('Unknown assimilation model')
    
    def allocate_carbon(self, nitfac):
        """Carbon allocation fractions.
        Allocations to foliage tends to decrease with stand age and wood stock 
        increases. In stressed (soil/nutrient) regions fine root allocations 
        increases.
        
        Parameters:
        -----------
        nitfac : float
            leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0) 
        
        Returns:
        --------
        allocfrac : object, float
            allocation frac for leaf, root, branch and stem
        
        """
       
        alleaf = (self.params.callocf + nitfac * 
                    (self.params.callocf_crit - self.params.callocf))
        alroot = (self.params.callocr + nitfac * 
                    (self.params.callocr_crit - self.params.callocr))
        
        albranch = (self.params.callocb + nitfac * 
                    (self.params.callocb_crit - self.params.callocb))
        
        alstem = 1.0 - alleaf - alroot - albranch
        
        # group stuff just to reduce func args
        allocfrac = Bunch(alleaf=alleaf, alroot=alroot, albranch=albranch, 
                          alstem=alstem)
        
        return allocfrac

    def nitrogen_distribution(self, wood_nc, fdecay, rdecay, allocfrac):
        """ Nitrogen distribution - allocate available N through system. 
        N is first allocated to the woody component, surplus N is then allocated
        to the shoot and roots with flexible ratios.
        
        Parameters:
        -----------
        wood_nc : float
            N:C ratio for wood
        fdecay : float
            foliage decay rate
        rdecay : float
            fine root decay rate
        allocfrac : object, float
            allocation frac for leaf, root, branch and stem    
        
        """
        
        # N retranslocated proportion from dying plant tissue
        retrans = self.nitrogen_retrans(fdecay, rdecay)    
        
        # N uptake from the soil
        if self.control.constuptake: 
            self.fluxes.nuptake = self.params.nuptakez
        else:
            if self.control.scale_nup_with_availb == 0:
                # evaluate nuptake : proportional to dynamic inorganic n pool
                self.fluxes.nuptake = self.params.rateuptake * self.state.inorgn
            else:
                # Assume N uptake depends on the rate at whioch soil mineral 
                # N is made available (self.params.Uo) and the value or root C 
                # at which 50% of the available N is taken up.
                # Dewar and McMurtrie, 1996, Tree Physiology, 16, 161-171.
                arg = (self.params.uo * self.state.inorgn * 
                        (self.state.root / (self.state.root + self.params.kr)))
                self.fluxes.nuptake = max(arg, 0.0)
        #print self.fluxes.nuptake
        # N lost from system through leaching and gaseous emissions
        # - a clear assumption using fixed rate constant across the year
        self.fluxes.nloss = self.params.rateloss * self.state.inorgn
        
        # total nitrogen to allocate (tn/ha/day)
        ntot = self.fluxes.nuptake + retrans  
        
        # allocate N to pools with fixed N:C ratios
        self.fluxes.npbranch = (self.fluxes.npp * allocfrac.albranch * 
                                    wood_nc.ncbnew)
        
        # N flux into new ring (immobile component -> structrual components)
        self.fluxes.npstemimm = (self.fluxes.npp * allocfrac.alstem * 
                                    wood_nc.ncwimm)
        
        # N flux into new ring (mobile component -> can be retrans for new
        # woody tissue)
        self.fluxes.npstemmob = (self.fluxes.npp * allocfrac.alstem * 
                                    (wood_nc.ncwnew - wood_nc.ncwimm))
        
        # If we have allocated more N than we have available - cut back N prodn
        arg = (self.fluxes.npstemimm + self.fluxes.npstemmob + 
                self.fluxes.npbranch)
        
        if float_gt(arg, ntot) and not self.control.fixleafnc:
            self.fluxes.npp *= (ntot / (self.fluxes.npstemimm + 
                                self.fluxes.npstemmob + self.fluxes.npbranch))
            self.fluxes.npbranch = (self.fluxes.npp * allocfrac.albranch * 
                                    wood_nc.ncbnew)
            self.fluxes.npstemimm = (self.fluxes.npp * allocfrac.alstem * 
                                    wood_nc.ncwimm)
            self.fluxes.npstemmob = (self.fluxes.npp * allocfrac.alstem * 
                                        (wood_nc.ncwnew - wood_nc.ncwimm))
        
        ntot -= (self.fluxes.npbranch + self.fluxes.npstemimm + 
                    self.fluxes.npstemmob) 
        
        # allocate remaining N to flexible-ratio pools
        self.fluxes.npleaf = (ntot * allocfrac.alleaf / 
                                (allocfrac.alleaf + allocfrac.alroot * 
                                self.params.ncrfac))
        self.fluxes.nproot = ntot - self.fluxes.npleaf
    
    def nitrogen_retrans(self, fdecay, rdecay):
        """ Nitrogen retranslocated from senesced plant matter.
        Constant rate of n translocated from mobile pool
        
        Parameters:
        -----------
        fdecay : float
            foliage decay rate
        rdecay : float
            fine root decay rate
        
        Returns:
        --------
        N retrans : float
            N retranslocated plant matter
            
        """
        arg1 = (self.params.fretrans * fdecay * self.state.shootn + 
                    self.params.rretrans * rdecay * self.state.rootn + 
                    self.params.bretrans * self.params.bdecay * 
                    self.state.branchn)
        arg2 = (self.params.wretrans * self.params.wdecay * 
                    self.state.stemnmob + self.params.retransmob * 
                    self.state.stemnmob)
            
        return arg1 + arg2
    
    def carbon_distribution(self, allocfrac, nitfac):
        """ C distribution - allocate available C through system 
        
        Parameters:
        -----------
        allocfrac : object, float
            allocation frac for leaf, root, branch and stem    
        nitfac : float
            leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0) 
        """
        self.fluxes.cpleaf = self.fluxes.npp * allocfrac.alleaf
        self.fluxes.cproot = self.fluxes.npp * allocfrac.alroot
        self.fluxes.cpbranch = self.fluxes.npp * allocfrac.albranch
        self.fluxes.cpstem = self.fluxes.npp * allocfrac.alstem
        
        # evaluate SLA of new foliage accounting for variation in SLA with tree 
        # and leaf age (Sands and Landsberg, 2002). Assume SLA of new foliage
        # is linearly related to leaf N:C ratio via nitfac
        sla_new = (self.params.slazero + nitfac * 
                    (self.params.slamax - self.params.slazero)) 
        
        # update leaf area [m2 m-2]
        self.state.lai += (self.fluxes.cpleaf * sla_new * const.M2_AS_HA / 
                            const.KG_AS_TONNES / self.params.cfracts - 
                            (self.fluxes.deadleaves + self.fluxes.ceaten) * 
                            self.state.lai / self.state.shoot)
        
        
    def update_plant_state(self, fdecay, rdecay):
        """ Daily change in C content 
        
        Parameters:
        -----------
        fdecay : float
            foliage decay rate
        rdecay : float
            fine root decay rate
            
        """
        
        self.state.shoot += (self.fluxes.cpleaf - self.fluxes.deadleaves - 
                                self.fluxes.ceaten)	
        self.state.root += self.fluxes.cproot - self.fluxes.deadroots
        self.state.branch += self.fluxes.cpbranch - self.fluxes.deadbranch
        self.state.stem += self.fluxes.cpstem - self.fluxes.deadstems                                         
        self.state.shootn += (self.fluxes.npleaf - fdecay * self.state.shootn - 
                                self.fluxes.neaten)
        self.state.rootn += self.fluxes.nproot - rdecay * self.state.rootn
        self.state.branchn += (self.fluxes.npbranch - self.params.bdecay * 
                                self.state.branchn)
        self.state.stemnimm += (self.fluxes.npstemimm - self.params.wdecay * 
                                self.state.stemnimm)
        self.state.stemnmob += (self.fluxes.npstemmob - self.params.wdecay * 
                                self.state.stemnmob - 
                                self.params.retransmob * self.state.stemnmob)
        self.state.stemn = self.state.stemnimm + self.state.stemnmob 
        
        # maximum leaf n:c ratio is function of stand age  
        #  - switch off age effect by setting ncmaxfyoung = ncmaxfold
        ncmaxf = (self.params.ncmaxfyoung - (self.params.ncmaxfyoung - 
                    self.params.ncmaxfold) * 
                    (self.params.age - self.params.ageyoung) / 
                    (self.params.ageold - self.params.ageyoung)) 
        
        if float_lt(ncmaxf, self.params.ncmaxfold):
            ncmaxf = self.params.ncmaxfold 
        
        if float_gt(ncmaxf, self.params.ncmaxfyoung):
            ncmaxf = self.params.ncmaxfyoung 
        
        # if foliage or root n:c ratio exceeds its max, then nitrogen uptake is 
        # cut back n.b. new ring n/c max is already set because it is related 
        # to leaf n:c
        extrar = 0.    
        extras = 0.    
        if float_gt(self.state.shootn, (self.state.shoot * ncmaxf)):
            extras = self.state.shootn - self.state.shoot * ncmaxf
            
            #n uptake cannot be reduced below zero. 
            if float_gt(extras, self.fluxes.nuptake):
                extras = self.fluxes.nuptake  
            
            self.state.shootn -= extras 
            self.fluxes.nuptake -= extras  
        
        ncmaxr = ncmaxf * self.params.ncrfac  # max root n:c 
        
        if float_gt(self.state.rootn, (self.state.root * ncmaxr)):
            extrar = self.state.rootn - self.state.root * ncmaxr
            
            #n uptake cannot be reduced below zero. 
            if float_gt((extras + extrar), self.fluxes.nuptake):
                extrar = self.fluxes.nuptake - extras 
            
            self.state.rootn -= extrar 
            self.fluxes.nuptake -= extrar #/ self.fluxes.deltay 
    
    
    
if __name__ == "__main__":
    
    from file_parser import ConfigFileParser
    import datetime
    
    # pylint: disable=C0103
    # pylint: disable=C0324
    cfg_fname = "/Users/mdekauwe/src/python/GDAY_model/params/gday.cfg"
    
    # read in user defined variables (stored in dictionaries)
    pars = ConfigFileParser(cfg_fname=cfg_fname) 
    
    (adj_control, adj_params, adj_state, adj_files, 
                                        adj_fluxes, forcing) = pars.main()
    
    adj_control.model_number = 3
    
    year = str(adj_control.startyear)
    month = str(adj_control.startmonth)
    day = str(adj_control.startday)
    datex = datetime.datetime.strptime((year + month + day), "%Y%m%d")
    
    adj_state.lai = (adj_params.slainit * const.M2_AS_HA / 
                            const.KG_AS_TONNES / adj_params.cfracts * 
                            adj_state.shoot)
    
    
   
    # Specific LAI (m2 onesided/kg DW)
    adj_state.sla = adj_params.slainit
    
    adj_control.model_number = 3
    
    # figure out photosynthesis
    PG = PlantGrowth(adj_control, adj_params, adj_state, adj_fluxes, forcing)
    
    
    
    for i in xrange(len(forcing.doy)):
        
        PG.grow(i, datex)
        print adj_fluxes.gpp / const.HA_AS_M2 * const.TONNES_AS_G
        
        
        
        # this is done in derive so do here
        # Specific LAI (m2 onesided/kg DW)
        adj_state.sla = (adj_state.lai / const.M2_AS_HA * 
                            const.KG_AS_TONNES * 
                            adj_params.cfracts / adj_state.shoot) 
        
        datex += datetime.timedelta(days=1)
    