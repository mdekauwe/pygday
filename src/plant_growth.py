""" Estimate daily Carbon fixed and pass around the aboveground portion of the
plant. """

from math import exp
import sys
import constants as const
from utilities import float_eq, float_lt, float_gt
from bewdy import Bewdy
from water_balance import WaterBalance, SoilMoisture
from mate import Mate
from optimal_root_model import RootingDepthModel

__author__  = "Martin De Kauwe"
__version__ = "1.0 (23.02.2011)"
__email__   = "mdekauwe@gmail.com"


class PlantGrowth(object):
    """ G'DAY plant growth module.

    Calls photosynthesis model, water balance and evolve plant state.
    Pools recieve C through allocation of accumulated photosynthate and N
    from both soil uptake and retranslocation within the plant.

    Key feedback through soil N mineralisation and plant N uptake

    * Note met_forcing is an object with radiation, temp, precip data, etc.
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
        
        self.mt = Mate(self.control, self.params, self.state, self.fluxes,
                       self.met_data)
        
        self.sm = SoilMoisture(self.control, self.params, self.state, 
                               self.fluxes)
        
        self.rm = RootingDepthModel(d0=0.35, r0=0.05, top_soil_depth=0.3)
   
    def calc_day_growth(self, project_day, fdecay, rdecay, daylen, doy, days_in_yr):
        """Evolve plant state, photosynthesis, distribute N and C"

        Parameters:
        -----------
        project_day : integer
            simulation day
        fdecay : float
            foliage decay rate
        rdecay : float
            fine root decay rate
        """
        # calculate NPP
        self.carbon_production(project_day, daylen)

        # calculate water balance
        self.wb.calculate_water_balance(project_day, daylen)
        
        # leaf N:C as a fraction of Ncmaxyoung, i.e. the max N:C ratio of
        # foliage in young stand
        nitfac = min(1.0, self.state.shootnc / self.params.ncmaxfyoung)
        
        if not self.control.deciduous_model:
            # figure out allocation fractions for C for the evergreen model. For
            # the deciduous model these are calculated at the annual time step.
            self.allocate_carbon(nitfac)
        
        # Distribute new C and N through the system
        self.carbon_distribution(nitfac, doy, days_in_yr)
        
        (ncbnew, ncwimm, ncwnew) = self.calculate_ncwood_ratios(nitfac)
        self.nitrogen_distribution(ncbnew, ncwimm, ncwnew, fdecay, rdecay, doy)
        
        
        self.update_plant_state(fdecay, rdecay, project_day)
         
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
        ncbnew : float
            N:C ratio of branch
        ncwimm : float
            N:C ratio of immobile stem
        ncwnew : float
            N:C ratio of mobile stem

        References:
        ----------
        * Jeffreys, M. P. (1999) Dynamics of stemwood nitrogen in Pinus radiata
          with modelled implications for forest productivity under elevated
          atmospheric carbon dioxide. PhD.
        """
        # n:c ratio of new branch wood
        ncbnew = (self.params.ncbnew + nitfac *
                 (self.params.ncbnew - self.params.ncbnewz))
        
        # fixed N:C in the stemwood
        if self.control.fixed_stem_nc == 1:
            # n:c ratio of stemwood - immobile pool and new ring
            ncwimm = (self.params.ncwimm + nitfac *
                     (self.params.ncwimm - self.params.ncwimmz))
            
            # New stem ring N:C at critical leaf N:C (mobile)
            ncwnew = (self.params.ncwnew + nitfac *
                     (self.params.ncwnew - self.params.ncwnewz))
            
        # vary stem N:C based on reln with foliage, see Jeffreys. Jeffreys 1999
        # showed that N:C ratio of new wood increases with foliar N:C ratio,
        # modelled here based on evidence as a linear function.
        else:
            ncwimm = max(0.0, (0.0282 * self.state.shootnc + 0.000234) * 
                         self.params.fhw)

            # New stem ring N:C at critical leaf N:C (mobile)
            ncwnew = max(0.0, 0.162 * self.state.shootnc - 0.00143)
        
        return (ncbnew, ncwimm, ncwnew)

    def carbon_production(self, project_day, daylen):
        """ Calculate GPP, NPP and plant respiration

        Parameters:
        -----------
        project_day : integer
            simulation day
        daylen : float
            daytime length (hrs)

        References:
        -----------
        * Jackson, J. E. and Palmer, J. W. (1981) Annals of Botany, 47, 561-565.
        """

        # leaf nitrogen content
        if self.state.lai > 0.0:
            # Leaf N content (g m-2)                       
            self.state.ncontent = (self.state.shootnc * self.params.cfracts /
                                   self.state.sla * const.KG_AS_G)
        else:
            self.state.ncontent = 0.0
            
        # fractional ground cover.
        if float_lt(self.state.lai, self.params.lai_cover):
            frac_gcover = self.state.lai / self.params.lai_cover
        else:
            frac_gcover = 1.0

        # Radiance intercepted by the canopy, accounting for partial closure
        # Jackson and Palmer (1981), derived from beer's law
        if self.state.lai > 0.0:
            self.state.light_interception = ((1.0 - exp(-self.params.kext *
                                             self.state.lai / frac_gcover)) *
                                             frac_gcover)
        else:
            self.state.light_interception = 0.0
        
        if self.control.water_stress:
            # Calculate the soil moisture availability factors [0,1] in the 
            # topsoil and the entire root zone
            (self.state.wtfac_tsoil, 
                self.state.wtfac_root) = self.sm.calculate_soil_water_fac()
        else:
            # really this should only be a debugging option!
            self.state.wtfac_tsoil = 1.0
            self.state.wtfac_root = 1.0
            
        # Estimate photosynthesis 
        if self.control.assim_model == "BEWDY":
            self.bw.calculate_photosynthesis(frac_gcover, project_day, daylen)
        elif self.control.assim_model == "MATE":
            self.mt.calculate_photosynthesis(project_day, daylen)
        else:
            raise AttributeError('Unknown assimilation model')
    
    
    
    def allocate_carbon(self, nitfac):
        """Carbon allocation fractions to move photosynthate through the plant.
        Allocations to foliage tends to decrease with stand age and wood stock
        increases (Makela and Hari, 1986; Cannell and Dewar, 1994; 
        Magnani et al, 2000). In stressed (soil/nutrient) regions fine root 
        allocations increases.

        Parameters:
        -----------
        nitfac : float
            leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0)

        Returns:
        --------
        alleaf : float
            allocation fraction for shoot
        alroot : float
            allocation fraction for fine roots
        albranch : float
            allocation fraction for branches
        alstem : float
            allocation fraction for stem
        alroot_exudate : float
            allocation fraction for root exudate 
       
        References:
        -----------
        McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.
        """
        self.state.alleaf = (self.params.callocf + nitfac *
                            (self.params.callocf - self.params.callocfz))
    
        self.state.alroot = (self.params.callocr + nitfac *
                            (self.params.callocr - self.params.callocrz))

        self.state.albranch = (self.params.callocb + nitfac *
                              (self.params.callocb - self.params.callocbz))
        
        # Remove some of the allocation to wood and instead allocate it to
        # root exudation. Following McMurtrie et al. 2000
        self.state.alroot_exudate = self.params.callocrx
        
        # allocate remainder to stem
        self.state.alstem = (1.0 - self.state.alleaf - self.state.alroot - 
                             self.state.albranch - self.state.alroot_exudate)
    
    def initialise_deciduous_model(self):
        """ Divide up NPP based on annual allocation fractions """
        self.state.c_to_alloc_shoot = self.state.alleaf * self.state.cstore
        self.state.c_to_alloc_root = self.state.alroot * self.state.cstore
        self.state.c_to_alloc_branch = self.state.albranch * self.state.cstore
        self.state.c_to_alloc_stem = self.state.alstem * self.state.cstore
        #self.state.c_to_alloc_rootexudate = (self.state.alroot_exudate *
        #                                        self.state.cstore)

        # annual available N for allocation to leaf
        self.state.n_to_alloc_shoot = (self.state.c_to_alloc_shoot *
                                        self.state.shootnc_yr)
    
    def allocate_stored_c_and_n(self):
        """
        At the end of the year allocate everything for the coming year
        based on stores from the previous year avaliable N for allocation
        """
        self.state.c_to_alloc_shoot = self.state.alleaf * self.state.cstore
        self.state.c_to_alloc_root = self.state.alroot * self.state.cstore
        self.state.c_to_alloc_branch = self.state.albranch * self.state.cstore
        self.state.c_to_alloc_stem = self.state.alstem * self.state.cstore
        self.state.n_to_alloc_root = (min(self.state.nstore,
                                          self.state.c_to_alloc_root *
                                          self.state.rootnc))

        # constant N:C of foliage during the growing season(kgN kg-1C)
        self.state.shootnc_yr = ((self.state.nstore -
                                  self.state.n_to_alloc_root) /
                                  (self.state.c_to_alloc_shoot))
        
        # if we want to put back a floating N:C then we need to have
        # self.state.c_to_alloc_shoot + self.state.c_to_alloc_stem * some factor

        # annual available N for allocation to leaf
        self.state.n_to_alloc_shoot = (self.state.c_to_alloc_shoot *
                                        self.state.shootnc_yr)

        
    
     
    def nitrogen_distribution(self, ncbnew, ncwimm, ncwnew, fdecay, rdecay, doy):
        """ Nitrogen distribution - allocate available N through system.
        N is first allocated to the woody component, surplus N is then allocated
        to the shoot and roots with flexible ratios.
        
        References:
        -----------
        McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.
        
        Parameters:
        -----------
        ncbnew : float
            N:C ratio of branch
        ncwimm : float
            N:C ratio of immobile stem
        ncwnew : float
            N:C ratio of mobile stem
        fdecay : float
            foliage decay rate
        rdecay : float
            fine root decay rate
        """
        # N retranslocated proportion from dying plant tissue and stored within
        # the plant
        self.fluxes.retrans = self.nitrogen_retrans(fdecay, rdecay)
        self.fluxes.nuptake = self.calculate_nuptake()
        
        # Ross's Root Model.
        # NOT WORKING YET
        
        if self.control.model_optroot == True:    
            # Attempt at floating rateuptake
            #slope = (20.0 - 0.1) / ((6.0) - 0.0)
            #y = slope * self.state.root + 0.0
            #self.fluxes.nuptake = (y/365.25) * self.state.inorgn
                
            # have you got all these units conversion right, they don't appear 
            # very consistent? you have grams and kg?!
            
            # convert t ha-1 d-1 to gN m-2 yr-1
            nsupply = self.calculate_nuptake() * 365.25 * 100.0 # t ha-1->gN m-2
            
            # covnert t ha-1 to kg m-2
            rtot = self.state.root * 0.1
            
            (root_depth, 
             self.fluxes.nuptake,
             self.fluxes.rabove) = self.rm.main(rtot, nsupply, depth_guess=1.0)
            
            # covert nuptake from gN m-2 yr-1 to t ha-1 d-1
            self.fluxes.nuptake = self.fluxes.nuptake * 0.01 / 365.25
            
            # covert from kg N m-2 to t ha-1
            self.fluxes.deadroots = self.params.rdecay * self.fluxes.rabove * 10
            self.fluxes.deadrootn = (self.state.rootnc * 
                                    (1.0 - self.params.rretrans) * 
                                    self.fluxes.deadroots)
            
            print  root_depth
            #print self.fluxes.gpp*100, self.state.lai, self.state.root*100, \
            #      self.fluxes.nuptake *100.
            #print self.fluxes.gpp*100, self.state.lai, self.state.root*100, \
            #      root_depth, self.fluxes.nuptake *100.
        
        
        # N lost from system is proportional to the soil inorganic N pool, 
        # where the rate constant empirically defines gaseous and leaching 
        # losses, see McMurtrie et al. 2001.
        self.fluxes.nloss = self.params.rateloss * self.state.inorgn
    
        # total nitrogen to allocate 
        ntot = self.fluxes.nuptake + self.fluxes.retrans
        
        # N flux into root exudation, see McMurtrie et al. 2000
        self.fluxes.nrootexudate = (self.fluxes.npp * self.state.alroot_exudate * 
                                    self.params.vxfix)
        
        if self.control.deciduous_model:
            # allocate N to pools with fixed N:C ratios
            # N flux into new ring (immobile component -> structrual components)
            self.fluxes.npstemimm = self.fluxes.cpstem * ncwimm
    
            # N flux into new ring (mobile component -> can be retrans for new
            # woody tissue)
            self.fluxes.npstemmob = self.fluxes.cpstem * (ncwnew - ncwimm)
            self.fluxes.npbranch = 0.0    
            self.fluxes.nproot = self.fluxes.cproot * self.state.rootnc
            self.fluxes.npleaf = (self.fluxes.lnrate * 
                                    self.state.growing_days[doy])
        else:
            # allocate N to pools with fixed N:C ratios
            
            # N flux into new ring (immobile component -> structrual components)
            self.fluxes.npstemimm = self.fluxes.npp * self.state.alstem * ncwimm
    
            # N flux into new ring (mobile component -> can be retrans for new
            # woody tissue)
            self.fluxes.npstemmob = (self.fluxes.npp * self.state.alstem * 
                                     (ncwnew - ncwimm))
            self.fluxes.npbranch = (self.fluxes.npp * self.state.albranch * 
                                     ncbnew)
            
            # If we have allocated more N than we have available 
            #  - cut back N prodn
            arg = (self.fluxes.npstemimm + self.fluxes.npstemmob +
                   self.fluxes.npbranch )
    
            if float_gt(arg, ntot) and self.control.fixleafnc == False:
                self.fluxes.npp *= (ntot / (self.fluxes.npstemimm +
                                    self.fluxes.npstemmob + 
                                    self.fluxes.npbranch ))
                self.fluxes.npbranch = (self.fluxes.npp * self.state.albranch * 
                                        ncbnew)
                self.fluxes.npstemimm = (self.fluxes.npp * self.state.alstem * 
                                         ncwimm)
                self.fluxes.npstemmob = (self.fluxes.npp * self.state.alstem * 
                                        (ncwnew - ncwimm))
                
                
            ntot -= (self.fluxes.npbranch + self.fluxes.npstemimm +
                        self.fluxes.npstemmob)
            
            # allocate remaining N to flexible-ratio pools
            self.fluxes.npleaf = (ntot * self.state.alleaf / 
                                 (self.state.alleaf + self.state.alroot *
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
        if self.control.deciduous_model:
            arg1 = (self.fluxes.leafretransn  +
                    self.params.rretrans * rdecay * self.state.rootn +
                    self.params.bretrans * self.params.bdecay * 
                    self.state.branchn)
            arg2 = (self.params.wretrans * self.params.wdecay *
                    self.state.stemnmob + self.params.retransmob *
                    self.state.stemnmob)
        else:
            arg1 = (self.params.fretrans * fdecay * self.state.shootn +
                    self.params.rretrans * rdecay * self.state.rootn +
                    self.params.bretrans * self.params.bdecay *
                    self.state.branchn)
            arg2 = (self.params.wretrans * self.params.wdecay *
                    self.state.stemnmob + self.params.retransmob *
                    self.state.stemnmob)
        
        return arg1 + arg2
    
    def calculate_nuptake(self):
        """ N uptake from the soil, note as it stands root biomass does not
        affect N uptake.
        
        Returns:
        --------
        nuptake : float
            N uptake
            
        References:
        -----------
        * Dewar and McMurtrie, 1996, Tree Physiology, 16, 161-171.    
            
        """
        if self.control.nuptake_model == 0:
            # Constant N uptake
            nuptake = self.params.nuptakez
        elif self.control.nuptake_model == 1:
            # evaluate nuptake : proportional to dynamic inorganic N pool
            nuptake = self.params.rateuptake * self.state.inorgn
        elif self.control.nuptake_model == 2:
            # Assume N uptake depends on the rate at which soil mineral
            # N is made available (self.params.Uo) and the value or root C
            # at which 50% of the available N is taken up (Dewar and McM).
            arg = (self.params.uo * self.state.inorgn *
                    (self.state.root / (self.state.root + self.params.kr)))
            nuptake = max(arg, 0.0)
        else:
            raise AttributeError('Unknown N uptake assumption')
        
        return nuptake
    
    def carbon_distribution(self, nitfac, doy, days_in_yr):
        """ C distribution - allocate available C through system

        Parameters:
        -----------
        nitfac : float
            leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0)
        
        References:
        -----------
        
        * Hale, M. G. et al. (1981) Factors affecting root exudation and 
          significance for the rhizosphere ecoystems. Biological and chemical 
          interactions in the rhizosphere. Stockholm. Sweden: Ecological 
          Research Committee of NFR. pg. 43--71.
        * Lambers, J. T. and Poot, P. (2003) Structure and Functioning of 
          Cluster Roots and Plant Responses to Phosphate Deficiency.
        * Martin, J. K. and Puckjeridge, D. W. (1982) Carbon flow through the 
          rhizosphere of wheat crops in South Australia. The cyclcing of carbon,
          nitrogen, sulpher and phosphorous in terrestrial and aquatic
          ecosystems. Canberra: Australian Academy of Science. pg 77--82.
        * McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.
        
        Also see:
        * Rovira, A. D. (1969) Plant Root Exudates. Botanical Review, 35, 
          pg 35--57.
        """
        if self.control.deciduous_model:
            self.fluxes.cpleaf = self.fluxes.lrate * self.state.growing_days[doy]
            self.fluxes.cpbranch = 0.0
            self.fluxes.cpstem = self.fluxes.wrate * self.state.growing_days[doy]
            self.fluxes.cproot = self.state.c_to_alloc_root * 1.0 / days_in_yr
        else:
            self.fluxes.cpleaf = self.fluxes.npp * self.state.alleaf
            self.fluxes.cproot = self.fluxes.npp * self.state.alroot
            self.fluxes.cpbranch = self.fluxes.npp * self.state.albranch
            self.fluxes.cpstem = self.fluxes.npp * self.state.alstem
        
        # C flux into root exudation, see McMurtrie et al. 2000. There is no 
        # reference given for the 0.15 in McM, however 14c work by Hale et al and
        # Martin and Puckeridge suggest values range between 10-20% of NPP. So
        # presumably this is where this values of 0.15 (i.e. the average) comes
        # from
        self.fluxes.cprootexudate = self.fluxes.npp * self.state.alroot_exudate
        #self.fluxes.cprootexudate = 0.0
        
        
        # evaluate SLA of new foliage accounting for variation in SLA 
        # with tree and leaf age (Sands and Landsberg, 2002). Assume 
        # SLA of new foliage is linearly related to leaf N:C ratio 
        # via nitfac
        sla_new = (self.params.slazero + nitfac *
                  (self.params.slamax - self.params.slazero))
        sla_new_tonnes_ha_C = (sla_new * const.M2_AS_HA / 
                             (const.KG_AS_TONNES * self.params.cfracts))
        
        
        if self.control.deciduous_model:
            if self.state.shoot == 0.0:
                self.state.lai = 0.0
            elif self.state.leaf_out_days[doy] > 0.0:
                self.state.lai += (self.fluxes.cpleaf * sla_new_tonnes_ha_C -
                                  (self.fluxes.deadleaves + self.fluxes.ceaten)*
                                   self.state.lai / self.state.shoot) 
            else:
                self.state.lai = 0.0
        else:
            # update leaf area [m2 m-2]
            self.state.lai += (self.fluxes.cpleaf * sla_new_tonnes_ha_C -
                               (self.fluxes.deadleaves + self.fluxes.ceaten) *
                                self.state.lai / self.state.shoot)
            
            
    def update_plant_state(self, fdecay, rdecay, project_day):
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
        
        if self.control.deciduous_model:
            self.state.shootn += (self.fluxes.npleaf - 
                                 (self.fluxes.deadleafn - self.fluxes.neaten))
        else:
            self.state.shootn += (self.fluxes.npleaf - 
                                  fdecay * self.state.shootn - 
                                  self.fluxes.neaten)
        self.state.shootn = max(0.0, self.state.shootn)
        
        self.state.rootn += self.fluxes.nproot - rdecay * self.state.rootn
        self.state.branchn += (self.fluxes.npbranch - self.params.bdecay *
                                self.state.branchn)
        self.state.stemnimm += (self.fluxes.npstemimm - self.params.wdecay *
                                self.state.stemnimm)
        self.state.stemnmob += (self.fluxes.npstemmob - self.params.wdecay *
                                self.state.stemnmob -
                                self.params.retransmob * self.state.stemnmob)
        self.state.stemn = self.state.stemnimm + self.state.stemnmob
        
        self.state.exu_pool += self.fluxes.cprootexudate
        self.fluxes.microbial_resp = self.calc_microbial_resp(project_day)
        self.state.exu_pool -= self.fluxes.microbial_resp
        #print self.fluxes.microbial_resp, self.fluxes.cprootexudate
        
        self.calculate_cn_store()
        
        if not self.control.deciduous_model:
            
            # maximum leaf n:c ratio is function of stand age
            #  - switch off age effect by setting ncmaxfyoung = ncmaxfold
            age_effect = ((self.state.age - self.params.ageyoung) / 
                            (self.params.ageold - self.params.ageyoung))
            
            ncmaxf = (self.params.ncmaxfyoung - (self.params.ncmaxfyoung -
                        self.params.ncmaxfold) * age_effect)
                   
            if float_lt(ncmaxf, self.params.ncmaxfold):
                ncmaxf = self.params.ncmaxfold
    
            if float_gt(ncmaxf, self.params.ncmaxfyoung):
                ncmaxf = self.params.ncmaxfyoung
    
            # if foliage or root n:c ratio exceeds its max, then nitrogen 
            # uptake is cut back n.b. new ring n/c max is already set because 
            # it is related to leaf n:c
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
                self.fluxes.nuptake -= extrar 
            
    def calculate_cn_store(self, tolerance=1.0E-05):        
        cgrowth = (self.fluxes.cpleaf + self.fluxes.cproot + 
                   self.fluxes.cpbranch + self.fluxes.cpstem)
        ngrowth = (self.fluxes.npleaf + self.fluxes.nproot + 
                   self.fluxes.npbranch + self.fluxes.npstemimm + 
                   self.fluxes.npstemmob)
        
        # C storage as TNC
        self.state.cstore += self.fluxes.npp - cgrowth
        self.state.nstore += self.fluxes.nuptake + self.fluxes.retrans - ngrowth

    def calc_microbial_resp(self, project_day):
        """ Based on LPJ-why
        
        References:
        * Wania (2009) Integrating peatlands and permafrost into a dynamic 
          global vegetation model: 1. Evaluation and sensitivity of physical 
          land surface processes. GBC, 23, GB3014.
        * Also part 2 and the geosci paper in 2010
        """
        tsoil = self.met_data['tsoil'][project_day]
        
        # Lloyd and Taylor, 1994
        if tsoil < -10.0:
            temp_resp = 0.0
        else:
            temp_resp = (exp(308.56 * ((1.0 / 56.02) - (1.0 / 
                        (tsoil + const.DEG_TO_KELVIN - 227.13)))))
        
        moist_resp = ((1.0 - exp(-1.0 * self.state.wtfac_tsoil)) /  
                        (1.0 - exp(-1.0)))
        
        # Pool turnover rate = every 2 weeks -> days. Check you have this right
        # and turn into a parameter if so
        k_exu10 = 0.0714128571 
        k_exu = k_exu10 * temp_resp * moist_resp
        
        return  self.state.exu_pool * (1.0 - exp(-k_exu))



        
if __name__ == "__main__":
    
    # timing...
    import sys
    import time
    start_time = time.time()
    
    from file_parser import initialise_model_data
    from utilities import float_lt, day_length
    import datetime

    fname = "/Users/mdekauwe/src/python/pygday/params/duke_testing.cfg"

    (control, params, state, files,
        fluxes, met_data,
            print_opts) = initialise_model_data(fname, DUMP=False)

    # figure out photosynthesis
    PG = PlantGrowth(control, params, state, fluxes, met_data)

    state.lai = (params.slainit * const.M2_AS_HA /
                            const.KG_AS_TONNES / params.cfracts *
                            state.shoot)

    # Specific LAI (m2 onesided/kg DW)
    state.sla = params.slainit


    year = str(control.startyear)
    month = str(control.startmonth)
    day = str(control.startday)
    datex = datetime.datetime.strptime((year + month + day), "%Y%m%d")

    #laifname = "/Users/mdekauwe/research/NCEAS_face/GDAY_duke_simulation/experiments/lai"
    #import numpy as np
    #laidata = np.loadtxt(laifname)

    fdecay = 0.5
    rdecay = 0.5
    fluxes.deadleaves = 0.0
    fluxes.ceaten = 0.0
    fluxes.neaten = 0.0
    fluxes.deadroots = 0.0
    fluxes.deadbranch = 0.0         
    fluxes.deadstems = 0.0
    for project_day in xrange(len(met_data['prjday'])):

        state.shootnc = state.shootn / state.shoot
        state.ncontent = (state.shootnc * params.cfracts /
                                state.sla * const.KG_AS_G)
        daylen = day_length(datex, params.latitude)
        state.wtfac_root = 1.0
        #state.lai = laidata[project_day]


        if float_lt(state.lai, params.lai_cover):
            frac_gcover = state.lai / params.lai_cover
        else:
            frac_gcover = 1.0

        state.light_interception = ((1.0 - exp(-params.kext *
                                            state.lai / frac_gcover)) *
                                            frac_gcover)



        PG.grow(project_day, datex, fdecay, rdecay)
        print fluxes.gpp / const.HA_AS_M2 * const.TONNES_AS_G



        datex += datetime.timedelta(days=1)
    end_time = time.time()
    sys.stderr.write("\nTotal simulation time: %.1f seconds\n\n" %
                                                    (end_time - start_time))
    
    