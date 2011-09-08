""" Carbon production module, call photosynthesis model """

import math

import constants as const
from utilities import float_eq, float_lt, float_gt, day_length
from bewdy import Bewdy
from plant_production_mcmurtrie import PlantProdModel
from water_balance import WaterBalance, WaterLimitedNPP
from mate import Mate

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

    def grow(self, day, date, fdecay, rdecay):
        """Evolve plant state, photosynthesis, distribute N and C"

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
            if self.control.assim_model != 7:
                self.wl.adjust_cproduction(self.control.water_model)

        # leaf N:C as a fraction of Ncmaxyoung, i.e. the max N:C ratio of
        # foliage in young stand
        nitfac = min(1.0, self.state.shootnc / self.params.ncmaxfyoung)

        # figure out allocation fractions for C
        (alleaf, alroot, albranch, alstem) = self.allocate_carbon(nitfac)

        # Distribute new C and N through the system
        (ncbnew, ncwimm, ncwnew) = self.calculate_ncwood_ratios(nitfac)
        self.nitrogen_distribution(ncbnew, ncwimm, ncwnew, fdecay, rdecay, 
                                    alleaf, alroot, albranch, alstem)
        self.carbon_distribution(alleaf, alroot, albranch, alstem, nitfac)
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
        return (ncbnew, ncwimm, ncwnew)

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
        (self.state.wtfac_tsoil, 
            self.state.wtfac_root) = self.wb.calculate_soil_water_fac()

        # Estimate photosynthesis using an empirical model
        if self.control.assim_model >=0 and self.control.assim_model <= 4:
            self.pp.calculate_photosynthesis(day)
        # Estimate photosynthesis using the mechanistic BEWDY model
        elif self.control.assim_model >=5 and self.control.assim_model <= 6:
            # calculate plant C uptake using bewdy
            self.bw.calculate_photosynthesis(frac_gcover, date, day, daylen)
        # Estimate photosynthesis using the mechanistic MATE model. Also need to
        # calculate a water availability scalar to determine Ci:Ca reln.
        elif self.control.assim_model ==7:
            self.mt.calculate_photosynthesis(day, daylen)
        else:
            raise AttributeError('Unknown assimilation model')

    def allocate_carbon(self, nitfac):
        """Carbon allocation fractions to move photosynthate through the plant.
        Allocations to foliage tends to decrease with stand age and wood stock
        increases. In stressed (soil/nutrient) regions fine root allocations
        increases.

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

        """

        alleaf = (self.params.callocf + nitfac *
                    (self.params.callocf_crit - self.params.callocf))
        alroot = (self.params.callocr + nitfac *
                    (self.params.callocr_crit - self.params.callocr))

        albranch = (self.params.callocb + nitfac *
                    (self.params.callocb_crit - self.params.callocb))

        alstem = 1.0 - alleaf - alroot - albranch

        return (alleaf, alroot, albranch, alstem)

    def nitrogen_distribution(self, ncbnew, ncwimm, ncwnew, fdecay, rdecay, 
                                alleaf, alroot, albranch, alstem):
        """ Nitrogen distribution - allocate available N through system.
        N is first allocated to the woody component, surplus N is then allocated
        to the shoot and roots with flexible ratios.

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
        alleaf : float
            allocation fraction for shoot
        alroot : float
            allocation fraction for fine roots
        albranch : float
            allocation fraction for branches
        alstem : float
            allocation fraction for stem

        """
        # N retranslocated proportion from dying plant tissue and stored within
        # the plant
        retrans = self.nitrogen_retrans(fdecay, rdecay)
        
        self.fluxes.nuptake = self.calculate_nuptake()
        
        # N lost from system through leaching and gaseous emissions
        self.fluxes.nloss = self.params.rateloss * self.state.inorgn
    
        # total nitrogen to allocate 
        ntot = self.fluxes.nuptake + retrans

        # allocate N to pools with fixed N:C ratios
        self.fluxes.npbranch = self.fluxes.npp * albranch * ncbnew

        # N flux into new ring (immobile component -> structrual components)
        self.fluxes.npstemimm = self.fluxes.npp * alstem * ncwimm

        # N flux into new ring (mobile component -> can be retrans for new
        # woody tissue)
        self.fluxes.npstemmob = self.fluxes.npp * alstem * (ncwnew - ncwimm)

        # If we have allocated more N than we have available - cut back N prodn
        arg = (self.fluxes.npstemimm + self.fluxes.npstemmob +
                self.fluxes.npbranch)

        if float_gt(arg, ntot) and not self.control.fixleafnc:
            self.fluxes.npp *= (ntot / (self.fluxes.npstemimm +
                                self.fluxes.npstemmob + self.fluxes.npbranch))
            self.fluxes.npbranch = self.fluxes.npp * albranch * ncbnew
            self.fluxes.npstemimm = self.fluxes.npp * alstem * ncwimm
            self.fluxes.npstemmob = self.fluxes.npp * alstem * (ncwnew - ncwimm)

        ntot -= (self.fluxes.npbranch + self.fluxes.npstemimm +
                    self.fluxes.npstemmob)

        # allocate remaining N to flexible-ratio pools
        self.fluxes.npleaf = (ntot * alleaf / 
                                (alleaf + alroot * self.params.ncrfac))
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
    
    def carbon_distribution(self, alleaf, alroot, albranch, alstem, nitfac):
        """ C distribution - allocate available C through system

        Parameters:
        -----------
        alleaf : float
            allocation fraction for shoot
        alroot : float
            allocation fraction for fine roots
        albranch : float
            allocation fraction for branches
        alstem : float
            allocation fraction for stem
        nitfac : float
            leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0)
        """
        self.fluxes.cpleaf = self.fluxes.npp * alleaf
        self.fluxes.cproot = self.fluxes.npp * alroot
        self.fluxes.cpbranch = self.fluxes.npp * albranch
        self.fluxes.cpstem = self.fluxes.npp * alstem

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

        state.light_interception = ((1.0 - math.exp(-params.kext *
                                            state.lai / frac_gcover)) *
                                            frac_gcover)



        PG.grow(project_day, datex, fdecay, rdecay)
        print fluxes.gpp / const.HA_AS_M2 * const.TONNES_AS_G



        datex += datetime.timedelta(days=1)
    end_time = time.time()
    sys.stderr.write("\nTotal simulation time: %.1f seconds\n\n" %
                                                    (end_time - start_time))
    
    