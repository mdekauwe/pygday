#!/usr/bin/env python
""" G'DAY is a process based model, which runs on a daily timestep and
simulates carbon, nutrient and water state and fluxes. See below for model
description.
"""

#import ipdb
import sys
import math
import constants as const
from file_parser import initialise_model_data
from plant_growth import PlantGrowth
from print_outputs import PrintOutput
from litter_production import LitterProduction
from soil_cnflows import CarbonFlows, NitrogenFlows
from update_pools import CarbonPools, NitrogenPools
from utilities import float_eq, calculate_daylength, uniq
from phenology import Phenology


__author__  = "Martin De Kauwe"
__version__ = "1.0 (15.02.2011)"
__email__   = "mdekauwe@gmail.com"


class Gday(object):
    """ The G'DAY (Generic Decomposition And Yield) model.

    GDAY simulates C, N and water cycling between the plant and the soil. The
    model is structured into three plant pools (foliage, wood and fine roots),
    four litter pools (above/below metabolic and structural litter) and three
    soil organic matter (SOM) pools with varying turnover rates (active, slow
    and passive). An adapted implementation of the CENTURY model simulates soil
    carbon and nutrient dynamics. There is an additional simple soil water
    balance module which can be used to limit growth. Model pools can be thought
    of as buckets as they don't have dimensions. This is in essence a port
    of the C++ code, but I am changing things at will (be warned)!

    References:
    ----------
    * Comins, H. N. and McMurtrie, R. E. (1993) Ecological Applications, 3,
      666-681.
    * Medlyn, B. E. et al (2000) Canadian Journal of Forest Research, 30,
      873-888.

    """
    def __init__(self, fname=None, DUMP=False, spin_up=False):

        """ Set up model

        Read meterological forcing file and user config file and adjust the
        model parameters, control or initial state attributes that are used
        within the code.

        Parameters:
        ----------
        fname : string
            filename of model parameters, including path
        chk_cmd_line : logical
            parse the cmd line?
        DUMP : logical
            dump a the default parameters to a file

        Returns:
        -------
        Nothing
            Controlling class of the model, runs things.

        """
        self.day_output = [] # store daily outputs
        
        (self.control, self.params,
            self.state, self.files,
            self.fluxes, self.met_data,
            self.print_opts) = initialise_model_data(fname, DUMP=DUMP)

        # printing stuff
        self.pr = PrintOutput(self.params, self.state, self.fluxes,
                              self.control, self.files, self.print_opts)
        
        # print model defaults
        if DUMP == True:
            self.pr.save_default_parameters()
            sys.exit(0)

        if self.control.print_options > 1:
            raise ValueError("Unknown output print option: %s  (try 0 or 1)" %
                             self.control.print_options)

        # set initial lai -> m2/m2
        self.state.lai = (self.params.slainit * const.M2_AS_HA /
                          const.KG_AS_TONNES / self.params.cfracts *
                          self.state.shoot)
        
        # Specific leaf area (m2 onesided/kg DW)
        self.state.sla = self.params.slainit

        self.time_constants = ['rateuptake', 'rateloss', 'retransmob',
                                'fdecay', 'fdecaydry', 'rdecay', 'rdecaydry',
                                'bdecay', 'wdecay', 'kdec1', 'kdec2', 'kdec3',
                                'kdec4', 'kdec5', 'kdec6', 'kdec7', 'nuptakez']
        self.correct_rate_constants(output=False)
        
        # class instances
        self.cf = CarbonFlows(self.control, self.params, self.state, 
                              self.fluxes, self.met_data)
        self.nf = NitrogenFlows(self.control, self.params, self.state, 
                                self.fluxes)
        self.lf = LitterProduction(self.control, self.params, self.state,
                                   self.fluxes)
        self.pg = PlantGrowth(self.control, self.params, self.state, 
                              self.fluxes, self.met_data)
        self.cpl = CarbonPools(self.control, self.params, self.state, 
                               self.fluxes)
        self.npl = NitrogenPools(self.control, self.params, self.state, 
                                 self.fluxes, self.met_data)
        
        if self.control.deciduous_model:
            self.initialise_deciduous_model()
            self.P = Phenology(self.fluxes, self.state, 
                                self.params.previous_ncd)
        
        # calculate initial C:N ratios and zero annual flux sums
        self.day_end_calculations(0, INIT=True)
        self.state.pawater_root = self.params.wcapac_root
        self.state.pawater_tsoil = self.params.wcapac_topsoil
        self.spin_up = spin_up
        
    def spin_up_pools(self, tolerance=1E-03, sequence=1000):
        """ Spin Up model plant, soil and litter pools.
        -> Examine sequences of 1000 years and check if C pools are changing
           or at steady state to 3 d.p.
           
        References:
        ----------
        Adapted from...
        * Murty, D and McMurtrie, R. E. (2000) Ecological Modelling, 134, 
          185-205, specifically page 196.
        """
        prev_plantc = -9999.9
        prev_soilc = -9999.9
        prev_litterc = -9999.9
        while (math.fabs(prev_plantc - self.state.plantc) > tolerance and 
               math.fabs(prev_soilc - self.state.soilc) > tolerance and 
               math.fabs(prev_litterc - self.state.litterc) > tolerance): 
            prev_plantc = self.state.plantc
            prev_soilc = self.state.soilc
            prev_litterc = self.state.litterc
            self.run_sim() # run the model...
            
            # Have we reached a steady state?
            sys.stderr.write("Nyears of spin: %d %f %f %F\n" % \
                            (sequence, self.state.plantc, self.state.soilc, \
                             self.state.litterc))
            sequence += 1000
        self.print_output_file()    
            
    def run_sim(self):
        """ Run model simulation! """
        project_day = 0
        for yr in uniq(self.met_data["year"]):
            days_in_year = len([x for x in self.met_data["year"] if x == yr])
            daylen = calculate_daylength(days_in_year, self.params.latitude)
            
            if self.control.deciduous_model:
                self.zero_annual_sums()
                self.P.calculate_phenology_flows(daylen, self.met_data, 
                                            days_in_year, project_day)
            for doy in xrange(days_in_year):   
                
                # litterfall rate: C and N fluxes
                (fdecay, rdecay) = self.lf.calculate_litter_flows(doy)
                
                # co2 assimilation, N uptake and loss
                self.pg.grow(project_day, fdecay, rdecay, daylen[doy], doy, 
                        float(days_in_year))
    
                # soil model fluxes
                self.cf.calculate_cflows(project_day)
                self.nf.calculate_nflows()
                
                # soil model - update pools
                (cact, cslo, cpas) = self.cpl.calculate_cpools()
                self.npl.calculate_npools(cact, cslo, cpas, project_day)
    
                # calculate C:N ratios and increment annual flux sums
                self.day_end_calculations(project_day, days_in_year)
                
                #if self.spin_up == False:
                #    print self.fluxes.gpp * 100, self.state.lai
                
                # save daily fluxes + state for daily output    
                if self.control.print_options == 0:
                    self.save_daily_outputs(yr, doy+1)
                project_day += 1
            # =============== #
            #   END OF YEAR   #                
            # =============== #
            if self.control.deciduous_model:
                self.allocate_stored_c_and_n() 
        
        if self.spin_up == False:
            self.print_output_file()
        
    def print_output_file(self):
        """ End of the simulation...either print the daily output file or
        print the final state + param file. """
        
        # print the daily output file.
        if self.control.print_options == 0:
            self.pr.write_daily_outputs_file(self.day_output)
        # print the final state
        elif self.control.print_options == 1:
            if not self.control.deciduous_model:
                # need to save initial SLA to current one!
                conv = const.M2_AS_HA * const.KG_AS_TONNES 
                self.params.slainit = (self.state.lai / const.M2_AS_HA * 
                                      const.KG_AS_TONNES *
                                       self.params.cfracts /self.state.shoot)
                
            self.correct_rate_constants(output=True)
            self.pr.save_state()

       
    
    def zero_annual_sums(self):
        self.state.shoot = 0.0
        self.state.shootn = 0.0
        self.state.shootnc = 0.0
        self.state.clabile_store = 0.0
        self.state.aroot_uptake = 0.0
        self.state.aretrans = 0.0
        self.state.anloss = 0.0
        self.state.lai = 0.0
        
    def correct_rate_constants(self, output=False):
        """ adjust rate constants for the number of days in years """
        if output == False:
            for i in self.time_constants:
                setattr(self.params, i, getattr(self.params, i) / 
                        const.NDAYS_IN_YR) 
        else:
            for i in self.time_constants:
                setattr(self.params, i, getattr(self.params, i) * 
                        const.NDAYS_IN_YR)
    
    def day_end_calculations(self, prjday, days_in_year=None, INIT=False):
        """Calculate derived values from state variables.

        Parameters:
        -----------
        day : integer
            day of simulation
        
        INIT : logical
            logical defining whether it is the first day of the simulation

        """
        self.fluxes.ninflow = self.met_data['ndep'][prjday]
       
        # update N:C of plant pools
        if self.control.deciduous_model:
            if float_eq(self.state.shoot, 0.0):
                self.state.shootnc = 0.0
            else:
                self.state.shootnc = max(self.state.shootn / self.state.shoot, 
                                         self.params.ncfmin)
            self.state.rootnc = 0.02
        else:
            self.state.shootnc = self.state.shootn / self.state.shoot 
            self.state.rootnc = self.state.rootn / self.state.root
        
        if self.state.lai > 0.0:
            # SLA (m2 onesided/kg DW) -> HA/tonnes C
            self.state.sla = (self.state.lai / const.M2_AS_HA * 
                              const.KG_AS_TONNES * self.params.cfracts / 
                              self.state.shoot)

        # total plant, soil & litter nitrogen
        self.state.soiln = (self.state.inorgn + self.state.activesoiln +
                            self.state.slowsoiln + self.state.passivesoiln)
        self.state.litternag = self.state.structsurfn + self.state.metabsurfn
        self.state.litternbg = self.state.structsoiln + self.state.metabsoiln
        self.state.littern = self.state.litternag + self.state.litternbg
        self.state.plantn = (self.state.shootn + self.state.rootn +
                             self.state.branchn + self.state.stemn)
        self.state.totaln = (self.state.plantn + self.state.littern +
                             self.state.soiln)

        # total plant, soil, litter and system carbon
        self.state.soilc = (self.state.activesoil + self.state.slowsoil +
                            self.state.passivesoil)
        self.state.littercag = self.state.structsurf + self.state.metabsurf
        self.state.littercbg = self.state.structsoil + self.state.metabsoil
        self.state.litterc = self.state.littercag + self.state.littercbg
        self.state.plantc = (self.state.root + self.state.shoot +
                             self.state.stem + self.state.branch)
        self.state.totalc = (self.state.soilc + self.state.litterc +
                             self.state.plantc)

        # optional constant passive pool
        if self.control.passiveconst != 0:
            self.state.passivesoil = self.params.passivesoilz
            self.state.passivesoiln = self.params.passivesoilnz

        if INIT == False:
            #Required so max leaf & root N:C can depend on Age 
            self.state.age += 1.0 / days_in_year
            
            # N Net mineralisation, i.e. excess of N outflows over inflows
            self.fluxes.nmineralisation = (self.fluxes.ninflow + 
                                            self.fluxes.ngross +
                                            self.fluxes.nrootexudate -
                                            self.fluxes.nimmob +
                                            self.fluxes.nlittrelease)

    def initialise_deciduous_model(self):
        # Divide up NPP based on annual allocation fractions
        
        self.state.c_to_alloc_shoot = (self.state.alleaf * 
                                        self.state.clabile_store)
        self.state.c_to_alloc_root = (self.state.alroot * 
                                        self.state.clabile_store)
        self.state.c_to_alloc_branch = (self.state.albranch * 
                                        self.state.clabile_store)
        self.state.c_to_alloc_stem = (self.state.alstem * 
                                        self.state.clabile_store)
        #self.state.c_to_alloc_rootexudate = (self.state.alroot_exudate *    
        #                                        self.state.clabile_store)
        
        # annual available N for allocation to leaf
        self.state.n_to_alloc_shoot = (self.state.c_to_alloc_shoot * 
                                        self.state.shootnc_yr)

    def allocate_stored_c_and_n(self):
        """ 
        At the end of the year allocate everything for the coming year
        based on stores from the previous year avaliable N for allocation
        """
        self.state.c_to_alloc_shoot = (self.state.alleaf * 
                                        self.state.clabile_store)
        
        Un = self.state.aroot_uptake + self.state.aretrans
        
        self.state.c_to_alloc_stem = (self.params.callocw * 
                                      (self.state.clabile_store - 
                                       self.state.c_to_alloc_stem))
        self.state.c_to_alloc_root = (self.state.clabile_store - 
                                      self.state.c_to_alloc_stem - 
                                      self.state.c_to_alloc_shoot)
        self.state.n_to_alloc_root = (min(Un, self.state.c_to_alloc_root * 
                                              self.state.rootnc))
        
        # constant N:C of foliage during the growing season(kgN kg-1C)
        self.state.shootnc_yr = ((Un - self.state.n_to_alloc_root) / 
                                 (self.state.c_to_alloc_shoot)) 
        # if we want to put back a floating N:C then we need to have
        # self.state.c_to_alloc_shoot + self.state.c_to_alloc_stem * some factor
        
        # annual available N for allocation to leaf
        self.state.n_to_alloc_shoot = (self.state.c_to_alloc_shoot * 
                                        self.state.shootnc_yr)    



    def save_daily_outputs(self, year, doy):
        """ Save the daily fluxes + state in a big list.
        
        This should be a more efficient way to write the daily output in a 
        single step at the end of the simulation. 

        Parameters:
        -----------
        project_day : integer
            simulation day
        """
        output = [year, doy]
        for var in self.print_opts:
            try:
                if hasattr(self.state, var):
                    value = getattr(self.state, var)
                    output.append(value)
                else:
                    value = getattr(self.fluxes, var)
                    output.append(value)
            except AttributeError:
                err_msg = "Error accessing var to print: %s" % var
                raise AttributeError, err_msg
        self.day_output.append(output)
        

def main():
    """ run a test case of the gday model """

    # pylint: disable=C0103
    # pylint: disable=C0324

    # timing...
    import time
    start_time = time.time()


    
    #fname = "/Users/mdekauwe/research/NCEAS_face/GDAY_duke_simulation/params/NCEAS_dk_youngforest.cfg"
    fname = "test.cfg"
    G = Gday(fname)
    G.run_sim()
    
    end_time = time.time()
    sys.stderr.write("\nTotal simulation time: %.1f seconds\n\n" %
                                                    (end_time - start_time))


def profile_main():
    """ profile code """
    import cProfile, pstats
    prof = cProfile.Profile()
    prof = prof.runctx("main()", globals(), locals())
    print "<pre>"
    stats = pstats.Stats(prof)
    stats.sort_stats("cumulative")  # Or cumulative
    stats.print_stats(500)  # 80 = how many to print
    # The rest is optional.
    # stats.print_callees()
    # stats.print_callers()
    print "</pre>"


if __name__ == "__main__":

    main()
    #profile_main()
