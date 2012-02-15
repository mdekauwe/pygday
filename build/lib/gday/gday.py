#!/usr/bin/env python
""" G'DAY is a process based model, which runs on a daily timestep and
simulates carbon, nutrient and water state and fluxes. See below for model
description.
"""

#import ipdb
import sys
import datetime
import calendar

import constants as const
from file_parser import initialise_model_data
from plant_growth import PlantGrowth
#from plant_growth_maestra import PlantGrowth
from print_outputs import PrintOutput
from litter_production import LitterProduction
from soil_cnflows import CarbonFlows, NitrogenFlows
from update_pools import CarbonPools, NitrogenPools
from utilities import float_eq, float_lt

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

    If site specific intial model pool information is not avaliable
    (e.g from the literature), then it is recommended that you "spin up" the
    model so that the initial pools reach a steady state.

    Model expects as input (on a daily timestep): (i) PAR, (ii) mean temp,
    (iii) rainfall and (iv) VPD. There is an external wrapper script to
    generate an input file from tmin, tmax and par.

    GPP ----->  ------------ --> Transpiration
                |  Foliage |
                ------------
                ------------
                |   Wood   |--------------------------->-----------------------
                ------------                                                  |
                ------------                                                  |
                |Fine Roots| <--<--                                           |
    CO2 <-----  ------------    | |                                           |
                                | |                                           |
                                | |                                           |
                                | |                                           |
    Nin  --> -------------- ----  |                                           |
             | Mineral N  | <-----|---  -----------------------------         |
    Nloss <- --------------       |     | Litter + Active Soil Pool | <-------
                                  | --> ----------------------------- ----> CO2
                                  | |   -----------------------------
                                  | --> |      Slow Soil Pool       | ----> CO2
    Rain -->  ---------- ----->---  -<- -----------------------------
              | Soil   |            |   -----------------------------
              | Water  |            --->|       Passive Soil Pool   | ----> CO2
    Drain <-- ---------- --> Evap       -----------------------------

    * 4 litter pools + active SOM pool are drawn together here but treated
      seperately in the model.


    References:
    ----------
    * Comins, H. N. and McMurtrie, R. E. (1993) Ecological Applications, 3,
      666-681.
    * Medlyn, B. E. et al (2000) Canadian Journal of Forest Research, 30,
      873-888.

    """
    def __init__(self, fname=None, chk_cmd_line=True, DUMP=False):

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

        # sweep the cmd line
        if chk_cmd_line == True:
            options, args = cmdline_parser()
            DUMP=options.DUMP

        (self.control, self.params,
            self.state, self.files,
            self.fluxes, self.met_data,
            print_opts) = initialise_model_data(fname, DUMP=DUMP)

        # printing stuff
        self.pr = PrintOutput(self.params, self.state, self.fluxes,
                                self.control, self.files, print_opts)

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

        # start date of simulation
        self.date = self.simulation_start_date()

        self.time_constants = ['rateuptake', 'rateloss', 'retransmob',
                                'fdecay', 'fdecaydry', 'rdecay', 'rdecaydry',
                                'bdecay', 'wdecay', 'kdec1', 'kdec2', 'kdec3',
                                'kdec4', 'kdec5', 'kdec6', 'kdec7', 'nuptakez']
        self.correct_rate_constants(output=False)
        
    def run_sim(self):
        """ Run model simulation! """

        # class instances
        cf = CarbonFlows(self.control, self.params, self.state, self.fluxes,
                            self.met_data)
        nf = NitrogenFlows(self.control, self.params, self.state, self.fluxes)
        lf = LitterProduction(self.control, self.params, self.state,
                                self.fluxes)
        pg = PlantGrowth(self.control, self.params, self.state, self.fluxes,
                            self.met_data)
        cpl = CarbonPools(self.control, self.params, self.state, self.fluxes)
        npl = NitrogenPools(self.control, self.params, self.state, self.fluxes,
                            self.met_data)
        
        # calculate initial C:N ratios and zero annual flux sums
        self.derive(1, self.date, INIT=True)
        
        for project_day in xrange(len(self.met_data['prjday'])):

            # litterfall rate: C and N fluxes
            (fdecay, rdecay) = lf.calculate_litter_flows()

            # co2 assimilation, N uptake and loss
            pg.grow(project_day, self.date, fdecay, rdecay)

            # soil model fluxes
            cf.calculate_cflows(project_day)
            nf.calculate_nflows()

            self.fluxes.nep = self.calculate_nep()
            
            # soil model - update pools
            (cact, cslo, cpas) = cpl.calculate_cpools()
            npl.calculate_npools(cact, cslo, cpas, project_day)

            # calculate C:N ratios and increment annual flux sums
            self.derive(project_day, self.date)
            
            #print self.state.plantc, self.state.soilc
            
            if self.control.print_options == 0:
                self.pr.save_daily_output(project_day + 1, self.date)

            self.increment_date()
            
            #print self.date.year, self.date.day, self.fluxes.rnet
            #print self.fluxes.gpp_gCm2# , self.fluxes.transpiration
            #print self.fluxes.transpiration
            
            #print self.fluxes.gpp_gCm2, self.met_data['amb_co2'][project_day]
            #print self.state.shootn, self.state.rootn, self.state.branchn, self.state.stemnimm, self.state.stemnmob
            #print self.fluxes.gpp_gCm2, self.fluxes.gs
            
            
            #print self.fluxes.gpp_gCm2
            
            #print self.date.year, self.date.month, self.state.pawater_root / self.params.wcapac_root * 100
            
            #print self.date.year, self.date.day, self.state.ncontent, self.state.lai
            #print self.fluxes.transpiration
            #print self.fluxes.gpp_gCm2
            #print self.params.g1 * self.state.wtfac_root, self.state.pawater_root / self.params.wcapac_root,self.fluxes.npp_gCm2 / self.fluxes.transpiration  
            #print self.state.stemn * 100.0
            
            #print self.params.g1 * self.state.wtfac_root, self.fluxes.npp_gCm2 / self.fluxes.transpiration
            
            #print self.state.soilc, self.state.plantc, self.fluxes.nep
                    
            
        if self.control.print_options == 1:
            # need to save initial SLA to current one!
            self.params.slainit = (self.state.lai / const.M2_AS_HA *
                                    const.KG_AS_TONNES * self.params.cfracts /
                                    self.state.shoot)
            self.correct_rate_constants(output=True)
            self.pr.save_state()

        # house cleaning, close ouput files
        self.pr.tidy_up()
       
    def simulation_start_date(self):
        """ figure out when the simulation starts

        Returns:
        -------
        start date : string (date format)
            Date string, year/month/day

        """
        year = str(self.control.startyear)
        month = str(self.control.startmonth)
        day = str(self.control.startday)
        return datetime.datetime.strptime((year + month + day), "%Y%m%d")

    def increment_date(self):
        """ move date object on a day, assumes daily date will need to adapt

        """
        # Total number of days in year
        if calendar.isleap(self.date.year):
            yr_days = 366.
        else:
            yr_days = 365.
        
        #Required so max leaf & root N:C can depend on Age 
        self.state.age += 1.0 / yr_days
        
        self.date += datetime.timedelta(days=1)

    def calculate_nep(self):
        """ carbon sink or source?

        Returns:
        --------
        NEP : float
            Net Ecosystem Productivity, C uptake
        """
        return (self.fluxes.npp - self.fluxes.hetero_resp -
                    self.fluxes.ceaten * (1. - self.params.fracfaeces))

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
    
    def derive(self, day, date, INIT=False):
        """Calculate derived values from state variables.

        Parameters:
        -----------
        day : integer
            day of simulation
        date : date format string
            date object yr/month/day
        INIT : logical
            logical defining whether it is the first day of the simulation

        """
        self.fluxes.ninflow = self.met_data['ndep'][day]
        
        # c/n ratios, most of these are just diagnostics, and not used.
        self.state.rootnc = nc_ratio(self.state.root, self.state.rootn, "Cr")
        self.state.shootnc = nc_ratio(self.state.shoot, self.state.shootn, "Cf")
        
        # Diagnostic N:C
        #branchnc = nc_ratio(self.state.branch, self.state.branchn)
        #stemnc = nc_ratio(self.state.stem, self.state.stemn)
        #structsurfnc = nc_ratio(self.state.structsurf, self.state.structsurfn)
        #metabsurfnc = nc_ratio(self.state.metabsurf, self.state.metabsurfn)
        #structsoilnc = nc_ratio(self.state.structsoil, self.state.structsoiln)
        #metabsoilnc = nc_ratio(self.state.metabsoil, self.state.metabsoiln)
        #activesoilnc = nc_ratio(self.state.activesoil, self.state.activesoiln)
        #slowsoilnc = nc_ratio(self.state.slowsoil, self.state.slowsoiln)
        #passivesoilnc = nc_ratio(self.state.passivesoil, self.state.passivesoiln)
        
        # SLA (m2 onesided/kg DW)
        self.state.sla = (self.state.lai / const.M2_AS_HA *
                            const.KG_AS_TONNES *
                            self.params.cfracts / self.state.shoot)

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
            # day of year 1-365/366
            doy = int(date.strftime('%j'))
            if doy == 1:
                self.state.nepsum = (self.fluxes.nep * const.TONNES_AS_G *
                                        const.M2_AS_HA)
                self.state.nppsum = (self.fluxes.npp * const.TONNES_AS_G *
                                        const.M2_AS_HA)
            else:
                self.state.nepsum += (self.fluxes.nep * const.TONNES_AS_G *
                                        const.M2_AS_HA)
                self.state.nppsum += (self.fluxes.npp * const.TONNES_AS_G *
                                        const.M2_AS_HA)

            # N Net mineralisation, i.e. excess of N outflows over inflows
            self.fluxes.nmineralisation = (self.fluxes.ninflow + 
                                            self.fluxes.ngross +
                                            self.fluxes.nrootexudate -
                                            self.fluxes.nimmob +
                                            self.fluxes.nlittrelease)

            # evaluate c input/output rates for mineral soil and soil+litter
            # Not used anyway so I have commented them out, diagnostics
            # mineral soil
            #cinsoil = sum(self.fluxes.cstruct) + sum(self.fluxes.cmetab)

            # litter + mineral soil
            #cinlitt = (self.fluxes.deadleaves + self.fluxes.deadroots +
            #            self.fluxes.deadbranch + self.fluxes.deadstems)

            # output from mineral soil
            #coutsoil = (self.fluxes.co2_to_air[4] + self.fluxes.co2_to_air[5] +
            #            self.fluxes.co2_to_air[6])

            # soil decomposition rate=flux/pool
            #soildecomp = coutsoil / self.state.soilc

def nc_ratio(carbon_val, nitrogen_val, pool):
    """Calculate nitrogen:carbon ratios

    Parameters:
    ----------
    carbon_val : float
        C value
    nitrogen_val: float
        N value#

    Returns:
    --------
    value : float
        N:C ratio
    """
    if float_lt(carbon_val, 0.0):
        # Note, previously the else branch was set to 1E6...presumably this 
        # was a hack to deal with the scenario for example where 
        # self.state.metabsurf and self.state.metabsurfn both start at zero.  
        # This was fine as this ratio isn't used in the code. Since I have 
        # commented out these diagnostics we shouldn't end up here unless there 
        # really is an error!!
        msg = "Dianostic for %s pool N:C has invalid values C:%s, N:%s" % \
                (pool, carbon_val, nitrogen_val)
        raise ValueError(msg)
    return nitrogen_val / carbon_val
      
def cmdline_parser():
    """ Parse the command line for user options

    Returns:
    --------
    options : object
        various cmd line options supplied by the user
    args : objects
        list of arguments

    """
    from optparse import OptionParser

    desc = """The G'DAY (Generic Decomposition And Yield) model)"""
    clp = OptionParser("Usage: %prog [options] filename", description = desc)
    clp.add_option("-d", "--dump", action="store_true", dest="DUMP",
                    default=False, help="Dump a default .INI file")
    options, args = clp.parse_args()
    return options, args


def main():
    """ run a test case of the gday model """

    # pylint: disable=C0103
    # pylint: disable=C0324

    # timing...
    import time
    start_time = time.time()


    fname = "/Users/mdekauwe/src/python/pygday/params/duke_testing.cfg"
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
