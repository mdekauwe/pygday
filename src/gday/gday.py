#!/usr/bin/env python
""" G'DAY is a process based model, which runs on a daily timestep and
simulates carbon, nutrient and water state and fluxes. See below for model
description.
"""

#import ipdb
import sys
import datetime

import constants as const
from file_parser import initialise_model_data
from plant_growth import PlantGrowth
from print_outputs import PrintOutput
from litter_production import LitterProduction
from soil_cnflows import CarbonFlows, NitrogenFlows
from update_pools import CarbonPools, NitrogenPools
from derive import Derive
from utilities import float_eq

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
                                'fdecay','fdecaydry', 'rdecay', 'rdecaydry',
                                'bdecay', 'wdecay', 'kdec1', 'kdec2', 'kdec3',
                                'kdec4', 'kdec5', 'kdec6', 'kdec7', 'nuptakez']
        self.correct_rate_constants(output=False)

    def run_sim(self):
        """ Run model simulation! """

        # class instances
        cf = CarbonFlows(self.control, self.params, self.state, self.fluxes,
                            self.met_data)
        nf = NitrogenFlows(self.control, self.params, self.state, self.fluxes)
        de = Derive(self.control, self.params, self.state, self.fluxes,
                        self.met_data)
        lf = LitterProduction(self.control, self.params, self.state,
                                self.fluxes)
        pg = PlantGrowth(self.control, self.params, self.state, self.fluxes,
                            self.met_data)
        cpl = CarbonPools(self.control, self.params, self.state, self.fluxes)
        npl = NitrogenPools(self.control, self.params, self.state, self.fluxes,
                            self.met_data)

        # calculate initial C:N ratios and zero annual flux sums
        de.derive_vals_from_state(1, self.date, INIT=True)

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
            de.derive_vals_from_state(project_day, self.date)

            if self.control.print_options == 0:
                self.pr.save_daily_output(project_day + 1, self.date)

            self.increment_date()

            #print self.fluxes.transpiration
            print self.fluxes.gpp_gCm2
            #print self.state.stemn * 100.0

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
