#!/usr/bin/env python
""" Model Any Terrestrial Ecosystem (MATE) model. Full description below """

import math
import constants as const
from utilities import float_eq, float_gt

__author__  = "Martin De Kauwe"
__version__ = "1.0 (04.08.2011)"
__email__   = "mdekauwe@gmail.com"


class Mate(object):
    """ Model Any Terrestrial Ecosystem (MATE) model

    Simulates photosyntehsis (GPP) based on Farquahar + von Caemmerer, using the
    canopy estimate of LUE derived from the methods of Sands. MATE is connected
    to G'DAY via LAI and leaf N content. Key feedback through soil N
    mineralisation and plant N uptake. Plant respiration is calculated via
    carbon-use efficiency (CUE=NPP/GPP). There is a further water limitation
    constraint on productivity through the ci:ca ratio.

    References:
    -----------
    * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    * McMurtrie, R. E. et al. (2008) Functional Change Biology, 35, 521-34.
    * Sands, P. J. (1995) Australian Journal of Plant Physiology, 22, 601-14.

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

        self.am = 0
        self.pm = 1

    def calculate_photosynthesis(self, day, daylen):
        """ Photosynthesis is calculated assuming GPP is proportional to APAR,
        a commonly assumed reln (e.g. Potter 1993, Myneni 2002). The slope of
        relationship btw GPP and APAR, i.e. LUE is modelled using the
        photosynthesis eqns from Sands.

        Assumptions:
        ------------
        (1) photosynthetic light response is a non-rectangular hyperbolic func
            of photon-flux density with a light-saturatred photosynthetic rate
            (Amax), quantum yield (alpha) and curvature (theta).
        (2) the canopy is horizontally uniform.
        (3) PAR distribution within the canopy obeys Beer's law.
        (4) light-saturated photosynthetic rate declines with canopy depth in
            proportion to decline in PAR
        (5) alpha + theta do not vary within the canopy
        (6) dirunal variation of PAR is sinusoidal.
        (7) The model makes no assumption about N within the canopy, however
            this version assumes N declines exponentially through the cnaopy.

        Parameters:
        ----------
        day : int
            project day.
        daylen : float
            length of day in hours.

        Returns:
        -------
        Nothing
            Method calculates GPP, NPP and Ra.

        """

        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        temp, par, vpd, ca = self.get_met_data(day)

        #jmax25 = self.params.jmaxn * self.state.ncontent
        #vcmax25 = self.params.vcmaxn * self.state.ncontent


        # Assumption leaf N declines exponentially through the canopy. Input N
        # is top of canopy (N0)
        N0 = (self.state.ncontent * self.state.lai * self.params.kext /
                (1.0 - math.exp(-self.params.kext * self.state.lai)))

        # calculate the maximum rate of Rubisco activity at 25 degC and the
        # maximum rate of electron transport at 25 degC
        jmax25 = self.params.jmaxn * N0
        vcmax25 = self.params.vcmaxn * N0

        (gamma_star, km, jmax,
            vmax) = self.calculate_mate_params(temp, jmax25, vcmax25)

        # calculate ratio of intercellular to atmospheric CO2 concentration.
        # Also allows productivity to be water limited through stomatal opening.
        cica = [self.calculate_ci_ca_ratio(vpd[k]) for k in am, pm]
        ci = [i * ca for i in cica]

        # store value as needed in water balance calculation
        self.fluxes.cica_avg = sum(cica) / len(cica)
        
        # Generally under low PAR photosynthesis is limited by the rate of 
        # electron transport in the light reactions. Whereas at high PAR 
        # photosynthesis is limited by rubisco. So leaves growing in the shade 
        # would achieve no gain investing in rubisco, therefore have low Vcmax. 
        # Sunlit leaves have high Vcmax to maximise the rate of photosynthsis.
        # Further, if rubisco is low, there is no need to have extra
        # chlorophyll to trap light, therefore low Vcmax is twinned with low
        # Jmax.
        
        # Rubisco-limited rate of photosynthesis
        ac = [self.ac(ci[k], gamma_star[k], km[k], vmax[k]) for k in am, pm]
        
        # Light-limited rate of photosynthesis allowed by RuBP regeneration
        aj = [self.aj(jmax[k], ci[k], gamma_star[k]) for k in am, pm]
        
        # Note that these are gross photosynthetic rates.
        asat = [min(aj[k], ac[k]) for k in am, pm]
        
        # LUE, calculation is performed for morning and afternnon periods
        lue = [self.epsilon(asat[k], par, daylen) for k in am, pm]
        
        # mol C mol-1 PAR - use average to simulate canopy photosynthesis
        lue_avg = sum(lue) / len(lue)

        if float_eq(self.state.lai, 0.0):
            self.fluxes.apar = 0.0
        else:
            par_mol = par * const.UMOL_TO_MOL
            self.fluxes.apar = par_mol * self.state.light_interception

        # gC m-2 d-1
        self.fluxes.npp_gCm2 = (self.fluxes.apar * lue_avg * self.params.cue *
                                    const.MOLE_C_TO_GRAMS_C)
        self.fluxes.gpp_gCm2 = self.fluxes.npp_gCm2 / self.params.cue


        # tonnes hectare-1 day-1
        self.fluxes.npp = (self.fluxes.npp_gCm2 * const.G_AS_TONNES /
                            const.M2_AS_HA)
        self.fluxes.gpp = self.fluxes.npp / self.params.cue

        # Plant respiration assuming carbon-use efficiency.
        self.fluxes.auto_resp = self.fluxes.gpp - self.fluxes.npp


    def get_met_data(self, day):
        """ Grab the days met data out of the structure and return day values.

        Parameters:
        ----------
        day : int
            project day.

        Returns:
        -------
        temp : float
            am/pm temp in a list [degC]
        vpd : float
            am/pm vpd in a list [kPa]
        par : float
            average daytime PAR [umol m-2 d-1]
        ca : float
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated. [umol mol-1]

        """
        temp = [ self.met_data['tam'][day], self.met_data['tpm'][day] ]
        vpd = [ self.met_data['vpd_am'][day], self.met_data['vpd_pm'][day] ]

        # if PAR is supplied by the user then use this data, otherwise use the
        # standard conversion factor. This was added as the NCEAS run had a
        # different scaling factor btw PAR and SW_RAD, i.e. not 2.3!
        if 'par' in self.met_data:
            par = self.met_data['par'][day] # umol m-2 d-1
        else:
            # convert MJ m-2 d-1 to -> umol m-2 day-1
            conv = const.RAD_TO_PAR * const.MJ_TO_MOL * const.MOL_TO_UMOL
            par = self.met_data['sw_rad'][day] * conv

        if self.control.co2_conc == 0:
            ca = self.met_data['amb_co2'][day]
        elif self.control.co2_conc == 1:
            ca = self.met_data['ele_co2'][day]

        return temp, par, vpd, ca

    def calculate_mate_params(self, temp, jmax25, vcmax25):
        """ Method to calculate the parameters used in the MATE model

        Parameters:
        ----------
        temp : float
            air temperature
        jmax25 : float
            Jmax at 25DegC, function of leaf N
        vcmax25 : float
            Vcmax at 25DegC, function of leaf N

        Returns:
        -------
        gamma_star : float, list [am, pm]
            CO2 compensation point in the abscence of mitochondrial respiration
        km : float, list [am, pm]
            effective Michaelis-Menten constant for Rubisco catalytic activity
        jmax : float, list [am, pm]
            maximum rate of electron transport
        vmax : float, list [am, pm]
            maximum rate of Rubisco activity

        References:
        -----------
        Rubisco kinetic parameter values are from:
        * Medlyn et al. (2002) PCE, 25, 1167-1179.

        """
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon

        # co2 compensation point in the absence of mitochondrial respiration
        gamma_star = [self.arrh(42.75, 37830.0, temp[k]) for k in am, pm]


        # effective Michaelis-Menten coefficent of Rubisco activity
        km = [self.arrh(404.9, 79430.0, temp[k]) *
                (1.0 + 205000.0 / self.arrh(278400.0, 36380.0, temp[k]))
                for k in am, pm]
        

        # max rate of electron transport and rubisco activity
        jmax = [self.jmaxt(temp[k], jmax25) for k in am, pm]
        vmax = [self.arrh(vcmax25, self.params.eav, temp[k]) for k in am, pm]
        
        return gamma_star, km, jmax, vmax


    def ac(self, ci, gamma_star, km, vmax):
        """Morning and afternoon calcultion of photosynthesis when Rubisco
        activity is limiting, Ac.

        Parameters:
        ----------
        ci : float
            intercellular CO2 concentration.
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration
        km : float
            effective Michaelis-Menten constant for Rubisco catalytic activity
            for CO2
        vmax : float
            maximum rate of Rubisco activity

        Returns:
        -------
        assimilation_rate : float
            assimilation rate when Rubisco activity is limiting

        """
        return max(0.0, (ci - gamma_star) * vmax) / (ci + km)

    def aj(self, jmax, ci, gamma_star):
        """Morning and afternoon calcultion of photosynthesis when
        ribulose-1,5-bisphosphate (RuBP)-regeneration is limiting

        Parameters:
        ----------
        jmax : float
            maximum rate of electron transport
        ci : float
            intercellular CO2 concentration.
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration

        Returns:
        -------
        assimilation_rate : float
            photosynthesis when RuBP regeneration is limiting

        """
        return (jmax * (ci - gamma_star)) / (4.0 * (ci + 2.0 * gamma_star))

    def calculate_ci_ca_ratio(self, vpd):
        """ Calculate the ratio of intercellular to atmos CO2 conc

        Formed by substituting gs = g0 + (1 + (g1/sqrt(D))) * A/Ca into
        A = gs / 1.6 * (Ca - Ci) and assuming intercept (g0) = 0.

        Parameters:
        ----------
        vpd : float
            vapour pressure deficit

        Returns:
        -------
        ci:ca : float
            ratio of intercellular to atmospheric CO2 concentration

        References:
        -----------
        * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
        """
        return (1.0 - ((1.6 * math.sqrt(vpd)) /
                (self.params.g1 * self.state.wtfac_root + math.sqrt(vpd))))


    def epsilon(self, amax, par, daylen):
        """ Integrate LUE using method from Sands 1995, 1996.

        Parameters:
        ----------
        amax : float
            light-saturated rate of photosynthesis
        par : float
            incident photosyntetically active radiation
        daylen : float
            length of day in hours.

        Returns:
        -------
        lue : float
            integrated light use efficiency over the canopy

        References:
        -----------
        See assumptions above...
        * Sands, P. J. (1995) Australian Journal of Plant Physiology, 22, 601-14.

        """
        if float_gt(amax, 0.0):
            q = (math.pi * self.params.kext * self.params.alpha * par /
                    (2.0 * daylen * const.HRS_TO_SECS * amax))

            # check sands but shouldn't it be 2 * q * sin x on the top?
            f = (lambda x: x / (1.0 + q * x + math.sqrt((1.0 + q * x)**2.0 -
                            4.0 * self.params.theta * q * x)))
            g = [f(math.sin(math.pi * i / 24.)) for i in xrange(1, 13, 2)]

            #6-point Gaussian quadrature (integration)
            #gg = (0.08566225 * (g[0 + g[5]) + 0.1803808 * (g[1] + g[4]) +
            #        0.233957 * (g[2] + g[3]))

            #Trapezoidal rule - seems more accurate
            gg = 0.16666666667 * sum(g)

            lue = self.params.alpha * gg * math.pi

        else:
            lue = 0.0

        return lue

    def arrh(self, kc, ea, tair):
        """ Arrhenius function: describes the effect of temperature on
        enzyme activity

        Parameters:
        ----------
        kc : float
            pre-exponential factor
        ea : float
            activation energy
        tair : float
            air temperature [degC]

        Returns:
        -------
        k : float
            rate constant

        """
        return (kc * math.exp(ea * (tair - 25.0) / const.RGAS /
                (tair + const.ABSZERO) / (25.0 + const.ABSZERO)))

    def jmaxt(self, tair, jmax25):
        """ Calculates the temperature dependences of Jmax

        Parameters:
        ----------
        tair : float
            air temperature [degC]
        jmax25 : float
            Jmax at 25DegC, function of leaf N

        Returns:
        -------
        jmaxt : float
            temperature dependence on Jmax

        """
        tref = 25.0 + const.ABSZERO
        tk = tair + const.ABSZERO
        arg1 = self.arrh(jmax25, self.params.eaj, tair)
        arg2 = (1 + math.exp((self.params.delsj * tref - self.params.edj) /
                const.RGAS / tref))
        arg3 = (1 + math.exp((self.params.delsj * tk - self.params.edj) /
                const.RGAS / tk))
        return arg1 * arg2 / arg3



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

    M = Mate(control, params, state, fluxes, met_data)

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



        #print state.shootn
        M.calculate_photosynthesis(project_day, daylen)

        print fluxes.gpp_gCm2




        datex += datetime.timedelta(days=1)
    end_time = time.time()
    sys.stderr.write("\nTotal simulation time: %.1f seconds\n\n" %
                                                    (end_time - start_time))