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

    Simulates photosynthesis (GPP) based on Sands (1995), accounting for diurnal
    variations in irradiance and temp (am [sunrise-noon], pm[noon to sunset]) 
    and the decline of irradiance with depth through the canopy.  
    
    MATE is connected to G'DAY via LAI and leaf N content. Key feedback through 
    soil N mineralisation and plant N uptake. Plant respiration is calculated 
    via carbon-use efficiency (CUE=NPP/GPP). There is a further water limitation
    constraint on productivity through the ci:ca ratio.

    References:
    -----------
    * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    * McMurtrie, R. E. et al. (2008) Functional Change Biology, 35, 521-34.
    * Sands, P. J. (1995) Australian Journal of Plant Physiology, 22, 601-14.

    Rubisco kinetic parameter values are from:
    * Bernacchi et al. (2001) PCE, 24, 253-259.
    * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.
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
        gamstar25 : float
            co2 compensation partial pressure in the absence of dark resp at 
            25 degC [umol mol-1]
        Oi : float
            intercellular concentration of O2 [umol mol-1]
        Kc25 : float
            Michaelis-Menten coefficents for carboxylation by Rubisco at 
            25degC [umol mol-1]
        Ko25: float
            Michaelis-Menten coefficents for oxygenation by Rubisco at 
            25degC [umol mol-1]. Note value in Bernacchie 2001 is in mmol!!
        Ec : float
            Activation energy for carboxylation [J mol-1]
        Eo : float
            Activation energy for oxygenation [J mol-1]
        Egamma : float
            Activation energy at CO2 compensation point [J mol-1]
        """
        self.params = params
        self.fluxes = fluxes
        self.control = control
        self.state = state
        self.met_data = met_data
        self.am = 0 # morning index
        self.pm = 1 # afternoon index
        self.gamstar25 = 42.75
        self.Oi = 205000.0
        self.Kc25 = 404.9 
        self.Ko25 = 278400.0 
        self.Ec = 79430.0
        self.Eo = 36380.0   # Note there is a typo in the R mate code here...  
        self.Egamma = 37830.0
        
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
        (8) Leaf temperature is the same as the air temperature.

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
        (am, pm) = self.am, self.pm # morning/afternoon
        (temp, par, vpd, ca) = self.get_met_data(day)
        Tk = [temp[k] + const.DEG_TO_KELVIN for k in am, pm]
        
        # calculate mate parameters, e.g. accounting for temp dependancy
        gamma_star = self.calculate_co2_compensation_point(Tk)
        km = self.calculate_michaelis_menten_parameter(Tk)
        N0 = self.calculate_leafn()
        jmax = self.calculate_jmax_parameter(Tk, N0)
        vcmax = self.calculate_vcmax_parameter(Tk, N0)
        alpha = self.calculate_quantum_efficiency(temp)
        
        # calculate ratio of intercellular to atmospheric CO2 concentration.
        # Also allows productivity to be water limited through stomatal opening.
        cica = [self.calculate_ci_ca_ratio(vpd[k]) for k in am, pm]
        ci = [i * ca for i in cica]
        
        # store value as needed in water balance calculation
        self.fluxes.cica_avg = sum(cica) / len(cica)
        
        # Rubisco-limited rate of photosynthesis
        ac = [self.aclim(ci[k], gamma_star[k], km[k], vcmax[k]) for k in am, pm]
        
        # Light-limited rate of photosynthesis allowed by RuBP regeneration
        aj = [self.ajlim(jmax[k], ci[k], gamma_star[k]) for k in am, pm]
        
        # Note that these are gross photosynthetic rates. Response to elevated
        # [CO2] is reduced if N declines, but increases as gs declines.
        asat = [min(aj[k], ac[k]) for k in am, pm]
        
        # GPP is assumed to be proportional to APAR, where the LUE defines the
        # slope of this relationship. LUE, calculation is performed for morning 
        # and afternnon periods.
        lue = [self.epsilon(asat[k], par, daylen, alpha[k]) for k in am, pm]
        
        # mol C mol-1 PAR - use average to simulate canopy photosynthesis
        lue_avg = sum(lue) / len(lue)
        
        if float_eq(self.state.lai, 0.0):
            self.fluxes.apar = 0.0
        else:
            par_mol = par * const.UMOL_TO_MOL
            self.fluxes.apar = par_mol * self.state.light_interception

        # gC m-2 d-1
        self.fluxes.gpp_gCm2 = (self.fluxes.apar * lue_avg * 
                                const.MOL_C_TO_GRAMS_C)
        self.fluxes.npp_gCm2 = self.fluxes.gpp_gCm2 * self.params.cue
        self.fluxes.gpp_am_gCm2 = ((self.fluxes.apar / 2.0) * lue[am] * 
                                    const.MOL_C_TO_GRAMS_C)
        self.fluxes.gpp_pm_gCm2 = ((self.fluxes.apar / 2.0) * lue[pm] * 
                                    const.MOL_C_TO_GRAMS_C)
        
        # tonnes hectare-1 day-1
        conv = const.G_AS_TONNES / const.M2_AS_HA
        self.fluxes.gpp = self.fluxes.gpp_gCm2 * conv
        self.fluxes.npp = self.fluxes.npp_gCm2 * conv
    
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
        temp = [self.met_data['tam'][day], self.met_data['tpm'][day]]
        vpd = [self.met_data['vpd_am'][day], self.met_data['vpd_pm'][day]]

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
            
        return (temp, par, vpd, ca)

    def calculate_co2_compensation_point(self, Tk):
        """ CO2 compensation point in the absence of mitochondrial respiration

        Parameters:
        ----------
        temp : float
            air temperature
        
        Returns:
        -------
        gamma_star : float, list [am, pm]
            CO2 compensation point in the abscence of mitochondrial respiration
        """
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        
        return [self.arrh(self.gamstar25, self.Egamma, Tk[k]) for k in am, pm]
    
    def calculate_quantum_efficiency(self, temp):
        """ Quantum efficiency for AM/PM periods following Sands 1996, it
        declines linearly with increasing temperature.

        Parameters:
        ----------
        temp : float
            air temperature
        alpha0 : float
            
        alpha1 : float
            
        
        Returns:
        -------
        alpha : float, list [am, pm]
            mol co2 mol-1 PAR
        """
        # BM - Ellsworth's table gives values 0.062 and 0.068 for upper canopy 
        # pine, so using 0.6
        alpha0 = 0.06  # quantum efficiency at 20 degC.
        alpha1 = 0.016 # characterises strength of the temp dependance of alpha
        
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        return [alpha0 * (1.0 - alpha1 * (temp[k]-20.0)) for k in am, pm]
    
    def calculate_michaelis_menten_parameter(self, Tk):
        """ Effective Michaelis-Menten coefficent of Rubisco activity

        Parameters:
        ----------
        temp : float
            air temperature
        
        Returns:
        -------
        value : float, list [am, pm]
            km, effective Michaelis-Menten constant for Rubisco catalytic activity
        
        References:
        -----------
        Rubisco kinetic parameter values are from:
        * Bernacchi et al. (2001) PCE, 24, 253-259.
        * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.
        
        """
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        
        # Michaelis-Menten coefficents for carboxylation by Rubisco
        Kc = [self.arrh(self.Kc25, self.Ec, Tk[k]) for k in am, pm]
        
        # Michaelis-Menten coefficents for oxygenation by Rubisco
        Ko = [self.arrh(self.Ko25, self.Eo, Tk[k]) for k in am, pm]
        
        # return effectinve Michaelis-Menten coeffeicent for CO2
        return [Kc[k] * (1.0 + self.Oi / Ko[k]) for k in am, pm]
                
    def calculate_leafn(self):  
        """ Assumption leaf N declines exponentially through the canopy. Input N
        is top of canopy (N0). See notes and Chen et al 93, Oecologia, 93,63-69. 
        """
        if self.state.ncontent > 0.0:
            # calculation for Leaf N content, top of the canopy (N0), [g m-2]                       
            N0 = (self.state.ncontent * self.params.kext /
                 (1.0 - math.exp(-self.params.kext * self.state.lai)))
        else:
            N0 = 0.0
        return N0
        
    def calculate_jmax_parameter(self, Tk, N0):
        """ Calculate the maximum RuBP regeneration rate for light-saturated 
        leaves at the top of the canopy (proportional to leaf-N content). 

        Parameters:
        ----------
        temp : float
            air temperature
        N0 : float
            leaf N

        Returns:
        -------
        jmax : float, list [am, pm]
            maximum rate of electron transport
        """
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        deltaS = self.params.delsj
        Ea = self.params.eaj
        Hd = self.params.edj
        
        # calculate the maximum rate of electron transport at 25 degC 
        if self.control.deciduous_model: 
            jmax25 = 40.462 * N0 + 13.691
        else:
            jmax25 = self.params.jmaxn * N0
        
        return [self.peaked_arrh(jmax25, Ea, Tk[k], deltaS, Hd) for k in am, pm]
        
    def calculate_vcmax_parameter(self, Tk, N0):
        """ Max rate of electron transport, Jmax. 

        Parameters:
        ----------
        temp : float
            air temperature
        N0 : float
            leaf N
            
        Returns:
        -------
        vcmax : float, list [am, pm]
            maximum rate of Rubisco activity
        """
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        
        # calculate the maximum rate of Rubisco activity at 25 degC 
        if self.control.deciduous_model: 
            vcmax25 = 20.497 * N0 + 8.403
        else:
            vcmax25 = self.params.vcmaxn * N0
        
        return [self.arrh(vcmax25, self.params.eav, Tk[k]) for k in am, pm]
    
    def aclim(self, ci, gamma_star, km, vcmax):
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
        vcmax : float
            maximum rate of Rubisco activity

        Returns:
        -------
        assimilation_rate : float
            assimilation rate when Rubisco activity is limiting

        """
        return max(0.0, (ci - gamma_star) * vcmax) / (ci + km)

    def ajlim(self, jmax, ci, gamma_star):
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

        Formed by substituting gs = g0 + 1.6 * (1 + (g1/sqrt(D))) * A/Ca into
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
        g1w = self.params.g1 * self.state.wtfac_root
        return g1w / (g1w + math.sqrt(vpd))

    def epsilon(self, amax, par, daylen, alpha):
        """ Canopy scale LUE using method from Sands 1995, 1996.

        Parameters:
        ----------
        amax : float
            photosynthetic rate at the top of the canopy
        par : float
            incident photosyntetically active radiation
        daylen : float
            length of day in hours.
        theta : float
            curvature of photosynthetic light response curve 
        alpha : float
            quantum yield of photosynthesis (mol mol-1)
            
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
            q = (math.pi * self.params.kext * alpha * par /
                    (2.0 * daylen * const.HRS_TO_SECS * amax))

            # check sands but shouldn't it be 2 * q * sin x on the top?
            f = (lambda x: x / (1.0 + q * x + math.sqrt((1.0 + q * x)**2.0 -
                            4.0 * self.params.theta * q * x)))
            g = [f(math.sin(math.pi * i / 24.)) for i in xrange(1, 13, 2)]
            
            #Trapezoidal rule - seems more accurate
            gg = 0.16666666667 * sum(g)

            lue = alpha * gg * math.pi
        else:
            lue = 0.0

        return lue
    
    def arrh(self, k25, Ea, Tk):
        """ Temperature dependence of kinetic parameters is described by an
        Arrhenius function

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC
        Ea : float
            activation energy for the parameter [kJ mol-1]
        Tk : float
            leaf temperature [deg K]

        Returns:
        -------
        kt : float
            temperature dependence on parameter 
        
        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.   
        """
        return k25 * math.exp((Ea * (Tk - 298.15)) / (298.15 * const.RGAS * Tk))
    
    def peaked_arrh(self, k25, Ea, Tk, deltaS, Hd):
        """ Temperature dependancy approximated by peaked Arrhenius eqn, 
        accounting for the rate of inhibition at higher temperatures. 

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC
        Ea : float
            activation energy for the parameter [kJ mol-1]
        Tk : float
            leaf temperature [deg K]
        deltaS : float
            entropy factor [J mol-1 K-1)
        Hd : float
            describes rate of decrease about the optimum temp [KJ mol-1]
        
        Returns:
        -------
        kt : float
            temperature dependence on parameter 
        
        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179. 
        
        """
        arg1 = self.arrh(k25, Ea, Tk)
        arg2 = 1.0 + math.exp((298.15 * deltaS - Hd) / 298.15 * const.RGAS)
        arg3 = 1.0 + math.exp((Tk * deltaS - Hd) / Tk * const.RGAS)
        
        return arg1 * arg2 / arg3

if __name__ == "__main__":
    
    import numpy as np
    # timing...
    import sys
    import time
    start_time = time.time()
    
    from file_parser import initialise_model_data
    from utilities import float_lt, day_length
    import datetime

    
    fname = "/Users/mdekauwe/research/NCEAS_face/GDAY_duke_simulation/params/NCEAS_dk_youngforest.cfg"
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

    control.co2_conc = 0
    
    npp_sum = np.zeros(0)
    
    #for project_day in xrange(365):
    for project_day in xrange(len(met_data['prjday'])):
        
        #state.shootn = 0.072422739989 
        #state.shoot = 6.54133760655 
        #state.lai = 6.01803059803
        #state.sla = 4.6
        #params.cfracts = 0.5
        #params.g1 = 4.8
        #params.jmaxn = 60.0
        #params.vcmaxn = 30.61
        #params.theta = 0.75
        
        
        
        
        
        #state.shootn = 0.071 # ornl val
        
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


        #daylen = 10.0
       
        M.calculate_photosynthesis(project_day, daylen)

        print fluxes.gpp_gCm2
        #print fluxes.gpp / state.shootn
        npp_sum = np.append(npp_sum, fluxes.gpp_gCm2*0.5) 



        datex += datetime.timedelta(days=1)
    
    #print npp_sum.sum() / (state.shootn *100)
    
    end_time = time.time()
    sys.stderr.write("\nTotal simulation time: %.1f seconds\n\n" %
                                                    (end_time - start_time))