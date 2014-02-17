#!/usr/bin/env python
""" Model Any Terrestrial Ecosystem (MATE) model. Full description below """

from math import exp, sqrt, sin, pi
import constants as const
from utilities import float_eq, float_gt, float_lt
import sys

__author__  = "Martin De Kauwe"
__version__ = "1.0 (04.08.2011)"
__email__   = "mdekauwe@gmail.com"


class MateC3(object):
    """ Model Any Terrestrial Ecosystem (MATE) model (C3)

    Simulates C3 photosynthesis (GPP) based on Sands (1995), accounting for 
    diurnal variations in irradiance and temp (am [sunrise-noon], 
    pm[noon to sunset]) and the decline of irradiance with depth through the 
    canopy.  
    
    MATE is connected to G'DAY via LAI and leaf N content. Plant autotrophic  
    respiration is calculated via carbon-use efficiency (CUE=NPP/GPP). 

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
        """
        self.params = params
        self.fluxes = fluxes
        self.control = control
        self.state = state
        self.met_data = met_data
        self.am = 0 # morning index
        self.pm = 1 # afternoon index        
        
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
        (Tair_K, par, vpd, ca) = self.get_met_data(day)
        
        # calculate mate params & account for temperature dependencies
        gamma_star = self.calculate_co2_compensation_point(Tair_K)
        Km = self.calculate_michaelis_menten_parameter(Tair_K)
        N0 = self.calculate_top_of_canopy_n()
        (jmax, vcmax) = self.calculate_jmax_and_vcmax(Tair_K, N0)
        ci = [self.calculate_ci(vpd[k], ca) for k in am, pm]
        
        # quantum efficiency calculated for C3 plants
        alpha = self.calculate_quantum_efficiency(ci, gamma_star)
        
        # Rubisco carboxylation limited rate of photosynthesis
        ac = [self.assim(ci[k], gamma_star[k], a1=vcmax[k], a2=Km[k]) \
              for k in am, pm]
        
        # Light-limited rate of photosynthesis allowed by RuBP regeneration
        aj = [self.assim(ci[k], gamma_star[k], a1=jmax[k]/4.0, \
              a2=2.0*gamma_star[k]) for k in am, pm]
        
        # light-saturated photosynthesis rate at the top of the canopy (gross)
        asat = [min(aj[k], ac[k]) for k in am, pm]
        
        # Assumption that the integral is symmetric about noon, so we average
        # the LUE accounting for variability in temperature, but importantly
        # not PAR
        lue = [self.epsilon(asat[k], par, daylen, alpha[k]) for k in am, pm]
        
        # mol C mol-1 PAR - use average to simulate canopy photosynthesis
        lue_avg = sum(lue) / 2.0

        if float_eq(self.state.lai, 0.0):
            self.fluxes.apar = 0.0
        else:
            par_mol = par * const.UMOL_TO_MOL
            # absorbed photosynthetically active radiation
            self.fluxes.apar = par_mol * self.state.fipar

        # gC m-2 d-1
        self.fluxes.gpp_gCm2 = (self.fluxes.apar * lue_avg * 
                                const.MOL_C_TO_GRAMS_C)
        self.fluxes.gpp_am_pm[am] = ((self.fluxes.apar / 2.0) * lue[am] * 
                                      const.MOL_C_TO_GRAMS_C)
        self.fluxes.gpp_am_pm[pm] = ((self.fluxes.apar / 2.0) * lue[pm] * 
                                      const.MOL_C_TO_GRAMS_C)
        
        self.fluxes.npp_gCm2 = self.fluxes.gpp_gCm2 * self.params.cue
        
        if self.control.nuptake_model == 3:
            self.fluxes.gpp_gCm2 *= self.params.ac
            self.fluxes.gpp_am_pm[am] *= self.params.ac
            self.fluxes.gpp_am_pm[pm] *= self.params.ac
            self.fluxes.npp_gCm2 = self.fluxes.gpp_gCm2 * self.params.cue
            
        # g C m-2 to tonnes hectare-1 day-1
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
        Tair_K : float
            am/pm air temperature in a list [Kelvin]
        vpd : float
            am/pm vpd in a list [kPa]
        par : float
            average daytime PAR [umol m-2 d-1]
        ca : float
            atmospheric co2 [umol mol-1]

        """
        Tair_K = [self.met_data['tam'][day] + const.DEG_TO_KELVIN, \
                self.met_data['tpm'][day] + const.DEG_TO_KELVIN]
        vpd = [self.met_data['vpd_am'][day], self.met_data['vpd_pm'][day]]
        ca = self.met_data["co2"][day]
        
        # if PAR is supplied by the user then use this data, otherwise use the
        # standard conversion factor. This was added as the NCEAS run had a
        # different scaling factor btw PAR and SW_RAD, i.e. not 2.3!
        if 'par' in self.met_data:
            par = self.met_data['par'][day] # umol m-2 d-1
        else:
            # convert MJ m-2 d-1 to -> umol m-2 day-1
            conv = const.RAD_TO_PAR * const.MJ_TO_MOL * const.MOL_TO_UMOL
            par = self.met_data['sw_rad'][day] * conv
        
        return (Tair_K, par, vpd, ca)

    def calculate_co2_compensation_point(self, Tk):
        """ CO2 compensation point in the absence of mitochondrial respiration
        Rate of photosynthesis matches the rate of respiration and the net CO2
        assimilation is zero.
        
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
        gamstar25 = self.params.gamstar25
        Egamma = self.params.Egamma
        
        return [self.arrh(gamstar25, Egamma, Tk[k]) for k in am, pm]
    
    def calculate_quantum_efficiency(self, ci, gamma_star):
        """ Quantum efficiency for AM/PM periods replacing Sands 1996 
        temperature dependancy function with eqn. from Medlyn, 2000 which is 
        based on McMurtrie and Wang 1993.
        
        References:
        -----------
        * Medlyn et al. (2000) Can. J. For. Res, 30, 873-888
        * McMurtrie and Wang (1993) PCE, 16, 1-13.
        
        """
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        
        return [self.assim(ci[k], gamma_star[k], a1=self.params.alpha_j/4.0, \
                a2=2.0*gamma_star[k]) for k in am, pm]
        
    
    def calculate_michaelis_menten_parameter(self, Tk):
        """ Effective Michaelis-Menten coefficent of Rubisco activity

        Parameters:
        ----------
        temp : float
            air temperature
        
        Returns:
        -------
        value : float, list [am, pm]
            Km, effective Michaelis-Menten constant for Rubisco catalytic 
            activity
        
        References:
        -----------
        Rubisco kinetic parameter values are from:
        * Bernacchi et al. (2001) PCE, 24, 253-259.
        * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.
        
        """
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        Ec = self.params.Ec
        Eo = self.params.Eo
        Kc25 = self.params.Kc25 
        Ko25 = self.params.Ko25 
        Oi = self.params.Oi 
        
        # Michaelis-Menten coefficents for carboxylation by Rubisco
        Kc = [self.arrh(Kc25, Ec, Tk[k]) for k in am, pm]
        
        # Michaelis-Menten coefficents for oxygenation by Rubisco
        Ko = [self.arrh(Ko25, Eo, Tk[k]) for k in am, pm]
        
        # return effectinve Michaelis-Menten coeffeicent for CO2
        return [Kc[k] * (1.0 + Oi / Ko[k]) for k in am, pm]
                
    def calculate_top_of_canopy_n(self):  
        """ Calculate the canopy N at the top of the canopy (g N m-2), N0.
        See notes and Chen et al 93, Oecologia, 93,63-69. 
        """
        
        if float_gt(self.state.lai, 0.0):
            # calculation for canopy N content at the top of the canopy                   
            N0 = (self.state.ncontent * self.params.kext /
                 (1.0 - exp(-self.params.kext * self.state.lai)))
        else:
            N0 = 0.0
        return N0
    
    
    def calculate_jmax_and_vcmax(self, Tk, N0):
        """ Calculate the maximum RuBP regeneration rate for light-saturated 
        leaves at the top of the canopy (Jmax) and the maximum rate of 
        rubisco-mediated carboxylation at the top of the canopy (Vcmax). 
        
        Parameters:
        ----------
        temp : float
            air temperature
        N0 : float
            leaf N
        """
    
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        deltaSj = self.params.delsj
        Eaj = self.params.eaj
        Eav = self.params.eav
        Hdj = self.params.edj
        
        if self.control.modeljm == True: 
            # the maximum rate of electron transport at 25 degC 
            jmax25 = self.params.jmaxna * N0 + self.params.jmaxnb
        
            jmax = [self.peaked_arrh(jmax25, Eaj, Tk[k], deltaSj, Hdj) \
                                    for k in am, pm]
            
            # the maximum rate of electron transport at 25 degC 
            vcmax25 = self.params.vcmaxna * N0 + self.params.vcmaxnb
            vcmax = [self.arrh(vcmax25, Eav, Tk[k]) for k in am, pm]
        else:
            jmax = [self.params.jmax, self.params.jmax]
            vcmax = [self.params.vcmax, self.params.vcmax]
        
        # reduce photosynthetic capacity with moisture stress
        jmax = [self.state.wtfac_root * jmax[k] for k in am, pm]
        vcmax = [self.state.wtfac_root * vcmax[k] for k in am, pm]  
    
        return jmax, vcmax
        
    def assim(self, ci, gamma_star, a1, a2):
        """Morning and afternoon calcultion of photosynthesis with the 
        limitation defined by the variables passed as a1 and a2, i.e. if we 
        are calculating vcmax or jmax limited.
        
        Parameters:
        ----------
        ci : float
            intercellular CO2 concentration.
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration
        a1 : float
            variable depends on whether the calculation is light or rubisco 
            limited.
        a2 : float
            variable depends on whether the calculation is light or rubisco 
            limited.

        Returns:
        -------
        assimilation_rate : float
            assimilation rate assuming either light or rubisco limitation.
        """
        if float_lt(ci, gamma_star):
            return 0.0
        else:
            return a1 * (ci - gamma_star) / (a2 + ci) 
   
    def calculate_ci(self, vpd, ca):
        """ Calculate the intercellular (Ci) concentration 

        Formed by substituting gs = g0 + 1.6 * (1 + (g1/sqrt(D))) * A/Ca into
        A = gs / 1.6 * (Ca - Ci) and assuming intercept (g0) = 0.

        Parameters:
        ----------
        vpd : float
            vapour pressure deficit
        ca : float
            ambient co2 concentration 
            
        Returns:
        -------
        ci:ca : float
            ratio of intercellular to atmospheric CO2 concentration

        References:
        -----------
        * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
        """
        g1w = self.params.g1 * self.state.wtfac_root
        cica = g1w / (g1w + sqrt(vpd))
        ci = cica * ca
        
        return ci
       
    def epsilon(self, asat, par, daylen, alpha):
        """ Canopy scale LUE using method from Sands 1995, 1996. 
        
        Sands derived daily canopy LUE from Asat by modelling the light response
        of photosysnthesis as a non-rectangular hyperbola with a curvature 
        (theta) and a quantum efficiency (alpha). 
        
        Assumptions of the approach are:
         - horizontally uniform canopy
         - PAR varies sinusoidally during daylight hours
         - extinction coefficient is constant all day
         - Asat and incident radiation decline through the canopy following 
           Beer's Law.
         - leaf transmission is assumed to be zero.
           
        * Numerical integration of "g" is simplified to 6 intervals. 

        Parameters:
        ----------
        asat : float
            photosynthetic rate at the top of the canopy
        par : float
            incident photosyntetically active radiation
        daylen : float
            length of day (hrs).
        theta : float
            curvature of photosynthetic light response curve 
        alpha : float
            quantum yield of photosynthesis (mol mol-1)
            
        Returns:
        -------
        lue : float
            integrated light use efficiency over the canopy (mol C mol-1 PAR)

        References:
        -----------
        See assumptions above...
        * Sands, P. J. (1995) Australian Journal of Plant Physiology, 
          22, 601-14.

        """
        delta = 0.16666666667 # subintervals scaler, i.e. 6 intervals
        h = daylen * const.HRS_TO_SECS 
        theta = self.params.theta # local var
        
        if float_gt(asat, 0.0):
            q = pi * self.params.kext * alpha * par / (2.0 * h * asat)
            integral_g = 0.0 
            for i in xrange(1, 13, 2):
                sinx = sin(pi * i / 24.)
                arg1 = sinx
                arg2 = 1.0 + q * sinx 
                arg3 = sqrt((1.0 + q * sinx)**2.0 - 4.0 * theta * q * sinx)
                integral_g += arg1 / (arg2 + arg3) * delta
            lue = alpha * integral_g * pi
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
            activation energy for the parameter [J mol-1]
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
        mt = self.params.measurement_temp + const.DEG_TO_KELVIN
        return k25 * exp((Ea * (Tk - mt)) / (mt * const.RGAS * Tk))
    
    def peaked_arrh(self, k25, Ea, Tk, deltaS, Hd):
        """ Temperature dependancy approximated by peaked Arrhenius eqn, 
        accounting for the rate of inhibition at higher temperatures. 

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]
        deltaS : float
            entropy factor [J mol-1 K-1)
        Hd : float
            describes rate of decrease about the optimum temp [J mol-1]
        
        Returns:
        -------
        kt : float
            temperature dependence on parameter 
        
        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179. 
        
        """
        mt = self.params.measurement_temp + const.DEG_TO_KELVIN
    
        arg1 = self.arrh(k25, Ea, Tk)
        arg2 = 1.0 + exp((mt * deltaS - Hd) / (mt * const.RGAS))
        arg3 = 1.0 + exp((Tk * deltaS - Hd) / (Tk * const.RGAS))
        
        return arg1 * arg2 / arg3



class MateC4(MateC3):
    """ Model Any Terrestrial Ecosystem (MATE) model (C4)

    Simulates C4 photosynthesis (GPP) based on Collatz (92) & Sands (1995), 
    accounting for diurnal variations in irradiance and temp (am [sunrise-noon], 
    pm[noon to sunset]) and the decline of irradiance with depth through the 
    canopy. The Collatz C4 model is a simplification, but functionally 
    equivalent model to the von Caemmerer model.
    
    MATE is connected to G'DAY via LAI and leaf N content. Plant autotrophic  
    respiration is calculated via carbon-use efficiency (CUE=NPP/GPP). 

    References:
    -----------
    * Collatz, G, J., Ribas-Carbo, M. and Berry, J. A. (1992) Coupled 
      Photosynthesis-Stomatal Conductance Model for Leaves of C4 plants. 
      Aust. J. Plant Physiol., 19, 519-38.
    * von Caemmerer, S. (2000) Biochemical Models of Leaf Photosynthesis. Chp 4. 
      Modelling C4 photosynthesis. CSIRO PUBLISHING, Australia. pg 91-122.

    Temperature dependancies:
    * Massad, R-S., Tuzet, A. and Bethenod, O. (2007) The effect of temperature 
      on C4-type leaf photosynthesis parameters. Plant, Cell and Environment, 
      30, 1191-1204.
      
    Intrinsic Quantum efficiency (mol mol-1), no Ci or temp dependancey 
    in c4 plants see: 
    * Ehleringer, J. R., 1978, Oecologia, 31, 255-267 or Collatz 1998.
    * Value taken from Table 1, Collatz et al.1998 Oecologia, 114, 441-454.
    """
    def __init__(self, control, params, state, fluxes, met_data):
        MateC3.__init__(self, control, params, state, fluxes, met_data)
        
        # curvature parameter, transition between light-limited and
        # carboxylation limited flux. Collatz table 2
        self.beta1 = 0.83	
        
        # curvature parameter, co-limitaiton between flux determined by
        # Rubisco and light and CO2 limited flux. Collatz table 2
        self.beta2 = 0.93
        
        # initial slope of photosynthetic CO2 response (mol m-2 s-1), 
        # Collatz table 2
        self.kslope = 0.7		    
        
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
        (Tair_K, par, vpd, ca) = self.get_met_data(day)

        ci = [self.calculate_ci(vpd[k], ca) for k in am, pm]
        N0 = self.calculate_top_of_canopy_n()
        alpha = self.params.alpha_c4
        
        # Temp dependancies from Massad et al. 2007
        (vcmax, vcmax25) = self.calculate_vcmax_parameter(Tair_K, N0)
        
        # Rubisco and light-limited capacity (Appendix, 2B)
        par_per_sec = par / (60.0 * 60.0 * daylen)
        M = [self.quadratic(a=self.beta1, b=-(vcmax[k] + alpha * par_per_sec), 
                            c=(vcmax[k] * alpha * par_per_sec)) for k in am, pm]

        # The limitation of the overall rate by M and CO2 limited flux:
        A = [self.quadratic(a=self.beta2, b=-(M[k] + self.kslope * ci[k]), 
                            c=(M[k] * self.kslope * ci[k])) for k in am, pm]

        # These respiration terms are just for assimilation calculations,
        # autotrophic respiration is stil assumed to be half of GPP
        (Rd) = self.calc_respiration(Tair_K, vcmax25)    

        # Net (saturated) photosynthetic rate, not sure if this
        # makes sense.
        Asat = [A[k] - Rd[k] for k in am, pm]
        
        # Assumption that the integral is symmetric about noon, so we average
        # the LUE accounting for variability in temperature, but importantly
        # not PAR
        lue = [self.epsilon(Asat[k], par, daylen, self.alphaf) for k in am, pm]

        # mol C mol-1 PAR - use average to simulate canopy photosynthesis
        lue_avg = sum(lue) / 2.0

        if float_eq(self.state.lai, 0.0):
            self.fluxes.apar = 0.0
        else:
            par_mol = par * const.UMOL_TO_MOL
            # absorbed photosynthetically active radiation
            self.fluxes.apar = par_mol * self.state.fipar
        
        # gC m-2 d-1
        self.fluxes.gpp_gCm2 = (self.fluxes.apar * lue_avg * 
                                const.MOL_C_TO_GRAMS_C)
        self.fluxes.gpp_am_pm[am] = ((self.fluxes.apar / 2.0) * lue[am] * 
                                      const.MOL_C_TO_GRAMS_C)
        self.fluxes.gpp_am_pm[pm] = ((self.fluxes.apar / 2.0) * lue[pm] * 
                                      const.MOL_C_TO_GRAMS_C)
        
        #print self.fluxes.gpp_gCm2
        self.fluxes.npp_gCm2 = self.fluxes.gpp_gCm2 * self.params.cue
        
        if self.control.nuptake_model == 3:
            self.fluxes.gpp_gCm2 *= self.params.ac
            self.fluxes.gpp_am_pm[am] *= self.params.ac
            self.fluxes.gpp_am_pm[pm] *= self.params.ac
            self.fluxes.npp_gCm2 = self.fluxes.gpp_gCm2 * self.params.cue
            
        # g C m-2 to tonnes hectare-1 day-1
        conv = const.G_AS_TONNES / const.M2_AS_HA
        self.fluxes.gpp = self.fluxes.gpp_gCm2 * conv
        self.fluxes.npp = self.fluxes.npp_gCm2 * conv
        
        # Plant respiration assuming carbon-use efficiency.
        self.fluxes.auto_resp = self.fluxes.gpp - self.fluxes.npp
        
    def calculate_vcmax_parameter(self, Tk, N0):
        """ Calculate the maximum rate of rubisco-mediated carboxylation at the
        top of the canopy
        
        # http://www.cesm.ucar.edu/models/cesm1.0/clm/CLM4_Tech_Note.pdf
        # Table 8.2 has PFT values...
        
        Parameters:
        ----------
        Tk : float
            air temperature (kelvin)
        N0 : float
            leaf N
            
        Returns:
        -------
        vcmax : float, list [am, pm]
            maximum rate of Rubisco activity
        """
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        #deltaS = self.params.delsv
        #Ea = self.params.eav
        #Hd = self.params.edv
        
        # Massad et al. 2007
        Ea = 67294.0
        Hd = 144568.0
        delS = 472.0
        
        # the maximum rate of electron transport at 25 degC 
        vcmax25 = self.params.vcmaxna * N0 + self.params.vcmaxnb
        vcmax = [self.peaked_arrh(vcmax25, Ea, Tk[k], delS, Hd) for k in am, pm]
        vcmax = [self.state.wtfac_root * vcmax[k] for k in am, pm] 
        
        return vcmax, vcmax25

    def calc_respiration(self, Tk, vcmax25, Tref=25.0):  
        """
        Mitochondrial respiration may occur in the mesophyll as well as in the 
        bundle sheath. As rubisco may more readily refix CO2 released in the 
        bundle sheath, Rd is described by its mesophyll and bundle-sheath 
        components: Rd = Rm + Rs
        
        Parameters:
        ----------
        Tk : float
            air temperature (kelvin)
        vcmax25 : float, list
            
        Returns:
        -------
        Rd : float, list [am, pm]
            (respiration in the light) 'day' respiration (umol m-2 s-1)
        
        
        References:
        -----------
        Tjoelker et al (2001) GCB, 7, 223-230.
        """
        am, pm = self.am, self.pm # morning/afternoon
        
        # scaling constant to Vcmax25, value = 0.015 after Collatz et al 1991. 
        # Agricultural and Forest Meteorology, 54, 107-136. But this if for C3
        # Value used in JULES for C4 is 0.025, using that one, see Clark et al.
        # 2011, Geosci. Model Dev, 4, 701-722.
        fdr = 0.025
        
        # specific respiration at a reference temperature (25 deg C)
        Rd25 = fdr * vcmax25
        
        # ratio between respiration rate at one temperature and the respiration 
        # rate at a temperature 10 deg C lower
        Q10 = 2.0 
        
        Rd = [Rd25 * Q10**(((Tk[k] - const.DEG_TO_KELVIN) - Tref) / 10.0) \
              for k in am, pm] 
       
        return Rd
        
    def quadratic(self, a=None, b=None, c=None):
        """ minimilist quadratic solution

        Parameters:
        ----------
        a : float
            co-efficient 
        b : float
            co-efficient
        c : float
            co-efficient

        Returns:
        -------
        root : float
    
        """
        d = b**2.0 - 4.0 * a * c # discriminant
        root = (-b - sqrt(d)) / (2.0 * a)	# Negative quadratic equation

        return root
       
         
if __name__ == "__main__":
    
    import numpy as np
    # timing...
    import sys
    import time
    start_time = time.time()
    
    from file_parser import initialise_model_data
    import datetime
    from utilities import float_eq, float_lt, float_gt, calculate_daylength, uniq
    
    met_header=4
    
    #fname = "/Users/mdekauwe/research/NCEAS_face/GDAY_ornl_simulation/params/NCEAS_or_youngforest.cfg"
    fname = "/Users/mdekauwe/research/FACE/GDAY_simulations/KSCO/experiment/params/NCEAS_KSCO_model_indust.cfg"
    (control, params, state, files,
        fluxes, met_data,
            print_opts) = initialise_model_data(fname, met_header, DUMP=False)
    
    control.ps_pathway = "C4"
    if control.ps_pathway == "C3":
        M = MateC3(control, params, state, fluxes, met_data)
    else:
        M = MateC4(control, params, state, fluxes, met_data)
   
    #flai = "/Users/mdekauwe/research/NCEAS_face/GDAY_ornl_simulation/experiments/silvias_LAI.txt"
    #lai_data = np.loadtxt(flai)
   

    # Specific LAI (m2 onesided/kg DW)
    state.sla = params.slainit

    
    project_day = 0
        
    # figure out the number of years for simulation and the number of
    # days in each year
    years = uniq(met_data["year"])
    #years = years[:-1] # dump last year as missing "site" LAI
    
    days_in_year = [met_data["year"].count(yr) for yr in years]
    
    for i, yr in enumerate(years):
        daylen = calculate_daylength(days_in_year[i], params.latitude)
        for doy in xrange(days_in_year[i]):
            state.wtfac_root = 1.0
            #state.lai = lai_data[project_day]
            state.lai = 2.0
            if state.lai > 0.0:
                state.shootnc = 0.03 #state.shootn / state.shoot
                state.ncontent = (state.shootnc * params.cfracts /
                                        state.sla * const.KG_AS_G)
            else:
                state.ncontent = 0.0        
        
        
            if float_lt(state.lai, params.lai_cover):
                frac_gcover = state.lai / params.lai_cover
            else:
                frac_gcover = 1.0

            state.fipar = ((1.0 - exp(-params.kext *
                                      state.lai / frac_gcover)) *
                                                frac_gcover)
            
            
            M.calculate_photosynthesis(project_day, daylen[doy])

            print fluxes.gpp_gCm2#, state.lai
        
            project_day += 1
