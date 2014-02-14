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
    """ Model Any Terrestrial Ecosystem (MATE) model

    Simulates photosynthesis (GPP) based on Sands (1995), accounting for diurnal
    variations in irradiance and temp (am [sunrise-noon], pm[noon to sunset]) 
    and the decline of irradiance with depth through the canopy.  
    
    MATE is connected to G'DAY via LAI and leaf N content. Key feedback through 
    soil N mineralisation and plant N uptake. Plant respiration is calculated 
    via carbon-use efficiency (CUE=NPP/GPP). There is a further water limitation
    constraint on productivity through the Ci:Ca ratio.

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
        alpha_j = 0.26 # initial slope
        
        # McM '08
        #return [0.07 * (ca - gamma_star[k]) / (cia  + (2.0*gamma_star[k])) \
        #        for k in am, pm]
        return [self.assim(ci[k], gamma_star[k], a1=alpha_j/4.0, \
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
    """ C4 photosynthesis pathway. Sands LUE implementation adjusted to make use
    of C4 equations from von Caemmerer.
    
    References:
    ===========
    * von Caemmerer, S. (2000) Biochemical Models of Leaf Photosynthesis. Chp 4. 
      Modelling C4 photosynthesis. CSIRO PUBLISHING, Australia. pg 91-122.
    * Massad, R-S., Tuzet, A. and Bethenod, O. (2007) The effect of temperature 
      on C4-type leaf photosynthesis parameters. Plant, Cell and Environment, 
      30, 1191-1204.
    
    """
    
    def __init__(self, control, params, state, fluxes, met_data):
        MateC3.__init__(self, control, params, state, fluxes, met_data)
        
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
        
        # Quantum efficiency (umol mol-1), no Ci or temp dependancey in c4
        # plants
        # Ehleringer, J. R., 1978: Implications of quantum yield differences 
        # on the distributions of C3 and C4 grasses.  Oecologia, 31, 255-267. 
        alpha = 0.04
        
        self.control.collatz = True
        if self.control.collatz:
            # C4 assimilation following from Collatz et al. 1992
            
            #?? Exponential factor in the equation defining kt
            #alpharf = 0.067			# mol/mol
            kslope = 0.7		    # initial slope of photosynthetic CO2 response (mol m-2 s-1), Collatz table 2
            theta = 0.83			# curvature parameter, Collatz table 2
            beta  = 0.93			# curvature parameter, Collatz table 2
            
            
            # http://www.cesm.ucar.edu/models/cesm1.0/clm/CLM4_Tech_Note.pdf
            # Table 8.2 has PFT values...
            vcmax25 = self.params.vcmaxna * N0 + self.params.vcmaxnb          
            
            
            # Massad et al. 2007
            Eav = 67294.0
            Hdv = 144568.0
            deltaSv = 472.0
            vcmax = [self.peaked_arrh(vcmax25, Eav, Tair_K[k], deltaSv, Hdv) for k in am, pm]
            
            
            Je = vcmax
            Jc = k / vmax * vcmax * ci
            Ji = alpha * I
            
            
            # Rubisco and light limited capacity
            par_per_sec = par / (60.0 * 60.0 * daylen)
            M = [self.quadratic(theta, -(alpha*par_per_sec+vcmax[k]), alpha*par_per_sec*vcmax[k]) for k in am, pm]

            # M and CO2 limitation
            A = [self.quadratic(beta, -(M[k]+kslope*ci[k]), M[k]*kslope*ci[k]) for k in am, pm]

            # These respiration terms are just for assimilation calculations,
            # autotrophic respiration is stil assumed to be half of GPP
            (Rd, Rm) = self.calc_respiration(Tair_K, vcmax)    


            # Net (saturated) photosynthetic rate, not sure if this
            # makes sense.
            Asat = [A[k] - Rd[k] for k in am, pm]
        else:    
        
            # calculate mate parameters, e.g. accounting for temp dependancy
            (Km, Kc, Ko, Kp) = self.calculate_michaelis_menten_parameter(Tair_K)
            gamma_star = self.calculate_co2_compensation_point(Tair_K)
        
            
    
        
            # Currently i dont have any information on how these depednancies 
            # vary with N, nor do I have site parameters so going to use
            # values at 25 degrees from the literature, hardwired till this is
            # resolved.
            # Values from table 4.1, in von Caemmerer 2000, pg 100.
            vcmax25 = 60.0
            vpmax25 = 120.0 
            jmax25 = 400.0
    
            if self.control.modeljm == True: 
                jmax = self.calculate_jmax_parameter(Tair_K, jmax25)
                vcmax = self.calculate_vcmax_parameter(Tair_K, vcmax25)
                vpmax = self.calculate_vpmax_parameter(Tair_K, vpmax25)
            else:
                jmax = [self.params.jmax, self.params.jmax]
                vcmax = [self.params.vcmax, self.params.vcmax]
                vpmax = [self.params.vpmax, self.params.vpmax]
        
        
        
            # reduce photosynthetic capacity with moisture stress
            #jmax = [self.state.wtfac_root * jmax[k] for k in am, pm]
            #vcmax = [self.state.wtfac_root * vcmax[k] for k in am, pm] 
            #vpmax = [self.state.wtfac_root * vpmax[k] for k in am, pm] 
    
            # These respiration terms are just for assimilation calculations,
            # autotrophic respiration is stil assumed to be half of GPP
            (Rd, Rm) = self.calc_respiration(Tair_K, vcmax)    
        
            # Calculate assimilation rates.
            Ac = self.calc_enzyme_limited_assim(Km, Kc, Ko, Kp, ci, vpmax, vcmax, 
                                                Rd, Rm)
            par_per_sec = par / (60.0 * 60.0 * daylen)
            Aj = self.calc_light_limited_assim(par_per_sec, jmax, ci, Rd, Rm)
            Asat = [min(Aj[k], Ac[k]) for k in am, pm]
    
        # Assumption that the integral is symmetric about noon, so we average
        # the LUE accounting for variability in temperature, but importantly
        # not PAR
        lue = [self.epsilon(Asat[k], par, daylen, alpha) for k in am, pm]

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
        
        print self.fluxes.gpp_gCm2
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
        * Massad et al. (2007) PCE, 30, 1191-1204.
        
        
        """
        # local var for tidyness
        am, pm = self.am, self.pm # morning/afternoon
        Tc = [Tk[k] - const.DEG_TO_KELVIN for k in am, pm]
        
        Kc25 = self.params.Kc25 
        Ko25 = self.params.Ko25 
        #Kp25 = self.params.Kp25 
        #Oi = self.params.Oi / 1000.0  # mbar
        Oi = 210.0
        # Massad et al 2007, Table 1.
        Kc25 = 650.0 # ubar
        Ko25 = 450.0 # mbar
        Kp25 = 80.0  # ubar
        Q10 = 2.0
        
        # Michaelis-Menten coefficents for carboxylation by Rubisco
        Kc = [Kc25 * Q10**((Tc[k] - 25.0) / 10.0) for k in am, pm]
        
        # Michaelis-Menten coefficents for oxygenation by Rubisco
        Ko = [Ko25 * Q10**((Tc[k] - 25.0) / 10.0) for k in am, pm]
        
        # Michaelis-Menten coefficents for PEPcase for CO2
        Kp = [Kp25 * Q10**((Tc[k] - 25.0) / 10.0) for k in am, pm]
        
        # return effectinve Michaelis-Menten coeffeicent for CO2
        return [Kc[k] * (1.0 + Oi / Ko[k]) for k in am, pm], Kc, Ko, Kp
 
    def calculate_jmax_parameter(self, Tk, jmax25):
        """ Calculate the maximum RuBP regeneration rate for light-saturated 
        leaves at the top of the canopy. 

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
        
        Ea = 77900.0
        deltaS = 627.0
        Hd = 191929.0
        
        
        # the maximum rate of electron transport at 25 degC 
        #jmax25 = self.params.jmaxna * N0 + self.params.jmaxnb
        
        #print jmax25, Ea, deltaS, Hd
        return [self.peaked_arrh(jmax25, Ea, Tk[k], deltaS, Hd) for k in am, pm]
    
        
    def calculate_vcmax_parameter(self, Tk, vcmax25):
        """ Calculate the maximum rate of rubisco-mediated carboxylation at the
        top of the canopy

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
        #deltaS = self.params.delsv
        #Ea = self.params.eav
        #Hd = self.params.edv
        
        Ea = 67294.0
        deltaS = 472.0
        Hd = 144568.0
        
        
        # the maximum rate of electron transport at 25 degC 
        #vcmax25 = self.params.vcmaxna * N0 + self.params.vcmaxnb
        
        
        return [self.peaked_arrh(vcmax25, Ea, Tk[k], deltaS, Hd) for k in am, pm]
        
    def calculate_vpmax_parameter(self, Tk, vpmax25):
        """ Calculate the carboxylation by PEPC, initial slope of the C4 A/ci
        curve

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
        #deltaS = self.params.delsvp
        #Ea = self.params.eavp
        #Hd = self.params.edvp
        
        Ea = 70373.0
        deltaS = 376.0
        Hd = 117910.0
        
        # the maximum rate of electron transport at 25 degC 
        #vpmax25 = self.params.vpmaxna * N0 + self.params.vpmaxnb
        
        
        return [self.peaked_arrh(vpmax25, Ea, Tk[k], deltaS, Hd) for k in am, pm]

    def calc_respiration(self, temp, vcmax):  
        """
        Mitochondrial respiration may occur in the mesophyll as well as in the 
        bundle sheath. As rubisco may more readily refix CO2 released in the 
        bundle sheath, Rd is described by its mesophyll and bundle-sheath 
        components: Rd = Rm + Rs
        
        Parameters:
        ----------
        temp : float
            air temperature
        vcmax : float, list
            
        Returns:
        -------
        Rd : float, list [am, pm]
            Day leaf respiration (umol m-2 s-1)
        Rm : float, list [am, pm]
            Maintanence respiration (umol m-2 s-1)
        """
        am, pm = self.am, self.pm # morning/afternoon
        # Day leaf respiration (umol m-2 s-1)
        Rd = [0.01 * vcmax[k] for k in am, pm]
        
        # Mesophyll mitochondrial respiration  (umol m-2 s-1)
        Rm = [0.5 * Rd[k] for k in am, pm]
        
        return (Rd, Rm) 

    def calc_pep_carboxylation_rate(self, ci, vpmax, Kp):
        """
        Parameters:
        ----------
        Vpr : float
            rate of PEP regeneration (mu mol m-2 s-1)
        """
        am, pm = self.am, self.pm # morning/afternoon
        vpr = self.params.vpr 
        
        return [min(ci[k] * vpmax[k] / (ci[k] + Kp[k]), vpr) for k in am, pm] 
    
    def calc_enzyme_limited_assim(self, Km, Kc, Ko, Kp, ci, Vpmax, Vcmax, Rd, Rm):
        """
        Calculate the AM/PM enzyme-limited CO2 assimilation rate.
        
        Parameters:
        ----------
        
        alpha : float
            Fraction of PSII activity in the bundle sheath [0-1]
        rub_sf : float 
            Half the reciprocal for Rubisco specificity 
        gbs : float
            bundle sheath conductance
        
        
        Returns:
        -------
        Ac : float, list [am, pm]
            enzyme-limited assimilation rate [umol m-2 s-1]
        """
    
        # local parameters
        am, pm = self.am, self.pm # morning/afternoon
        alpha = self.params.alpha_psii		
        #Oi = self.params.Oi / 1000.0 # mbar
        Oi = 210.0
        gbs = self.params.gbs
        rub_sf = self.params.rub_sf 
        
        # rate of PEP carboxylation
        Vp = self.calc_pep_carboxylation_rate(ci, Vpmax, Kp)
        
        
        
        # Cm assumed to equal Ci
        Cm = ci
       
        # Quadratic solution for enzyme limited C4 assimilation
        Ac = [None] * 2
        for k in am, pm:
            a = 1 - (alpha * Kc[k]) / (0.047 * Ko[k])
            b = -((Vp[k] - Rm[k] + gbs * Cm[k]) + 
                  (Vcmax[k] - Rd[k]) + gbs * Km[k] +
                  (alpha / 0.047 * (rub_sf * Vcmax[k] + Rd[k] * Kc[k] / Ko[k])))
            
            c = ( (Vcmax[k] - Rd[k]) * (Vp[k] - Rm[k] + gbs * Cm[k]) - 
                  (Vcmax[k] * gbs * rub_sf * Oi + Rd[k] * gbs * Km[k]) )
        
            Ac[k] = (-b - sqrt(b**2 - 4.0 * a * c)) / (2.0 * a)
		
        return Ac
    
    
    def calc_light_limited_assim(self, par, jmax, ci, Rd, Rm):
        """
        Calculate the AM/PM electron-transport-limited CO2 assimilation rate.
        
        Parameters:
        ----------
        x : float
            Partitioning factor of electron transport
        rub_sf : float 
            Half the reciprocal for Rubisco specificity 
        gbs : float
            bundle sheath conductance
        theta : float
             curvature factor
        alpha : float
            Fraction of PSII activity in the bundle sheath [0-1]
        f : float
            correction factor for spectral quality of light
        
        Returns:
        -------
        Aj : float, list [am, pm]
            Light-limited assimilation rate [umol m-2 s-1]
        """
        
        # local parameters
        am, pm = self.am, self.pm # morning/afternoon
        alpha = self.params.alpha_psii		
        theta = self.params.theta
        x = self.params.xpart_j
        #Oi = self.params.Oi / 1000.0 # mbar
        Oi = 210.0
        rub_sf = self.params.rub_sf 
        gbs = self.params.gbs 
        f = self.params.fspec
        labs = self.params.labs
        
        # useful light absorbed by photosystem II (PSII)
        light_abs = par * labs * (1.0 - f) / 2.0
        
        # Cm assumed to equal Ci
        Cm = ci
        
        # Quadratic solution for light-limited C4 assimilation
        Aj = [None] * 2
        for k in am, pm:
            # Total electron transport rat (umol electrons m-2 s-1)
            Jt = ((1.0 / (2.0 * theta)) * (light_abs + jmax[k] - 
                   sqrt((light_abs + jmax[k])**2.0 - 
                   4.0 * theta * light_abs * jmax[k])))
        
        
            a = 1.0 - ((7.0 * rub_sf * alpha) / (3.0 * 0.047))
            b = -( ((x * Jt) / 2.0 - Rm[k] + gbs * Cm[k]) + 
                   ((1.0 - x) * Jt / 3.0 - Rd[k]) + 
                   (gbs * (7.0 * rub_sf * Oi / 3.0)) + 
                   (alpha * rub_sf / 0.047) * 
                   (((1.0 - x) * Jt / 3.0) + (7.0 * Rd[k] / 3.0)) )
            c = ( ((x * Jt / 2.0 - Rm[k] + gbs * Cm[k]) *
                   ((1.0 - x) * Jt / 3.0 - Rd[k])) -
                   (gbs * rub_sf * Oi) * 
                   ((1.0 - x) * Jt / 3.0 + (7.0 * Rd[k] / 3.0)) )
                   
            Aj[k] = (-b - sqrt(b**2 - 4.0 * a * c)) / (2.0 * a)
        
        return Aj
    
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
