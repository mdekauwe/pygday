#!/usr/bin/env python
""" Photosynthesis model """

import sys
import datetime
import math
import constants as const
from utilities import float_eq, float_lt, float_le, float_gt, float_ge
from misc_funcs import day_length

__author__  = "Martin De Kauwe"
__version__ = "1.0 (08.03.2011)"
__email__   = "mdekauwe@gmail.com"


class Bewdy(object): 
    """ BEWDY - calculates plant C assimilation.
        
    Mechanistic model of gross canopy photosynthesis (GPP) as a function of 
    LAI, intensity of direct radiation, intensity of diffuse radiation at 
    the top of the canopy and leaf N content at the top of the canopy.
    See Medlyn et al for more details.
    
    The model by using a sunlit/shaded approach is a 
    compromise between a more detailed canopy model that models vertical 
    (and horizontal?) incident radiation and a big-leaf (no diffuse/direct
    distinction). I guess this is therefore effectively a "two-leaf" model.
    
    * should we have 2 different temps for sunlit and shaded leaves?
    
    References:
    -----------
    * Medlyn, B. E. et al (2000) Canadian Journal of Forest Research,30,873-888.
    
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
        
    def calculate_photosynthesis(self, frac_gcover, date, day, daylen):
        """
        Parameters:
        -----------
        frac_gcover : float
            fraction of ground cover
        date : date object string
            date string (yr/mth/day)
        day : integer
            day of simulation
        daylen : float
            daylength in hours
        
        Returns:
        --------
        Nothing 
            Method calculates GPP, NPP and Ra.
        """
        temp, sw_rad, ca = self.get_met_data(day) 
        
        # presumably conversion to seconds 
        daylength = daylen * const.HRS_TO_SECS * 1E-6
        
        # calculated from the canopy-averaged leaf N
        leaf_absorptance = ((self.state.ncontent / 2.8) / 
                            (self.state.ncontent / 2.8 + 0.076))
        
        direct_rad = sw_rad / daylength / 0.235 * self.params.direct_frac
        diffuse_rad = (sw_rad / daylength / 0.235 * 
                        (1.0 - self.params.direct_frac))
        
        (quantum_yield, rho) = self.calculate_bewdy_params(temp, ca)
        
        b = quantum_yield * self.params.kext * direct_rad * leaf_absorptance
        s = quantum_yield * self.params.kext * diffuse_rad * leaf_absorptance
        
        # effect of incomplete ground cover - modifies lai to lai/cover 
        # (jackson & palmer)
        lai_mod = self.state.lai / frac_gcover
        n = (rho * self.params.kext * 
                (self.state.ncontent - self.params.nmin) * lai_mod / 
                (1.0 - math.exp(-self.params.kext * lai_mod)))
        
        self.fluxes.gpp = ((self.sunlit_contribution(b, s, n, lai_mod) + 
                            self.shaded_contribution(b, s, n, lai_mod))) 
        
        # rescale gpp 
        self.fluxes.gpp *= daylength 
        
        self.fluxes.npp = self.calculate_npp(temp, frac_gcover)
        self.fluxes.npp_gCm2 = (self.fluxes.npp * const.M2_AS_HA / 
                                const.G_AS_TONNES)
        self.fluxes.gpp_gCm2 = self.fluxes.npp_gCm2 / self.params.cue
        #print self.fluxes.npp
        self.fluxes.apar = -999.9
        
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
        sw_rad : float
            SW down radiation [mj/m2/day]
        ca : float
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated.
        
        """
        temp = [ self.met_data['tam'][day], self.met_data['tpm'][day] ]
        sw_rad = self.met_data['sw_rad'][day]   
        
        if self.control.co2_conc == 0:
            ca = self.met_data['amb_co2'][day]
        elif self.control.co2_conc == 1:
            ca = self.met_data['ele_co2'][day]
        
        return temp, sw_rad, ca
        
    def calculate_bewdy_params(self, temp, ca):                  
        """ Calculates BEWDY model parameters
        
        Estimate the quantum yield (alpha) of absorbed radiation and rho, the 
        linear relationship between photosynthesis and leaf N content, using 
        the Farquhar and von Caemmerer (1982) model of leaf photosynthesis 
        In this model leaf photosysnthesis is given by the minimum of the rate 
        of carboxylation when Rubiso is limiting (ac) and the rate of 
        carboxylation when RUBP regeneration is limiting (aj).
        
        Temperature is assumed to affect km, gamma_star, jmax and vcmax.
        may want to use harley - type self.fluxes.temperature functions? However 
        these temperature dependences are taken from Brooks and Farquahr (1985) 
        for gamma_star, McMurtrie and Wang (1993) for km, and Kirschbaum (1986) 
        for jmax and vcmax.
        
        Parameters:
        -----------
        temp : float
            air temperature [degC]
        ca : float
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated. [umol mol-1]
        
        Returns:
        --------
        quantum_yield : float
            quantum yield of absorbed radiation
        rho : float
            model parameter
        """    
        # co2 compensation point in the absence of mitochondrial respiration
        gamma_star = (42.7 + 1.68 * (temp - 25.0) + 0.012 * (temp - 25.0)**2)
        
        # effective Michaelis-Menen coefficent of Rubisco activity
        km = 39.05 * math.exp(0.085 * temp) + 9.58 * gamma_star
        
        # max rate of electron transport and rubisco activity
        jmax, vcmax = self.jmax_and_vcmax_func(temp)
        
        # co2 concentration of intercellular air spaces
        ci = self.intercellular_co2_conc(gamma_star, ca)
        
        # quantum yield of absorbed radiation
        quantum_yield = self.calculate_quantum_yield(ci, gamma_star)
        
        # rate of carboxylation when rubiusco is limiting (Ac) 
        aj = ((jmax / 4.0) * ((ci - gamma_star) / (ci + 2. * gamma_star)) * 
                const.MOLE_C_TO_GRAMS_C)
        
        # rate of carboxylation when RUBP regeneration is limiting (Aj)
        ac = vcmax * ((ci - gamma_star) / (ci + km)) * const.MOLE_C_TO_GRAMS_C
        rho = min(ac, aj)
        
        return quantum_yield, rho
   
    def sunlit_contribution(self, b, s, n, lai_mod):
        """Calculate contribution from sunlit foliage
        
        Parameters:
        -----------
        b : float
            model parameter
        s : float
            model parameter
        n : float
            model parameter
        lai_mod : float
            effect of incomplete ground cover - modifies lai to lai/cover         
        
        Returns:
        --------
        sunlit_contribution : float
            contribution from sunlit foliage to GPP
        """
    
        arg1 = (1.0 / self.params.kext * 
                (1.0 - math.exp(-self.params.kext * lai_mod)))
        arg2 = (n * s * (n + s) + b * n**2) / (n + s)**2

        return arg1 * arg2
    
    def shaded_contribution(self, b, s, n, lai_mod):
        """Calculate contribution from shaded foliage
        
        Parameters:
        -----------
        b : float
            model parameter
        s : float
            model parameter
        n : float
            model parameter
        lai_mod : float
            effect of incomplete ground cover - modifies lai to lai/cover         
        
        Returns:
        --------
        shaded_contribution : float
            contribution from shaded foliage to GPP
        """
       
        arg1 = 1.0 / self.params.kext * (b**2 * n**2) / (n + s)**3.0
        arg2 = (math.log(((n + s) * math.exp(-self.params.kext * 
                    lai_mod) + b) / (n + s + b)))
       
        return arg1 * arg2
       
    
    def calculate_npp(self, temp, frac_gcover):
        """ figure out net photosynthesis 
        
        Parameters:
        -----------
        temp : float
            air temperature
        frac_gcover : float
            fraction of ground cover
        
        Returns:
        --------
        npp : float
            net primary productivity
        """
        
        if self.control.model_number == 5:
            # use dependence on nitrogen and temperature
            self.fluxes.auto_resp = self.calc_autotrophic_respiration(temp)
            
            npp = (self.params.growth_efficiency * 
                    (self.fluxes.gpp * frac_gcover * const.G_M2_2_TON_HEC - 
                    self.fluxes.auto_resp))
        else:
            # use proportionality with GPP
            npp = (self.params.cue * self.fluxes.gpp * frac_gcover * 
                    const.G_M2_2_TON_HEC)
        
        return npp
        
    def calc_autotrophic_respiration(self, temp):
        """Calculate respiration with dependence on N and temperature
        
        Parameters:
        -----------
        temp : float
            air temperature
        Returns:
        --------
        ra : float
            autotrophic respiration
        
        """
        plantn = self.state.shootn + self.state.rootn + self.state.stemnmob
        ra = (0.0106 * plantn * 12.0 / 14.0 * 
                                math.exp(self.params.kq10 * (temp - 15.0)))    
        return ra 

    def jmax_and_vcmax_func(self, temp):
        """ Maximum rate of electron transport (jmax) and of rubisco activity
        
        Parameters:
        -----------
        temp : float
            air temperature
        
        Returns:
        --------
        jmax : float
            maximum rate of electron transport
        vcmax : float
            maximum rate of Rubisco activity
        """
    
        if float_gt(temp, 10.0):
            jmax = (self.params.jmaxn * (1. + (temp - 25.0) * (0.05 + 
                    (temp - 25.0) * (-1.81 * 1E-3 + (temp - 25.0) * 
                    (-1.37 * 1E-4)))))       
            vcmax = (self.params.vcmaxn * (1.0 + (temp - 25.0) * 
                    (0.0485 + (temp - 25.0) * (-6.93 * 1E-4 + (temp - 25.0) * 
                    (-3.9 * 1E-5)))))
        elif float_gt(temp, 0.0):
            jmax = self.params.jmaxn * 0.0305 * temp
            vcmax = self.params.vcmaxn * 0.0238 * temp
        else:
            jmax = 0.0
            vcmax = 0.0
        
        return jmax, vcmax
   
    
    def intercellular_co2_conc(self, gamma_star, ca):
        """ Calculate Ci, intercellular CO2 concentration 
        
        Parameters:
        -----------
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration
        ca : float
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated.
            
        Returns:
        --------
        ci : float
            intercellular CO2 concentration. 
        """
        if self.control.use_leuning == 1:
            ci = (self.params.ambient_co2 - 
                    (self.params.ambient_co2 - gamma_star) * 
                    (1.0 + self.params.vpd / self.params.d0) * 1.6 / 
                    self.params.a1)
        else:
            # assume CO2 conc in the intercellular air spaces, ci is a constant
            # fraction of the atmospheric CO2, Ca.
            ci = self.params.ci_ca_ratio * ca
            
        return ci
    
    def calculate_quantum_yield(self, ci, gamma_star):
        """co2 fixed / photons absorbed 
        Parameters:
        -----------
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration
        ci : float
            intercellular CO2 concentration.
        
        Returns:
        --------
        alpha_a : float
            model_parameter
        
        """
        
        arg1 = self.params.alpha_j / 4.0
        arg2 = ((ci - gamma_star) / (ci + 2. * gamma_star))
        
        return arg1 * arg2 * const.MOLE_C_TO_GRAMS_C
        
    
    

 
if __name__ == "__main__":
    
    
    from file_parser import initialise_model_data
    from misc_funcs import day_length
    from utilities import float_lt
    import datetime
    
    default_dir = "/Users/mdekauwe/research/NCEAS_face/GDAY_duke_simulation/params"
    fname = "dk_varyco2_varyndep_grassequilib_then_forest_dukegrass_youngforest"
    
    (control, params, 
            state, files, 
            fluxes, met_data,
            print_opts) = initialise_model_data(fname, 
                                            default_dir=default_dir, DUMP=False) 
    
    B = Bewdy(control, params, state, fluxes, met_data)
    
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
    
    params.co2_offset = 0.0
    
    for project_day in xrange(len(met_data['prjday'])):
        
        state.shootnc = state.shootn / state.shoot 
        state.ncontent = (state.shootnc * params.cfracts / 
                                state.sla * const.KG_AS_G)
        daylen = day_length(datex, params.latitude)
        
        if float_lt(state.lai, params.lai_cover):
            frac_gcover = state.lai / params.lai_cover
        else:
            frac_gcover = 1.0
        
        B.calculate_photosynthesis(frac_gcover, datex, project_day, daylen)
        
        
        datex += datetime.timedelta(days=1)