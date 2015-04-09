#!/usr/bin/env python
""" 
Plot responses from the MATE photosynthesis model to PAR, Ca, SW and temp
"""

from math import exp, sqrt, sin, pi
from gday import constants as const
from gday.utilities import float_eq, float_gt, float_lt
from gday.mate import MateC3
import sys

__author__ = "Martin De Kauwe"
__version__ = "1.0 (04.03.2014)"
__email__  = "mdekauwe@gmail.com"


class Mate(MateC3):
    
    def calculate_photosynthesis(self, tair, par, vpd, ca, daylen, sw):
        
        
        Tk_am = tair + const.DEG_TO_KELVIN
        Tk_pm = tair + const.DEG_TO_KELVIN
        vpd_am = vpd
        vpd_pm = vpd
        
        state.wtfac_root = calc_sw_modifier(sw, self.params.ctheta_root, 
                                            self.params.ntheta_root)
        
       # calculate mate params & account for temperature dependencies
        N0 = self.calculate_top_of_canopy_n()
        
        gamma_star_am = self.calculate_co2_compensation_point(Tk_am)
        gamma_star_pm = self.calculate_co2_compensation_point(Tk_pm)
        
        Km_am = self.calculate_michaelis_menten_parameter(Tk_am)
        Km_pm = self.calculate_michaelis_menten_parameter(Tk_pm)
        
        (jmax_am, vcmax_am) = self.calculate_jmax_and_vcmax(Tk_am, N0)
        (jmax_pm, vcmax_pm) = self.calculate_jmax_and_vcmax(Tk_pm, N0)
        
        ci_am = self.calculate_ci(vpd_am, ca) 
        ci_pm = self.calculate_ci(vpd_pm, ca) 
        
        # quantum efficiency calculated for C3 plants
        alpha_am = self.calculate_quantum_efficiency(ci_am, gamma_star_am)
        alpha_pm = self.calculate_quantum_efficiency(ci_pm, gamma_star_pm)
        
        # Reducing assimilation if we encounter frost. Frost is assumed to 
        # impact on the maximum photosynthetic capacity and alpha_j
        # So there is only an indirect effect on LAI, this could be changed...
        if self.control.frost:
            Tmax = self.met_data['tmax'][day]
            Tmin = self.met_data['tmin'][day]
            
            Thard = self.calc_frost_hardiness(daylen, Tmin, Tmax)
            (total_alpha_limf, 
            total_amax_limf) = self.calc_frost_impact_factors(Thard, Tmin, Tmax)
            alpha_am *= total_alpha_limf
            alpha_pm *= total_alpha_limf
     
        # Rubisco carboxylation limited rate of photosynthesis
        ac_am = self.assim(ci_am, gamma_star_am, a1=vcmax_am, a2=Km_am) 
        ac_pm = self.assim(ci_pm, gamma_star_pm, a1=vcmax_pm, a2=Km_pm) 
        
        # Light-limited rate of photosynthesis allowed by RuBP regeneration
        aj_am = self.assim(ci_am, gamma_star_am, a1=jmax_am/4.0,
                           a2=2.0*gamma_star_am)
        aj_pm = self.assim(ci_pm, gamma_star_pm, a1=jmax_pm/4.0,
                           a2=2.0*gamma_star_pm)
        
        # light-saturated photosynthesis rate at the top of the canopy (gross)
        asat_am = min(aj_am, ac_am) 
        asat_pm = min(aj_pm, ac_pm) 
        if self.control.frost:
            asat_am *= total_amax_limf
            asat_pm *= total_amax_limf
            
        # LUE (umol C umol-1 PAR)
        lue_am = self.epsilon(asat_am, par, daylen, alpha_am)
        lue_pm = self.epsilon(asat_pm, par, daylen, alpha_pm)
        # use average to simulate canopy photosynthesis
        lue_avg = (lue_am + lue_pm) / 2.0 
        
        if float_eq(self.state.lai, 0.0):
            self.fluxes.apar = 0.0
        else:
            # absorbed photosynthetically active radiation (umol m-2 s-1)
            self.fluxes.apar = par * self.state.fipar
        apar_half_day = self.fluxes.apar / 2.0
        
        
        # convert umol m-2 d-1 -> gC m-2 d-1
        self.fluxes.gpp_gCm2 = self.fluxes.apar * lue_avg * const.UMOL_2_GRAMS_C
        self.fluxes.gpp_am = apar_half_day * lue_am * const.UMOL_2_GRAMS_C
        self.fluxes.gpp_pm = apar_half_day * lue_pm * const.UMOL_2_GRAMS_C
        
        # g C m-2 to tonnes hectare-1 day-1
        self.fluxes.gpp = self.fluxes.gpp_gCm2 * const.GRAM_C_2_TONNES_HA
        # photosynthesis -> umol m-2 s-1
        A = self.fluxes.apar * lue_avg / (daylen * 3600.)
        return A
            
def calc_sw_modifier(theta, c_theta, n_theta):
    """ From Landsberg and Waring """
    return 1.0  / (1.0 + ((1.0 - theta) / c_theta)**n_theta)

       
if __name__ == "__main__":
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    from gday.file_parser import read_met_forcing
    import gday.default_control as control
    import gday.default_files as files
    import gday.default_params as params
    import gday.default_fluxes as fluxes
    import gday.default_state as state
    
    # set up parameters
    params.alpha_j           = 0.26    
    params.cfracts           = 0.5     
    params.cue               = 0.5     
    params.delsj             = 644.4338
    params.eac               = 79430.0 
    params.eao               = 36380.0 
    params.eag               = 37830.0 
    params.eaj               = 43790.0 
    params.eav               = 51560.0 
    params.edj               = 2e+05   
    params.gamstar25         = 42.75   
    params.jmaxna            = 40.462  
    params.jmaxnb            = 13.691  
    params.kc25              = 404.9   
    params.ko25              = 278400.0
    params.measurement_temp  = 25.0    
    params.oi                = 205000.0
    params.theta             = 0.7     
    params.vcmaxna           = 20.497  
    params.vcmaxnb           = 8.403   
    params.g1                = 4.8 
    params.ctheta_root       = 0.4
    params.ntheta_root       = 3.0
    
    
    daylen = 12.0
    N = 100  
    par = np.linspace(0, 1500 * 3600. * daylen, N)
    tair = np.linspace(2, 40, N)
    vpd = np.linspace(0.5, 6, N)
    sw = np.linspace(0.0, 1.0, N)
    ca = np.linspace(250, 1000, N)

    met_fname = "met_data/DUKE_met_data_amb_co2.csv"
    met_data = read_met_forcing(met_fname, met_header=4)
    M = Mate(control, params, state, fluxes, met_data)
    
    # A vs PAR
    Alow_par = np.zeros(N)
    Ahigh_par = np.zeros(N)
    for i in xrange(N):
        
        state.ncontent = 4.5
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
        
        Alow_par[i] = M.calculate_photosynthesis(tair=25.0, par=par[i], vpd=1.0, 
                                          ca=380.0, daylen=daylen, sw=1.0)
        
        
        state.ncontent = 6.0
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
                                     
        Ahigh_par[i] = M.calculate_photosynthesis(tair=25.0, par=par[i], vpd=1.0, 
                                          ca=380.0, daylen=daylen, sw=1.0)
    
    
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(par/2600/daylen, Alow_par, "b-", label="Ncontent=4.5")
    ax.plot(par/2600/daylen, Ahigh_par, "r-", label="Ncontent=6.0")
    ax.set_ylabel('A ($\mu$mol m$^{-2}$ s$^{-1}$)')
    ax.set_xlabel('PAR')
    ax.legend(numpoints=1, loc='best', shadow=True).draw_frame(True)
    plt.show()
    
    # A vs temp
    par = 1500.0 * 3600.0 * daylen    
    Alow_tair = np.zeros(N)
    Ahigh_tair = np.zeros(N)
    for i in xrange(N):
        
        state.ncontent = 4.5
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
        
        Alow_tair[i] = M.calculate_photosynthesis(tair=tair[i], par=par, vpd=1.0, 
                                          ca=380.0, daylen=daylen, sw=1.0)
        
        
        state.ncontent = 6.0
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
                                     
        Ahigh_tair[i] = M.calculate_photosynthesis(tair=tair[i], par=par, vpd=1.0, 
                                          ca=380.0, daylen=daylen, sw=1.0)
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(tair, Alow_tair, "b-", label="Ncontent=4.5")
    ax.plot(tair, Ahigh_tair, "r-", label="Ncontent=6.0")
    ax.set_ylabel('A ($\mu$mol m$^{-2}$ s$^{-1}$)')
    ax.set_xlabel('Tair (deg C)')
    ax.legend(numpoints=1, loc='best', shadow=True).draw_frame(True)
    plt.show()
    
    # A vs Ca
    par = 1500.0 * 3600.0 * daylen    
    Alow_ca = np.zeros(N)
    Ahigh_ca = np.zeros(N)
    for i in xrange(N):
        
        state.ncontent = 4.5
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
        
        Alow_ca[i] = M.calculate_photosynthesis(tair=25.0, par=par, vpd=1.0, 
                                          ca=ca[i], daylen=daylen, sw=1.0)
        
        
        state.ncontent = 6.0
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
                                     
        Ahigh_ca[i] = M.calculate_photosynthesis(tair=25.0, par=par, vpd=1.0, 
                                          ca=ca[i], daylen=daylen, sw=1.0)
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(ca, Alow_ca, "b-", label="Ncontent=4.5")
    ax.plot(ca, Ahigh_ca, "r-", label="Ncontent=6.0")
    ax.set_ylabel('A ($\mu$mol m$^{-2}$ s$^{-1}$)')
    ax.set_xlabel('Ca ')
    ax.legend(numpoints=1, loc='best', shadow=True).draw_frame(True)
    plt.show()
    
    
    # A vs moisture
    par = 1500.0 * 3600.0 * daylen    
    Alow_sw = np.zeros(N)
    Ahigh_sw = np.zeros(N)
    for i in xrange(N):
        
        state.ncontent = 4.5
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
        
        Alow_sw[i] = M.calculate_photosynthesis(tair=25.0, par=par, vpd=1.0, 
                                          ca=380., daylen=daylen, sw=sw[i])
        
        
        state.ncontent = 6.0
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
                                     
        Ahigh_sw[i] = M.calculate_photosynthesis(tair=25.0, par=par, vpd=1.0, 
                                          ca=380., daylen=daylen, sw=sw[i])
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(sw*100, Alow_sw, "b-", label="Ncontent=4.5")
    ax.plot(sw*100, Ahigh_sw, "r-", label="Ncontent=6.0")
    ax.set_ylabel('A ($\mu$mol m$^{-2}$ s$^{-1}$)')
    ax.set_xlabel('Soil moisture availability (%)')
    ax.legend(numpoints=1, loc='best', shadow=True).draw_frame(True)
    plt.show()
    
    
    # A vs vpd
    par = 1500.0 * 3600.0 * daylen    
    Alow_vpd = np.zeros(N)
    Ahigh_vpd = np.zeros(N)
    for i in xrange(N):
        
        state.ncontent = 4.5
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
        
        Alow_vpd[i] = M.calculate_photosynthesis(tair=25.0, par=par, vpd=vpd[i], 
                                          ca=380., daylen=daylen, sw=1.0)
        
        
        state.ncontent = 6.0
        state.lai = 3.0
        state.fipar = (1.0 - exp(-params.kext * state.lai))
                                     
        Ahigh_vpd[i] = M.calculate_photosynthesis(tair=25.0, par=par, vpd=vpd[i], 
                                          ca=380., daylen=daylen, sw=1.0)
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(vpd, Alow_vpd, "b-", label="Ncontent=4.5")
    ax.plot(vpd, Ahigh_vpd, "r-", label="Ncontent=6.0")
    ax.set_ylabel('A ($\mu$mol m$^{-2}$ s$^{-1}$)')
    ax.set_xlabel('VPD (kPa)')
    ax.legend(numpoints=1, loc='best', shadow=True).draw_frame(True)
    plt.show()
    
    par = np.linspace(0, 1500, N)
    f = open("mate_stuff.csv", "w")
    print >> f, "par tair vpd sw ca alow_ca ahigh_ca alow_sw ahigh_sw alow_par ahigh_par alow_tair ahigh_tair" 
    for i in xrange(N):
        print >> f, par[i], tair[i], vpd[i], \
                    calc_sw_modifier(sw[i],params.ctheta_root,params.ntheta_root),\
                    ca[i], \
                    Alow_ca[i], Ahigh_ca[i],\
                    Alow_sw[i], Ahigh_sw[i],\
                    Alow_par[i], Ahigh_par[i],\
                    Alow_tair[i], Ahigh_tair[i]
    f.close()
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    