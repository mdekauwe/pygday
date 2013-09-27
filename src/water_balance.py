# -*- coding: UTF-8 -*-

from math import log, exp, sqrt, pi

from utilities import float_gt, float_eq, clip
import constants as const
import sys

__author__  = "Martin De Kauwe"
__version__ = "1.0 (02.05.2012)"
__email__   = "mdekauwe@gmail.com"

          
class WaterBalance(object):
    """Dynamic water balance model.

    Contains a few extra routinues to do with WUE calculation from MATE

    References:
    ===========
    * McMurtrie, R. (1990) Water/nutrient interactions affecting the
        productivity of stands of Pinus radiata. Forest Ecology and Management,
        30, 415-423.

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
        self.P = PenmanMonteith(canht=self.params.canht, 
                                dz0v_dh=self.params.dz0v_dh,
                                displace_ratio=self.params.displace_ratio,
                                z0h_z0m=self.params.z0h_z0m)
        self.am = 0
        self.pm = 1
            
    def calculate_water_balance(self, day, daylen):
        """ Calculate water balance

        Parameters:
        ----------
        day : int
            project day.
        daylen : float
            length of day in hours.

        """
        half_day = daylen/2.0
        (am, pm) = self.am, self.pm # morning/afternoon
        
        # met forcing
        (tair_ampm, tair_day, rain, sw_rad_ampm, 
         sw_rad_day, vpd_ampm, vpd_day, wind_ampm, 
         wind_day, ca, press) = self.get_met_data(day, daylen)
        
        net_rad_day = self.calc_radiation(tair_day, sw_rad_day, daylen)
        net_rad_ampm = [self.calc_radiation(tair_ampm[am], sw_rad_ampm[am], \
                                            half_day), \
                        self.calc_radiation(tair_ampm[pm], sw_rad_ampm[pm], \
                                            half_day)]
        # calculate water fluxes
        if self.control.trans_model == 0:
            # transpiration calculated from WUE...
            self.calc_transpiration()
        elif self.control.trans_model == 1:
        
            if self.control.assim_model == "BEWDY":
                self.calc_transpiration_penmon(vpd_day, net_rad_day, tair_day, 
                                               wind_day, ca, daylen, press)
            
            elif self.control.assim_model == "MATE":
                self.calc_transpiration_penmon_am_pm(net_rad_ampm, wind_ampm, 
                                                     ca, daylen, press, 
                                                     vpd_ampm, tair_ampm)
        
        elif self.control.trans_model == 2:
            self.calc_transpiration_priestay(net_rad_avg, tair_day, press)
    
        self.calc_infiltration(rain)
        self.fluxes.soil_evap = self.calc_soil_evaporation(tair_day, 
                                                           net_rad_day,
                                                           press, daylen, 
                                                           sw_rad_day)
        self.fluxes.et = (self.fluxes.transpiration + self.fluxes.soil_evap +
                          self.fluxes.interception)
        self.fluxes.runoff = self.update_water_storage()
       
    def get_met_data(self, day, daylen):
        """ Grab the days met data out of the structure and return day values.

        Parameters:
        ----------
        day : int
            project day.
        daylen : float
            length of day in hours.

        Returns:
        -------
        tavg : float
            average daytime temperature [degC]
        rain : float
            rainfall [mm d-1]
        sw_rad : float
            sw down radiation [mj m-2 day-1]
        vpd : float
            average daily vpd [kPa]
        ca : float
            atmospheric co2 [umol mol-1]
        wind : float
            average daily wind speed [mm s-1]
        press : float
            average daytime pressure [kPa]

        """
        am, pm = self.am, self.pm
        ca = self.met_data['co2'][day]
        tair_day = self.met_data['tair'][day]
        tair_ampm = [self.met_data['tam'][day], self.met_data['tpm'][day]]
        sw_rad_day = self.met_data['sw_rad'][day]
        sw_rad_ampm = [self.met_data['sw_rad_am'][day], \
                       self.met_data['sw_rad_pm'][day]]
        rain = self.met_data['rain'][day]
        vpd_day = self.met_data['vpd_avg'][day] # daytime average
        vpd_ampm= [self.met_data['vpd_am'][day], self.met_data['vpd_pm'][day]]
        wind_ampm = [self.met_data['wind_am'][day], \
                     self.met_data['wind_pm'][day]]
        wind_day = self.met_data['wind'][day]
        
        if ('atmos_press' in self.met_data and not
            self.met_data['atmos_press'] is None):
            press = self.met_data['atmos_press'][day]
        else:
            press = None # use method below to calculate pressure

        return (tair_ampm, tair_day, rain, sw_rad_ampm, sw_rad_day, vpd_ampm, 
                vpd_day,  wind_ampm, wind_day, ca, press)

    def calc_infiltration(self, rain):
        """ Estimate "effective" rain, or infiltration I guess.

        Simple assumption that infiltration relates to leaf area
        and therefore canopy storage capacity (wetloss). Interception is
        likely to be ("more") erroneous if a canopy is subject to frequent daily
        rainfall I would suggest.

        Parameters:
        -------
        rain : float
            rainfall [mm d-1]

        """
        if self.state.lai > 0.0:
            self.fluxes.erain = max(0.0, rain * self.params.rfmult -
                                    self.state.lai * self.params.wetloss)
            self.fluxes.interception = (rain * self.params.rfmult - 
                                        self.fluxes.erain)
        else:
            self.fluxes.erain = max(0.0, rain)
            self.fluxes.interception = 0.0
                
    def calc_transpiration(self):
        """ units mm/day """
        if float_gt(self.fluxes.wue, 0.0):
            self.fluxes.transpiration = self.fluxes.gpp_gCm2 / self.fluxes.wue
        else:
            self.fluxes.transpiration = 0.0

    def calc_transpiration_priestay(self, net_rad, tavg, press):
        """ Calculate canopy transpiration using the Priestley Taylor eqn
        units (mm/day)

        Parameters:
        -----------
        tavg : float
            average daytime temp [degC]
        net_rad : float
            net radiation [mj m-2 s-1]
        press : float
            average daytime pressure [kPa]

        """
        P = PriestleyTaylor()
        self.fluxes.transpiration = P.calc_evaporation(net_rad, tavg, press,
                                                        pt_coeff=1.26)

    def calc_transpiration_penmon(self, vpd, net_rad, tavg, wind, ca, daylen, 
                                  press):
        """ Calculate canopy transpiration using the Penman-Monteith equation.
        units mm/day

        Parameters:
        -----------
        vpd : float
            average daily vpd [kPa]
        net_rad : float
            net radiation [mj m-2 s-1]
        tavg : float
            average daytime temp [degC]
        wind : float
            average daily wind speed [m s-1]
        ca : float
            atmospheric co2 [umol mol-1]
        daylen : float
            daylength in hours
        press : float
            average daytime pressure [kPa]

        """
        
        gs, gs_mol = self.calc_stomatal_conductance(vpd, ca, daylen, 
                                                    self.fluxes.gpp_gCm2, 
                                                    press, tavg)
        transp, omegax = self.P.calc_evaporation(vpd, wind, gs, net_rad, tavg, 
                                                 press)
        
        self.fluxes.gs_mol_m2_sec = gs_mol
        ga = self.P.canopy_boundary_layer_conductance(wind)
        self.fluxes.ga_mol_m2_sec = ga / const.CONV_CONDUCT
       
        tconv = 60.0 * 60.0 * daylen # seconds to day
        self.fluxes.transpiration = transp * tconv
        
    def calc_transpiration_penmon_am_pm(self, net_rad, wind, ca, 
                                        daylen, press, vpd, tair):
        """ Calculate canopy transpiration using the Penman-Monteith equation
        using am and pm data [mm/day]
        
        Parameters:
        -----------
        vpd : float
            average daily vpd [kPa]
        net_rad_am : float
            net radiation [mj m-2 s-1] (morning)
        net_rad_pm : float
            net radiation [mj m-2 s-1] (afternoon)
        tair : float
            AM/PM air temp [degC] am/pm
        wind : float
            daily wind speed [m s-1]
        ca : float
            atmospheric co2 [umol mol-1]
        daylen : float
            daylength in hours
        press : float
            average daytime pressure [kPa]

        """
        gs_mol = make_empty_list_of_zeros(2)
        ga = make_empty_list_of_zeros(2)
        gs = make_empty_list_of_zeros(2) # m s-1 (but per half day)
        trans = make_empty_list_of_zeros(2)
        omegax = make_empty_list_of_zeros(2)
        gpp = self.fluxes.gpp_am_pm # list
        half_day = daylen / 2.0
        for i in self.am, self.pm:
            (gs[i], 
             gs_mol[i]) = self.calc_stomatal_conductance(vpd[i], ca, half_day, 
                                                         gpp[i], press, tair[i])
            
            ga[i] = self.P.canopy_boundary_layer_conductance(wind[i])
            
            (trans[i], 
             omegax[i]) = self.P.calc_evaporation(vpd[i], wind[i], gs[i], 
                                                  net_rad[i], tair[i], press, 
                                                  ga=ga[i])
        # print out pre-noon values
        #print self.fluxes.omega = omegax[0]
        self.fluxes.omega = sum(omegax) / 2.0                                  
        self.fluxes.gs_mol_m2_sec = sum(gs_mol) #/ 2.0
        
        #ga = P.canopy_boundary_layer_conductance(wind_avg)
        self.fluxes.ga_mol_m2_sec = (sum(ga) / 2.0) / const.CONV_CONDUCT
        
        # seconds to day, note daylength divided by 2, because fluxes were
        # calculated per half day
        tconv = 60.0 * 60.0 * half_day
        self.fluxes.transpiration = sum(trans) * tconv
        
    def calc_stomatal_conductance(self, vpd, ca, daylen, gpp, press, temp):
        """ Calculate stomatal conductance, note assimilation rate has been
        adjusted for water availability at this point.
        
        gs = g0 + 1.6 * (1 + g1/sqrt(D)) * A / Ca 

        units: m s-1 (conductance)
        References:
        -----------
        For conversion factor for conductance see...
        * Jones (1992) Plants and microclimate, pg 56 + Appendix 3
        * Diaz et al (2007) Forest Ecology and Management, 244, 32-40.

        Parameters:
        -----------
        vpd : float
            average daily vpd [kPa]
        ca : float
            atmospheric co2 [umol mol-1]
        daylen : float
            daylength in hours, passing half day as AM and PM

        Returns:
        --------
        gs : float
            stomatal conductance [m s-1]
        """
        # convert conductance to water vapour units
        g1_c = self.params.g1 / const.RATIO_DIFF_H2O_TO_CO2
       
        # time unit conversion day-1 -> seconds-1
        tconv =  1.0 / (60.0 * 60.0 * daylen)
        gpp_umol_m2_sec = (gpp * const.GRAMS_C_TO_MOL_C * const.MOL_TO_UMOL * 
                           tconv)
        
        arg1 = 1.6 * (1.0 + (g1_c * self.state.wtfac_root) / sqrt(vpd))
        arg2 = gpp_umol_m2_sec / ca # umol mol-1
        gs_mol_m2_sec = arg1 * arg2 * const.RATIO_DIFF_H2O_TO_CO2
        
        # convert to mm s-1 and then to m s-1
        #return (gs_mol_m2_sec * const.MOL_TO_MILLIMOLES * const.CONV_CONDUCT * 
        #        const.MM_TO_M)
        
        # See Jones, 1992, appendix
        tk = temp + const.DEG_TO_KELVIN
        conv = const.MM_TO_M / (press / (const.RGAS * tk))
        
        # convert to mm s-1 to m s-1
        return (gs_mol_m2_sec * conv, gs_mol_m2_sec) 
   
    
    def calc_radiation(self, tavg, sw_rad, daylen):
        """
        Estimate net radiation assuming 'clear' skies...

        References:
        -----------
        * Ritchie, 1972, Water Resources Research, 8, 1204-1213.
        * Monteith and Unsworth (1990) Principles of Environmental Physics.

        Parameters:
        -----------
        tavg : float
            average daytime temp [degC]
        sw_rad : float
            sw down radiation [mj m-2 d-1]
        daylen : float
            daylength in hours

        Returns:
        --------
        net_rad : float
            net radiation [mj m-2 s-1]

        """
        # Net loss of longwave radiation
        # Monteith and Unsworth '90, pg. 52, 54.
        net_lw = (107.0 - 0.3 * tavg) * daylen * const.WATT_HR_TO_MJ
        net_rad = max(0.0, sw_rad * (1.0 - self.params.albedo) - net_lw)
       
        # convert units for met data
        tconv = 1.0 / (60.0 * 60.0 * daylen)  # day-1 to seconds-1
        
        return net_rad * tconv # MJ m-2 s-1

    def calc_soil_evaporation(self, tavg, net_rad, press, daylen, sw_rad):
        """ Use Penman eqn to calculate top soil evaporation flux at the
        potential rate.

        Soil evaporation is dependent upon soil wetness and plant cover. The net
        radiation term is scaled for the canopy cover passed to this func and
        the impact of soil wetness is accounted for in the wtfac term. As the
        soil dries the evaporation component reduces significantly.

        Key assumptions from Ritchie...

        * When plant provides shade for the soil surface, evaporation will not
        be the same as bare soil evaporation. Wind speed, net radiation and VPD
        will all belowered in proportion to the canopy density. Following
        Ritchie role ofwind, VPD are assumed to be negligible and are therefore
        ignored.

        These assumptions are based on work with crops and whether this holds
        for tree shading where the height from the soil to the base of the
        crown is larger is questionable.

        units = (mm/day)

        References:
        -----------
        * Ritchie, 1972, Water Resources Research, 8, 1204-1213.

        Parameters:
        -----------
        tavg : float
            average daytime temp [degC]
        net_rad : float
            net radiation [mj m-2 day-1]
        press : float
            average daytime pressure [kPa]

        Returns:
        --------
        soil_evap : float
            soil evaporation [mm d-1]

        """
        P = Penman()
        soil_evap = P.calc_evaporation(net_rad, tavg, press)
        
        # Surface radiation is reduced by overstory LAI cover. This empirical
        # fit comes from Ritchie (1972) and is formed by a fit between the LAI
        # of 5 crops types and the fraction of observed net radiation at the
        # surface. Whilst the LAI does cover a large range, nominal 0–6, there
        # are only 12 measurements and only three from LAI > 3. So this might
        # not hold as well for a forest canopy?
        # Ritchie 1972, Water Resources Research, 8, 1204-1213.
        soil_evap *= exp(-0.398 * self.state.lai)
        
        # reduce soil evaporation if top soil is dry
        soil_evap *= self.state.wtfac_topsoil
        tconv = 60.0 * 60.0 * daylen # seconds to day
        
        return soil_evap * tconv
        
        
    def update_water_storage(self, tolerance=1E-08):
        """ Calculate root and top soil plant available water and runoff.
        
        Soil drainage is estimated using a "leaky-bucket" approach with two
        soil layers. In reality this is a combined drainage and runoff 
        calculation, i.e. "outflow". There is no drainage out of the "bucket" 
        soil. 
        
        Returns:
        --------
        outflow : float
            outflow [mm d-1]
        """
        # reduce transpiration from the top soil if it is dry
        trans_frac = (self.params.fractup_soil * self.state.wtfac_topsoil)
        
        # Total soil layer
        self.state.pawater_topsoil += (self.fluxes.erain -
                                    (self.fluxes.transpiration *
                                     trans_frac) -
                                     self.fluxes.soil_evap)
        
        self.state.pawater_topsoil = clip(self.state.pawater_topsoil, min=0.0,
                                        max=self.params.wcapac_topsoil) 
        
        # Total root zone
        previous = self.state.pawater_root
        self.state.pawater_root += (self.fluxes.erain -
                                    self.fluxes.transpiration -
                                    self.fluxes.soil_evap)
        
        # calculate runoff and remove any excess from rootzone
        if self.state.pawater_root > self.params.wcapac_root:
            runoff = self.state.pawater_root - self.params.wcapac_root
            self.state.pawater_root -= runoff 
        else:
            runoff = 0.0
        
        self.state.pawater_root = clip(self.state.pawater_root, min=0.0,
                                       max=self.params.wcapac_root)
        
        self.state.delta_sw_store = self.state.pawater_root - previous
        
        return runoff 


class SoilMoisture(object):
    """ Estimate current soil moisture factor 
    
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
    
    References:
    -----------
    * Cosby et al. (1984) Water Resources Research, 20, 682-690.
    """
    def __init__(self, control, params, state, fluxes):
        
        self.params = params
        self.fluxes = fluxes
        self.control = control
        self.state = state
        self.silt_index = 0
        self.sand_index = 1
        self.clay_index = 2
    
    def initialise_parameters(self):
        # initialise parameters, if these are not known for the site use
        # values derived from Cosby et al to calculate the amount of plant
        # available water.
        
        # local variable
        topsoil_type = self.params.topsoil_type
        rootsoil_type = self.params.rootsoil_type

        if self.control.calc_sw_params:
            fsoil_top = self.get_soil_fracs(topsoil_type)
            fsoil_root = self.get_soil_fracs(rootsoil_type)  
            
            # topsoil
            (self.theta_fc_tsoil, 
             self.theta_wp_tsoil) = self.calc_soil_params(fsoil_top)
            
            # Plant available water in top soil (mm)
            self.params.wcapac_topsoil = (self.params.topsoil_depth * 
                                         (self.theta_fc_tsoil - 
                                          self.theta_wp_tsoil))
            
            # Rootzone
            (self.theta_fc_root, 
             self.theta_wp_root) = self.calc_soil_params(fsoil_root)
            
            # Plant available water in rooting zone (mm)
            self.params.wcapac_root = (self.params.rooting_depth * 
                                      (self.theta_fc_root - 
                                       self.theta_wp_root))
        
        # calculate Landsberg and Waring SW modifier parameters if not
        # specified by the user based on a site calibration
        if (self.params.ctheta_topsoil is None and 
            self.params.ntheta_topsoil is None and
            self.params.ctheta_root is None and 
            self.params.ntheta_root is None):      
           
            (self.params.ctheta_topsoil, 
             self.params.ntheta_topsoil) = self.get_soil_params(topsoil_type)
            
            (self.params.ctheta_root, 
             self.params.ntheta_root) = self.get_soil_params(rootsoil_type)  
        """
        print self.params.rooting_depth
        print self.params.wcapac_topsoil
        print self.params.wcapac_root
        print 
        
        print "theta_wilt", self.theta_wp_root
        print "theta_field_capacity", self.theta_fc_root
        
        print
    
        print "theta_wilt", self.theta_wp_root * self.params.wcapac_root
        print "theta_field_capacity", self.theta_fc_root * self.params.wcapac_root
        
        print 
        print self.params.ctheta_topsoil
        print self.params.ntheta_topsoil
        print self.params.ctheta_root
        print self.params.ntheta_root
        
        sys.exit()
        """
        
    def get_soil_params(self, soil_type):
        """ For a given soil type, get the parameters for the soil
        moisture availability based on Landsberg and Waring.
        
        Reference
        ---------
        * Landsberg and Waring (1997) Forest Ecology & Management, 95, 209-228.
        * updated with additional values based on Ward et al. 2000 cited in
          Feikema et al (2010) Description of the 3PG+ forest growth model
         """
        soil_types = ["sand", "sandy_loam", "clay_loam", "clay"]
        fsoil = None
        if soil_type == "sand":
            c_theta = 0.75
            n_theta = 10.0
        elif soil_type == "loamy_sand":
            c_theta = 0.7
            n_theta = 9.0
        elif soil_type == "sandy_loam":
            c_theta = 0.65
            n_theta = 8.0
        elif soil_type == "loam":
            c_theta = 0.6
            n_theta = 7.0
        elif soil_type == "silty_loam":
            c_theta = 0.55
            n_theta = 6.0
        elif soil_type == "silty_clay_loam":
            c_theta = 0.5
            n_theta = 5.0
        elif soil_type == "clay_loam":
            c_theta = 0.45
            n_theta = 4.0
        elif soil_type == "sandy_clay":
            c_theta = 0.4
            n_theta = 3.0
        elif soil_type == "silty_clay":
            c_theta = 0.35
            n_theta = 2.0
        elif soil_type == "clay":
            c_theta = 0.3
            n_theta = 1.0
        else:
            print 'There are no parameters for your soil type. Either use the'
            print 'other soil water stress model or specify the parameters.'
            sys.exit()
        return c_theta, n_theta   
        
       
    def get_soil_fracs(self, soil_type):
        """ Based on Table 2 in Cosby et al 1984, page 2."""
        fsoil = None
        if soil_type == "sand":
            fsoil = [0.05, 0.92, 0.03]
        elif soil_type == "loamy_sand":
            fsoil = [0.12, 0.82, 0.06]
        elif soil_type == "sandy_loam":
            fsoil = [0.32, 0.58, 0.1]
        elif soil_type == "loam":
            fsoil = [0.39, 0.43, 0.18]
        elif soil_type == "silty_loam":
            fsoil = [0.70, 0.17, 0.13]
        elif soil_type == "sandy_clay_loam":
            fsoil = [0.15, 0.58, 0.27]
        elif soil_type == "clay_loam":
            fsoil = [0.34, 0.32, 0.34]
        elif soil_type == "silty_clay_loam":
            fsoil = [0.56, 0.1, 0.34]
        elif soil_type == "sandy_clay":
            fsoil = [0.06, 0.52, 0.42]
        elif soil_type == "silty_clay":
            fsoil = [0.47, 0.06, 0.47]
        elif soil_type == "clay":
            fsoil = [0.2, 0.22, 0.58]
        else:
            print 'Could not understand soil type', soil_type
            sys.exit()
        return fsoil
    
    def calc_soil_params(self, fsoil):
        """ Cosby parameters for use within the Clapp Hornberger soil hydraulics
        scheme are calculated based on the texture components of the soil.
        
        NB: Cosby et al were ambiguous in their paper as to what log base to 
        use.  The correct implementation is base 10, as below.
        
        Parameters:
        ----------
        fsoil : list
            fraction of silt, sand, and clay (in that order
        
        Returns:
        --------
        theta_fc : float
            volumetric soil water concentration at field capacity
        theta_wp : float
            volumetric soil water concentration at the wilting point
            
        """
        pressure_head_wilt = 152.9
        pressure_head_crit = 3.364
        
        # Clapp Hornberger exponent
        b = 3.1 + 15.7 * fsoil[self.clay_index] - 0.3 * fsoil[self.sand_index] 
        
        # saturated soil water suction kg m-2 s-1
        sathh = (0.01 * 10.0**(2.17 - 0.63 * fsoil[self.clay_index] - 1.58 * 
                               fsoil[self.sand_index]))
        
        # volumetric soil moisture concentrations at the saturation point
        theta_sp = (0.505 - 0.037 * fsoil[self.clay_index] - 0.142 * 
                     fsoil[self.sand_index])
        
        # volumetric soil moisture concentrations at the wilting point
        # assumed to = to a suction of -1.5 MPa or a depth of water of 152.9 m
        theta_wp = theta_sp * (sathh / pressure_head_wilt)**(1.0 / b)
        
        # volumetric soil moisture concentrations at the critical point (field
        # capacity) assumed to equal a suction of -0.033 MPa or a 
        # depth of water of 3.364 m
        theta_fc = theta_sp * (sathh / pressure_head_crit)**(1.0 / b)
        
        return (theta_fc, theta_wp)
    
    def calculate_soil_water_fac(self):
        """ Estimate a relative water availability factor [0..1]

        A drying soil results in physiological stress that can induce stomatal
        closure and reduce transpiration. Further, N mineralisation depends on 
        top soil moisture.
        
        self.params.qs = 0.2 in SDGVM
        
        References:
        -----------
        * Landsberg and Waring (1997) Forest Ecology and Management, 95, 209-228.
          See --> Figure 2.
        * Egea et al. (2011) Agricultural Forest Meteorology, 151, 1370-1384.
          
        But similarly see:
        * van Genuchten (1981) Soil Sci. Soc. Am. J, 44, 892--898.
        * Wang and Leuning (1998) Ag Forest Met, 91, 89-111.
        
        * Pepper et al. (2008) Functional Change Biology, 35, 493-508
       
        Returns:
        --------
        wtfac_topsoil : float
            water availability factor for the top soil [0,1]
        wtfac_root : float
            water availability factor for the root zone [0,1]    
        """
        # turn into fraction...
        smc_topsoil = self.state.pawater_topsoil / self.params.wcapac_topsoil
        smc_root = self.state.pawater_root / self.params.wcapac_root
        
        if self.control.sw_stress_model == 0:
            wtfac_topsoil = smc_topsoil**self.params.qs  
            wtfac_root = smc_root**self.params.qs  
            
        elif self.control.sw_stress_model == 1:
            wtfac_topsoil = self.calc_sw_modifier(smc_topsoil, 
                                                  self.params.ctheta_topsoil, 
                                                  self.params.ntheta_topsoil)
  
            wtfac_root = self.calc_sw_modifier(smc_root, 
                                               self.params.ctheta_root, 
                                               self.params.ntheta_root)
        return (wtfac_topsoil, wtfac_root) 
        
    def calc_sw_modifier(self, theta, c_theta, n_theta):
        """ From Landsberg and Waring """
        return 1.0  / (1.0 + ((1.0 - theta) / c_theta)**n_theta)

class PenmanMonteith(object):

    """ Water loss from a canopy (ET), representing surface as a big "leaf".
    The resistance to vapour transfer from the canopy to the atmosphere is
    determined by the aerodynamic canopy conductance (ga) and the stomatal 
    conductance (gs). If the surface is wet then there is a further water vapour
    flux from the soil/surface (calculated elsewhere!).

    Assumption is that calculation is for the entire stand (one surface), e.g. 
    the single-layer approach. Second major assumption is that soil heat is
    zero over the course of a day and is thus ignored.

    Value for cp comes from Allen et al 1998.

    units: mm day-1

    References:
    -----------
    * Monteith and Unsworth (1990) Principles of Environmental
      Physics, pg. 247. Although I have removed the soil heat flux as G'DAY calculates soil evaporation seperately.
    * Allen et al. (1989) Operational estimates of reference evapotranspiration.
      Agronomy Journal, 81, 650-662.
    * Allen et al. (1998) Crop evapotranspiration - Guidelines for computing
      crop water requirements - FAO Irrigation and drainage paper 56.
      http://www.fao.org/docrep/X0490E/x0490e00.htm#Contents. PDF in bibtex lib.
    * Harrison (1963) Fundamentals concepts and definitions relating to
      humidity. In Wexler, A. (Ed.) Humidity and moisture. Vol 3, Reinhold
      Publishing Co., New York, NY, USA.
    * Dawes and Zhang (2011) Waves - An integrated energy and water balance model
      http://www.clw.csiro.au/products/waves/downloads/chap3.pdf
    """

    def __init__(self, cp=1.013E-3, vk=0.41, epsilon=0.6222, zele_sea=125.0,
                    canht=20.0, dz0v_dh=0.1, displace_ratio=0.67, z0h_z0m=1.0):

        """
        Parameters:
        -----------
        cp : float
            specific heat of dry air [MJ kg-1 degC-1]
        vk : float
            von Karman's constant [unitless]
        epsilon : float
            ratio molecular weight of water vap/dry air
        zele_sea : float
            elevation above sea level [m]
        canht : float
            canopy height [m]
        dz0v_dh : float
            rate change of roughness for momentum with height
        displace_ratio : float
            zero plain displacement height
        z0h_z0m : float
            Ratio of the roughness length for heat to the roughness length for 
            momentum, see comment in method below!!!
        """

        self.cp = cp
        self.vk = vk
        self.epsilon = epsilon
        self.zele_sea = zele_sea
        self.J_TO_MJ = 1.0E-6
        self.C_TO_K = 273.15
        self.canht = canht
        self.dz0v_dh = dz0v_dh
        self.displace_ratio = displace_ratio # zero plan displacement height
        self.z0h_z0m = z0h_z0m
        
    def calc_evaporation(self, vpd, wind, gs, net_rad, tavg, press, ga=None):

        """
        Parameters:
        -----------
        vpd : float
            vapour pressure def [kPa]
        wind : float
            average daytime wind speed [m s-1]
        gs : float
            stomatal conductance [m s-1]
        net_rad : float
            net radiation [mj m-2 s-1] 
        tavg : float
            daytime average temperature [degC]
        press : float
            average daytime pressure [kPa]

        Returns:
        --------
        et : float
            evapotranspiration [mm d-1]

        """
        # if not read from met file calculate atmospheric pressure from sea lev
        if press == None:
            press = self.calc_atmos_pressure()
        
        lambdax = self.calc_latent_heat_of_vapourisation(tavg)
        gamma = self.calc_pyschrometric_constant(lambdax, press)
        slope = self.calc_slope_of_saturation_vapour_pressure_curve(tavg)
        rho = self.calc_density_of_air(tavg)
        if ga is None:
            ga = self.canopy_boundary_layer_conductance(wind)
       
        if float_gt(gs, 0.0):
            # decoupling coefficent, Jarvis and McNaughton, 1986
            # when omega is close to zero, it is said to be well coupled and
            # gs is the dominant controller of water loss (gs<ga).
            e = slope / gamma # chg of latent heat relative to sensible heat of air
            omega = (e + 1.0) / (e + 1.0 + (ga / gs))
            
            arg1 = ((slope * net_rad ) + (rho * self.cp * vpd * ga))
            arg2 = slope + gamma * (1.0 + (ga / gs))
            et = (arg1 / arg2) / lambdax
        else:
            et = 0.0
            omega = 0.0
        
        return et, omega

    def canopy_boundary_layer_conductance(self, wind):
        """  Canopy boundary layer conductance, ga or 1/ra

        Characterises the heat/water vapour from evaporating surface, but does 
        not account for leaf boundary layer conductance, which is the parellel 
        sum of single leaf boundary layer conductances for all leaves in the 
        canopy.

        Notes:
        ------
        'Estimates of ga for pine canopies from LAI of 3 to 6 vary from
        3.5 to 1.1 mol m-2 s-1  (Kelliher et al., 1993; Juang et al., 2007).'
        Drake et al, 2010, 17, pg. 1526.

        References:
        ------------
        * Jones 1992, pg. 67-8.
        * Monteith and Unsworth (1990), pg. 248. Note this in the inverted form
          of what is in Monteith (ga = 1 / ra)
        * Allen et al. (1989) pg. 651.
        * Gash et al. (1999) Ag forest met, 94, 149-158.

        Parameters:
        -----------
        wind : float
            average daytime wind speed [m s-1]

        Returns:
        --------
        ga : float
            canopy boundary layer conductance [m s-1]
        """
        # z0m roughness length governing momentum transfer [m]
        z0m = self.dz0v_dh * self.canht
    
        # z0h roughness length governing transfer of heat and vapour [m]
        # *Heat tranfer typically less efficent than momentum transfer. There is
        #  a lot of variability in values quoted for the ratio of these two...
        #  JULES uses 0.1, Campbell and Norman '98 say z0h = z0m / 5. Garratt 
        #  and Hicks, 1973/ Stewart et al '94 say z0h = z0m / 7. Therefore for 
        #  the default I am following Monteith and Unsworth, by setting the 
        #  ratio to be 1, the code below is identical to that on page 249, 
        #  eqn 15.7
        z0h = self.z0h_z0m * z0m
        
        # zero plan displacement height [m]
        d = self.displace_ratio * self.canht
        
        arg1 = self.vk**2 * wind
        arg2 = log((self.canht - d) / z0m)
        arg3 = log((self.canht - d) / z0h) 
        
        return arg1 / (arg2 * arg3)
        
        
    def calc_slope_of_saturation_vapour_pressure_curve(self, tavg):
        """ Eqn 13 from FAO paper, Allen et al. 1998.

        Parameters:
        -----------
        tavg : float
            average daytime temperature

        Returns:
        --------
        slope : float
            slope of saturation vapour pressure curve [kPa degC-1]

        """
        t = tavg + 237.3
        arg1 = 4098.0 * (0.6108 * exp((17.27 * tavg) / t))
        arg2 = t**2
        return (arg1 / arg2)

    def calc_pyschrometric_constant(self, lambdax, press):
        """ Psychrometric constant ratio of specific heat of moist air at
        a constant pressure to latent heat of vaporisation.

        References:
        -----------
        * Eqn 8 from FAO paper, Allen et al. 1998.

        Parameters:
        -----------
        lambdax : float
             latent heat of water vaporization [MJ kg-1]
        press : float
            average daytime pressure [kPa]

        Returns:
        --------
        gamma : float
            pyschrometric_constant [kPa degC-1]

        """
        return (self.cp * press) / (self.epsilon * lambdax)

    def calc_atmos_pressure(self):
        """ Pressure exerted by the weight of earth's atmosphere.

        References:
        -----------
        * Eqn 7 from FAO paper, Allen et al. 1998.

        Returns:
        --------
        press : float
            modelled average daytime pressure [kPa]

        """
        return (101.3 * ((293.0 - 0.0065 * self.zele_sea) / (293.0))**5.26)

    def calc_latent_heat_of_vapourisation(self, tavg):
        """ After Harrison (1963), should roughly = 2.45 MJ kg-1

        Returns:
        -----------
        lambdax : float
             latent heat of water vaporization [MJ kg-1]
        """
        return 2.501 - 0.002361 * tavg

    def calc_density_of_air(self, tavg):
        """ Found in lots of places but only reference I could find it in that
        wasn't an uncited equation is Dawes and Zhang (2011). No doubt there
        is a better reference

        Parameters:
        -----------
        tavg : float
            average daytime temperature [degC]

        Returns:
        --------
        density : float
            density of air [kg m-3]
        """
        return 1.292 - (0.00428 * tavg)

class Penman(PenmanMonteith):
    """
    Evaporation at the potential/equilibrium rate, where aerodynamic conductance
    is zero (i.e. winds are calm).

    References
    ----------
    * Monteith and Unsworth (1990) Principles of Environmental
      Physics, pg. 185-187.
    """

    def calc_evaporation(self, net_rad, tavg, press):
        """ Equilibrium evaporation

        Parameters:
        -----------
        net_rad : float
            net radiation [mj m-2 day-1]
        tavg : float
            daytime average temperature [degC]
        press : float
            average daytime pressure [kPa]

        Returns:
        --------
        soil_evap : float
            bare soil evaporation [mm day-1]

        """
        if press == None:
            press = self.calc_atmos_pressure()

        lambdax = self.calc_latent_heat_of_vapourisation(tavg)
        gamma = self.calc_pyschrometric_constant(lambdax, press)
        slope = self.calc_slope_of_saturation_vapour_pressure_curve(tavg)

        return ((slope / (slope + gamma)) * net_rad) / lambdax


class PriestleyTaylor(PenmanMonteith):

    """
    Calculate ET using Priestley Taylor, "potenial evaporation", i.e.
    simplified Penman method (radiation, temperature are the only inputs).
    Justification is that ET is generally determined by Rnet, rather than
    wind and air dryness.

    Key assumption is that the role of the soil heat flux is ignored at daily
    time scales. Not sure this holds...

    Penman-Monteith eqn aerodynamic term replaced by empirical multiplier, 1.26.
    Quoted range from literature for value is 1.2-1.3, although I have seen
    papers with lower values e.g. Viswanadham et al. 1991, Forest Ecology and
    Management, 38, 211-225.


    References:
    -----------
    * Priestley and Taylor (1972) On the assessment of surface heat flux and
      evaporation using large-scale parameters. Monthly Weather Review, 100,
      81-82.
    """

    def calc_evaporation(self, net_rad, tavg, press, pt_coeff=1.26):
        """
        Parameters:
        -----------
        net_rad : float
            net radiation [mj m-2 day-1]
        tavg : float
            daytime average temperature [degC]
        press : float
            average daytime pressure [kPa]
        pt_coeff : float, optional
            Priestley-Taylor coefficient

        Returns:
        --------
        transpiration : float
            transpiration [mm day-1]
        """
        lambdax = self.calc_latent_heat_of_vapourisation(tavg)
        gamma = self.calc_pyschrometric_constant(lambdax, press)
        slope = self.calc_slope_of_saturation_vapour_pressure_curve(tavg)

        return (pt_coeff / lambdax) * (slope / (slope + gamma)) * net_rad


def make_empty_list_of_zeros(size):
    """ create an empty (zero'd) list of a given size """
    return [0.0]*size 
