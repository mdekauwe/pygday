# -*- coding: UTF-8 -*-

import math

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
        self.delta_store = 0.0
        
    def calculate_water_balance(self, day, daylen):
        """ Calculate water balance

        Parameters:
        ----------
        day : int
            project day.
        daylen : float
            length of day in hours.

        """
        # met forcing
        (tavg, rain, sw_rad, vpd, wind, net_rad, ca, 
         press, amb_co2, net_rad_am, net_rad_pm, vpd_am,
         vpd_pm, temp_am, temp_pm, 
         wind_am, wind_pm) = self.get_met_data(day, daylen)
        
        # calculate water fluxes
        if self.control.trans_model == 0:
            # transpiration calculated from WUE...
            self.calc_wue(vpd, ca, amb_co2)
            self.calc_transpiration()
        elif self.control.trans_model == 1:
            
            #self.calc_transpiration_penmon(vpd, net_rad, tavg, wind,
            #                                    ca, daylen, press)
            self.calc_transpiration_penmon_am_pm(net_rad_am, net_rad_pm, wind, ca, daylen, 
                                                 press, vpd_am, vpd_pm, 
                                                 temp_am, temp_pm, 
                                                 wind_am, wind_pm)
            self.calc_wue(vpd, ca, amb_co2)
        elif self.control.trans_model == 2:
            self.calc_transpiration_priestay(net_rad, tavg, press)
            self.calc_wue(vpd, ca, amb_co2)

        self.calc_infiltration(rain)
        self.fluxes.soil_evap = self.calc_soil_evaporation(tavg, net_rad,
                                                            press)
        self.fluxes.et = (self.fluxes.transpiration + self.fluxes.soil_evap +
                            self.fluxes.interception)

        self.fluxes.runoff = self.update_water_storage()
       
        #water_balance = rain - self.fluxes.et- self.fluxes.runoff - self.delta_store
        
        
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
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated. [umol mol-1]
        wind : float
            average daily wind speed [mm s-1]
        press : float
            average daytime pressure [kPa]

        """
        tavg = self.met_data['tair'][day]
        temp_am, temp_pm = self.met_data['tam'][day], self.met_data['tpm'][day]
        sw_rad_am, sw_rad_pm = self.met_data['sw_rad_am'][day], self.met_data['sw_rad_pm'][day]
        rain = self.met_data['rain'][day]
        sw_rad = self.met_data['sw_rad'][day]
        vpd = self.met_data['vpd_avg'][day] # daytime average
        wind = self.met_data['wind'][day]
        net_rad = self.calc_radiation(tavg, sw_rad, daylen)
        net_rad_am = self.calc_radiation(temp_am, sw_rad_am, daylen/2.0)
        net_rad_pm = self.calc_radiation(temp_pm, sw_rad_pm, daylen/2.0)
        vpd_am, vpd_pm =self.met_data['vpd_am'][day], self.met_data['vpd_pm'][day]
        temp_am, temp_pm = self.met_data['tam'][day], self.met_data['tpm'][day]
        wind_am, wind_pm = self.met_data['wind_am'][day], self.met_data['wind_pm'][day]
        
        if self.control.co2_conc == 0:
            #ca = 385.0
            ca = self.met_data['amb_co2'][day]
        elif self.control.co2_conc == 1:
            #ca = 550.0
            ca = self.met_data['ele_co2'][day]
        amb_co2 = self.met_data['amb_co2'][day]
        
        
        if ('atmos_press' in self.met_data and not
            self.met_data['atmos_press'] is None):
            press = self.met_data['atmos_press'][day]
        else:
            press = None # use method below to calculate pressure

        return (tavg, rain, sw_rad, vpd, wind, net_rad, ca, press, amb_co2,
                net_rad_am, net_rad_pm, vpd_am, vpd_pm, temp_am, temp_pm, 
                wind_am, wind_pm)

    def calc_wue(self, vpd, ca, amb_co2):
        """water use efficiency

        Not sure of units conversions here, have to ask BM

        Parameters:
        -----------
        vpd : float
            average daily vpd [kPa]
        ca : float
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated. [umol mol-1]

        """
        if self.control.wue_model == 0:
            # Gday original implementation
            # (gC / kg H20)
            if float_gt(vpd, 0.0):
                # WUE Power law dependence on co2, Pepper et al 2005.
                
                co2_ratio = (ca / amb_co2)
                co2_adjustment = co2_ratio**self.params.co2_effect_on_wue

                # wue inversely proportional to daily mean vpd
                self.fluxes.wue = self.params.wue0 * co2_adjustment / vpd
            else:
                self.fluxes.wue = 0.0
        elif self.control.wue_model == 1 and self.control.assim_model == 7:
            conv = const.MOL_C_TO_GRAMS_C / const.MOL_WATER_TO_GRAMS_WATER
            self.fluxes.wue = (conv * 1000.0 * (ca * const.UMOL_TO_MOL *
                                (1.0 - self.fluxes.cica_avg) /
                                (1.6 * vpd / 101.0)))
                        #if self.fluxes.wue > 20.0: self.fluxes.wue = 20.0    # FIX THIS!!!
            # what is this? ask BM

        elif self.control.wue_model == 2:
            self.fluxes.wue = (self.params.wue0 * 0.27273 / vpd *
                                ca / amb_co2)
        elif self.control.wue_model == 3 :
            if float_eq(self.fluxes.transpiration, 0.0):
                self.fluxes.wue = 0.0
            else:
                self.fluxes.wue = (self.fluxes.gpp_gCm2 /
                                    self.fluxes.transpiration)
        else:
            raise AttributeError('Unknown WUE calculation option')

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
        
        self.fluxes.erain = max(0.0, rain * self.params.rfmult -
                                self.state.lai * self.params.wetloss)
        self.fluxes.interception = rain * self.params.rfmult - self.fluxes.erain
        

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

    def calc_transpiration_penmon(self, vpd, net_rad, tavg, wind, ca,
                                        daylen, press):
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
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated. [umol mol-1]
        daylen : float
            daylength in hours
        press : float
            average daytime pressure [kPa]

        """
        
        P = PenmanMonteith(canht=self.params.canht, dz0v_dh=self.params.dz0v_dh,
                            displace_ratio=self.params.displace_ratio)
        
        gs = self.calc_stomatal_conductance(vpd, ca, daylen, self.fluxes.gpp_gCm2,
                                            press, tavg)
        
        transp = P.calc_evaporation(vpd, wind, gs, net_rad, tavg, press)
        
        tconv = (60.0 * 60.0 * daylen) # seconds to day
        self.fluxes.transpiration = transp * tconv
       
       
        
    def calc_transpiration_penmon_am_pm(self, net_rad_am, net_rad_pm, wind, ca, 
                                        daylen, press, vpd_am, vpd_pm, temp_am, 
                                        temp_pm, wind_am, wind_pm):
        """ Calculate canopy transpiration using the Penman-Monteith equation
        using am and pm data.
        units mm/day

        Parameters:
        -----------
        vpd : float
            average daily vpd [kPa]
        net_rad_am : float
            net radiation [mj m-2 s-1] (morning)
        net_rad_pm : float
            net radiation [mj m-2 s-1] (afternoon)
        tavg : float
            average daytime temp [degC]
        wind : float
            average daily wind speed [m s-1]
        ca : float
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated. [umol mol-1]
        daylen : float
            daylength in hours
        press : float
            average daytime pressure [kPa]

        """
        P = PenmanMonteith(canht=self.params.canht, dz0v_dh=self.params.dz0v_dh,
                            displace_ratio=self.params.displace_ratio)
        
        gs_am = self.calc_stomatal_conductance(vpd_am, ca, daylen/2.0, 
                                                self.fluxes.gpp_am_gCm2, 
                                                press, temp_am)
        gs_am_mol_m2_sec = self.fluxes.gs_mol_m2_sec
        
        gs_pm = self.calc_stomatal_conductance(vpd_pm, ca, daylen/2.0, 
                                                self.fluxes.gpp_pm_gCm2,
                                                press, temp_am)
        gs_pm_mol_m2_sec = self.fluxes.gs_mol_m2_sec
        gs = (gs_am + gs_pm) / 2.0
        self.fluxes.gs_mol_m2_sec = (gs_am_mol_m2_sec + gs_pm_mol_m2_sec) / 2.0
        
        
        transp_am = P.calc_evaporation(vpd_am, wind_am, gs_am, net_rad_am, 
                                        temp_am, press)
        transp_pm = P.calc_evaporation(vpd_pm, wind_pm, gs_pm, net_rad_pm, 
                                        temp_pm, press)
        
        ga = P.calc_atmos_boundary_layer_conductance(wind)
        self.fluxes.ga_mol_m2_sec = ga / const.CONV_CONDUCT
        
        tconv = (60.0 * 60.0 * daylen) # seconds to day
        self.fluxes.transpiration = (transp_am + transp_pm) / 2.0 * tconv
        
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
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated. [umol mol-1]
        daylen : float
            daylength in hours

        Returns:
        --------
        gs : float
            stomatal conductance [m s-1]
        """
        # convert conductance to water vapour units
        g1_c = self.params.g1 / const.RATIO_DIFF_H2O_TO_CO2
       
        # time unit conversion day-1 -> seconds-1
        tconv =  1.0 / (60.0 * 60.0 * daylen)
        gpp_umol_m2_sec = (gpp * const.GRAMS_C_TO_MOL_C * const.MOL_TO_UMOL * tconv)
        
        arg1 = 1.6 * (1.0 + (g1_c * self.state.wtfac_root) / math.sqrt(vpd))
        arg2 = gpp_umol_m2_sec / ca # umol mol-1
        gs_mol_m2_sec = arg1 * arg2 * const.RATIO_DIFF_H2O_TO_CO2
        
        
        self.fluxes.gs_mol_m2_sec = gs_mol_m2_sec
        
        # convert to mm s-1 and then to m s-1
        #return (gs_mol_m2_sec * const.MOL_TO_MILLIMOLES * const.CONV_CONDUCT * 
        #        const.MM_TO_M)
        
        # See Jones, 1992, appendix
        tk = temp + const.DEG_TO_KELVIN
        conv = const.MM_TO_M / (press / (const.RGAS * tk))
        
        # convert to mm s-1 to m s-1
        return (gs_mol_m2_sec * conv) 
   
    
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

    def calc_soil_evaporation(self, tavg, net_rad, press):
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
        
        soil_evap *= math.exp(-0.398 * self.state.lai)
        
        # if the available soil moisture is low the soil evaporation needs to
        # be reduced as well
        soil_evap *= self.state.wtfac_tsoil
        
        return soil_evap

    def update_water_storage(self):
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
        # Total root zone
        
        prev = self.state.pawater_root
        self.state.pawater_root += (self.fluxes.erain -
                                    self.fluxes.transpiration -
                                    self.fluxes.soil_evap)
       
        if self.state.pawater_root > self.params.wcapac_root:
            runoff = self.state.pawater_root - self.params.wcapac_root 
            
        else:
            runoff = 0.0
            
        self.state.pawater_root = clip(self.state.pawater_root, min=0.0,
                                        max=self.params.wcapac_root)
        
        self.delta_store = self.state.pawater_root - prev
       
        # Total soil layer
        self.state.pawater_tsoil += (self.fluxes.erain -
                                     self.fluxes.transpiration *
                                     self.params.fractup_soil -
                                     self.fluxes.soil_evap)

        self.state.pawater_tsoil = clip(self.state.pawater_tsoil, min=0.0,
                                        max=self.params.wcapac_topsoil) 
        
        return runoff

    def calculate_soil_water_fac(self):
        """ Estimate a relative water availability factor [0..1]

        A drying soil results in physiological stress that can induce stomatal
        closure and reduce transpiration. Further N mineralisation depends on 
        top soil moisture.

        References:
        -----------
        * Pepper et al. (2008) Functional Change Biology, 35, 493-508

        But similarly see:
        * van Genuchten (1981) Soil Sci. Soc. Am. J, 44, 892--898.
        * Wang and Leuning (1998) Ag Forest Met, 91, 89-111.

        Returns:
        --------
        wtfac_tsoil : float
            water availability factor for the top soil [0,1]
        wtfac_root : float
            water availability factor for the root zone [0,1]    
        """
        # turn into fraction...
        smc_root = self.state.pawater_root / self.params.wcapac_root
        smc_topsoil = self.state.pawater_tsoil / self.params.wcapac_topsoil

        # Calculate a soil moisture availability factor, used to adjust
        # ci/ca ratio in the face of limited water supply.
        arg = self.params.fwpmax - self.params.fwpmin
        wtfac_tsoil = (smc_topsoil - self.params.fwpmin) / arg
        wtfac_root = (smc_root - self.params.fwpmin) / arg
        
        #return 1.0, 1.0
        return (clip(wtfac_tsoil, min=0.0, max=1.0), 
                clip(wtfac_root, min=0.0, max=1.0))
        
        
        


class WaterLimitedNPP(object):
    """ Adjust carbon uptake for water limitations

    Water limited growth depends on the maximum rate of which water can
    be extracted by the roots
    """
    def __init__(self, control, params, state, fluxes):
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

        """
        self.params = params
        self.fluxes = fluxes
        self.control = control
        self.state = state

    def adjust_cproduction(self, option):
        """ select model?

        It seems hybrid is the same as biomass, so check this and remove. I
        have set it up so that a call to hybrid just calls biomass

        Parameters:
        -----------
        option : integer
            model option

        """
        if option == 1 or option == 3:
            self.monteith_rescap_model()
        elif option == 2 and float_gt(self.params.fwpmax, self.params.fwpmin):
            self.biomass_model()
        elif option == 3 and float_gt(self.params.fwpmax, self.params.fwpmin):
            self.hybrid_model() # same as calling biomass model!
        else:
            err_msg = "Unknown water bal model (try 1-3): %s\n" % option
            raise RuntimeError, err_msg

    def monteith_rescap_model(self):
        """Growth on a given day will either be

        light or water limited and this depends on the plants ability to
        capture resources (leaves to intercept light + roots to supply
        water)

        """
        if self.control.fixroot:
            rootz = self.params.fixed_rootc
        else:
            rootz = self.state.root

        growth_water_lim = (self.fluxes.wue * self.params.extraction /
                                100.0 * rootz * self.state.pawater_root)

        # growth is the lesser of two limiting (water and light) rates
        # light limited growth is the already calculated npp
        self.fluxes.npp = min(self.fluxes.npp, growth_water_lim)

    def biomass_model(self):
        """ Water limit on carbon production: Biomass model """

        # water limitation factor for NPP, which depends on total plant water.
        # Used to limit growth (i.e. water stress)
        self.fluxes.npp *= self.state.wtfac_root

    def hybrid_model(self):
        """Rescap/biomass hybrid water balance model

        water-limited transpirationn is proportional to npp*wtfac

        """
        # the original code has this setup identical to the biomass model, are
        # they meant to be different? For the moment I will just call biomass
        # model from here but really this should be remvoed?
        self.biomass_model()


class PenmanMonteith(object):

    """ Water loss from a canopy (ET), representing surface as a big "leaf".
    The resistance to vapour transfer from the canopy to the atmosphere is
    determined by the aerodynamic resistance (1/ga) and the canopy resistance
    (1/gc). If the surface is wet then there is a further water vapour flux
    from the soil/surface (calculated elsewhere!).

    Assumes inputs are measured at a 2m screen height, though you can change
    this with the the class initialisation. I have ignored the soil heat
    flux (G) assuming it balances to zero over the course of the day.

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
                    canht=20.0, dz0v_dh=0.1, displace_ratio=0.67):

        """
        Parameters:
        -----------
        cp : float
            specific heat of dry air [MJ kg-1 degC-1]
        vk : float
            von karman's constant [unitless]
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

        # ratio of the roughness length for heat to the length for momentum
        self.z0h_z0m = 0.1

    def calc_evaporation(self, vpd, wind, gs, net_rad, tavg, press):

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
        
        ga = self.calc_atmos_boundary_layer_conductance(wind)
       
        # our model is a big leaf, so canopy conductance, gc = gs
        gc = gs
        
        if float_gt(gc, 0.0):
            # decoupling coefficent, Jarvis and McNaughton, 1986
            e = slope / gamma # chg of latent heat relative to sensible heat of air
            omega = (e + 1.0) / (e + 1.0 + (ga / gc))
            
            arg1 = ((slope * net_rad ) + (rho * self.cp * vpd * ga))
            arg2 = slope + gamma * (1.0 + ga / gc)
            et = (arg1 / arg2) / lambdax
        else:
            et = 0.0
            omega = 0.0
        
        return et

    def calc_atmos_boundary_layer_conductance(self, wind):
        """ atmospheric boundary layer conductance, i.e. 1/ra

        Transfer of heat/water vapour from evaporating surface into air
        above the canopy is determined by aerodynamic conductance. Key
        assumption is that roughness length for momentum and for heat are
        identical.

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
        #wind = self.adj_wind_speed_2_screen(wind, canht)
        # roughness length [m]
        z0 = self.dz0v_dh * self.canht
        
        # zero plan displacement height [m]
        d = self.displace_ratio * self.canht

        arg1 = self.vk**2 * wind
        arg2 = (math.log((self.canht - d) / z0))**2
        
        return arg1 / arg2
    
    
    
    def adj_wind_speed_2_screen(self, wind):
        """ Standard PM expects wind speed at a 2m height, need to adjust
        observations for the height difference

        References:
        -----------
        * Bos et al (2009) Water Requirements for Irrigation and the
          Environment, pg. 29
        * http://www.apesimulator.it/help/models/evapotranspiration/\
          Aerodynamic_resistance.html

        Parameters:
        -----------
        wind : float
            wind speed [m s-1]

        Returns:
        --------
        wind_speed : float
            adjusted wind speed for height [m s-1]
        """
        wind *= 4.87 / (math.log(67.8 * self.canht - 5.42))
        return wind

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
        arg1 = 4098.0 * (0.6108 * math.exp((17.27 * tavg) / t))
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


