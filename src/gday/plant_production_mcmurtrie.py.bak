""" Original empirical G'DAY photosynthesis model """

__author__  = "Martin De Kauwe"
__version__ = "1.0 (23.02.2011)"
__email__   = "mdekauwe@gmail.com"

import sys
import math 
import constants as const
from utilities import float_eq, float_lt, float_gt

class PlantProdModel(object):
    """ 'simple' plant production model based on McMurtrie (1991).
        
    Calculate potential net photosynthesis based on radiation, efficiency with 
    which radiation is utilised (lue), light interception and leaf N concontrol. 
    
    References:
    ===========
    * McMurtrie, R. E. (2000) Relationship of forest productivity to 
        nutrient and carbon supply - a modelling analysistate. 
        Tree Physiology, 9, 87-99.
    
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
        
    def calculate_photosynthesis(self, day):
        """
        Parameters:
        -----------
        day : integer
            simulation day
        
        """
        option = self.control.model_number
    
        if option == 0:
            net_production = self.model_zero(day)
        elif option == 1:
            net_production = self.model_one(day)
        elif option == 2:
            net_production = self.model_two(day)
        elif option == 3:    
            net_production = self.model_three(day)
        elif option == 4:
            net_production = self.model_four(day)
        else:
            raise ValueError('Unknown assimilation model option %s\n' % 
                                                    self.control.model_number)
        
        self.fluxes.gpp = net_production / self.params.cue
        self.fluxes.npp = net_production
        self.fluxes.npp_gCm2 = (self.fluxes.npp * const.M2_AS_HA / 
                                    const.G_AS_TONNES)
        self.fluxes.gpp_gCm2 = self.fluxes.npp_gCm2 / self.params.cue
        
        
    def model_zero(self, day):
        """
        References:
        -----------
        * Kirschbaum et al 1994 pc & e 
        
        Parameters:
        -----------
        day : integer
            simulation day
        
        Returns: 
        --------
        npp : float
            net primary productivity
            
        """
        (sw_rad, ca, temp) = self.get_met_data(day)   
        if float_eq(sw_rad, 0.0):
            sw_rad = 0.000001
        
        gpp_max = ((self.params.cfracts * sw_rad / const.M2_AS_HA * 
                    self.params.epsilon * const.G_AS_TONNES) * temp_dep(temp) / 
                    temp_dep(14.0))
        
        lue = (1.21) * self.state.shootnc  / (0.0076 + self.state.shootnc)
        
        return lue * gpp_max * self.state.light_interception
        
    def model_one(self, day):       
        """ old g'day n:c and temperature factors
        
        References:
        -----------
        * mcmurtrie et al 1992 aust j bot
        
        Parameters:
        -----------
        day : integer
            simulation day
        
        Returns: 
        --------
        npp : float
            net primary productivity
            
        """
        (sw_rad, ca, temp) = self.get_met_data(day)
    
        if float_eq(sw_rad, 0.0):
            sw_rad = 0.000001
        rco2 = 1.632 * (ca - 60.9) / (ca + 121.8)
        gpp_max = ((self.params.cfracts * sw_rad / const.M2_AS_HA * 
                    self.params.epsilon * const.G_AS_TONNES) * rco2)
        
        # leaf n:c effect on photosynthetic efficiency (theta = 1).
        if float_gt(self.state.shootnc, self.params.n_crit):
            lue = 1.0
        else:
            rat = self.state.shootnc / self.params.n_crit
            if float_eq(self.params.ncpower, 1.0):
                lue = rat
            else:
                lue = math.pow(rat, self.params.ncpower)
        
        return lue * gpp_max * self.state.light_interception
        
    def model_two(self, day):
        """
        Parameters:
        -----------
        day : integer
            simulation day
        
        Returns: 
        --------
        npp : float
            net primary productivity
        
        """
        (sw_rad, ca, temp) = self.get_met_data(day)
        
        if float_eq(sw_rad, 0.0):
            sw_rad = 0.000001              
        rco2 = 1.82 * (ca - 50.0) / (ca + 197.0)
        gpp_max = sw_rad / const.M2_AS_HA * const.G_AS_TONNES
        lue = ((1.45 * self.state.ncontent - 1.29) / 
                    (self.state.ncontent + 1.65) * rco2)
        
        return lue * gpp_max * self.state.light_interception
        
    def model_three(self, day):
        """
        Parameters:
        -----------
        day : integer
            simulation day
        
        Returns: 
        --------
        npp : float
            net primary productivity
        
        """
        (sw_rad, ca, temp) = self.get_met_data(day)
        
        if float_eq(sw_rad, 0.0):
            sw_rad = 0.000001           
        rco2 = 1.82 * (ca - 50.0) / (ca + 197.0)
        lue = (1.45 * self.state.ncontent) / (self.state.ncontent + 2.21) * rco2
        gpp_max = sw_rad / const.M2_AS_HA * const.G_AS_TONNES
        
        
        return lue * gpp_max * self.state.light_interception
           
            
    def model_four(self, day):
        """
        Parameters:
        -----------
        day : integer
            simulation day
        
        Returns: 
        --------
        npp : float
            net primary productivity
        
        """
        (sw_rad, ca, temp) = self.get_met_data(day)
        
        rco2 = 2.34 * (ca - 50.0) / (ca + 354.0)
        
        gpp_max = sw_rad / const.M2_AS_HA * const.G_AS_TONNES
        lue = (self.state.ncontent + 0.2) / (self.state.ncontent + 2.67) * rco2
        
        return lue * gpp_max * self.state.light_interception
    
    def get_met_data(self, day):
        """ Grab the days met data out of the structure and return day values. 
        
        Parameters:
        ----------
        day : int
            project day. 
        
        Returns:
        -------
        sw_rad : float
            sw down radiation [mj m-2 d-1]
        ca : float
            atmospheric co2, depending on flag set in param file this will be
            ambient or elevated. [umol mol-1]
        
        """
        temp = self.met_data['tair'][day]   
        sw_rad = self.met_data['sw_rad'][day] 
        if self.control.co2_conc == 0:
            ca = self.met_data['amb_co2'][day]
        elif self.control.co2_conc == 1:
            ca = self.met_data['ele_co2'][day]
        
        return sw_rad, ca, temp
        
def temp_dep(x):
    """Temperature dependency func
    
    Parameters:
    -----------
    x : float
        temperature
    
    Returns:
        temp dependence
    
    """
    return 1.0 / (1.0 + math.exp(1.693 - 0.1047 * x))            



def cmdline_parser():
    """ Parse the command line for user options """
    from optparse import OptionParser
    
    desc = """G'DAY photosynthesis submodel"""
    clp = OptionParser("Usage: %prog [options] filename", description = desc)
    options, args = clp.parse_args()
    return options, args
     
    
def main():
    
    # sweep the cmd line
    options, args = cmdline_parser()
    
    from file_parser import initialise_model_data
    # pylint: disable=C0103
    # pylint: disable=C0324
    # pylint: disable=C0103
    
    fname = "gday"
    fdir = "/Users/mdekauwe/src/python/GDAY_model/params"
    (adj_control, adj_params, 
        adj_state, adj_files, 
        adj_fluxes, met_data) = initialise_model_data(fname, default_dir=fdir)
                                            
    # figure out photosynthesis
    P = PlantProdModel(adj_control, adj_params, adj_state, adj_fluxes, met_data)
    
    adj_state.lai = (adj_params.slainit * const.M2_AS_HA / 
                            const.KG_AS_TONNES / adj_params.cfracts * 
                            adj_state.shoot)
    
    # Specific LAI (m2 onesided/kg DW)
    adj_state.sla = adj_params.slainit
    
    
    
    adj_control.model_number = 3
    
    num_days = len(met_data['doy'])
    for i in xrange(num_days):
        
        if float_lt(adj_state.lai, adj_params.lai_cover):
            gcover = adj_state.lai / adj_params.lai_cover
        else:
            gcover = 1.0
        
        adj_state.fapar = ((1.0 - math.exp(-adj_params.kext * adj_state.lai / 
                    gcover)) * gcover) 
        
        adj_state.shootnc = adj_state.shootn / adj_state.shoot 
        
        P.run_sim(i)
        
       
        print adj_fluxes.gpp / const.HA_AS_M2 * const.TONNES_AS_G
            
        #print adj_state.lai    
            
    
        # this is done in derive so do here
        # Specific LAI (m2 onesided/kg DW)
        adj_state.sla = (adj_state.lai / const.M2_AS_HA * 
                            const.KG_AS_TONNES * 
                            adj_params.cfracts / adj_state.shoot)
   
    
    return 
    
    
    
if __name__ == "__main__":
    sys.exit(main())









