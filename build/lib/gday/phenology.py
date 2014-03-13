from math import floor, exp
from utilities import day_length


class Phenology(object):
    """ Phenology scheme based on Botta calibration against EO data. The 
    estimated leaf out/off dates are used to define the allocation of C&N in the
    deciduous model as a function of soil+air temperature and daylength.
    
    References:
    -----------
    Botta, A. et al. (2000) GCB, 6, 709-725.
    """
    
    def __init__(self, fluxes, state, control, previous_ncd=None, pa=-68., 
                 pb=638., pc=-0.01, store_transfer_len=None):
        """
        Parameters:
        ----------
        pa : float
            (days) Leaf flush params following Botta.
        pb : float
            (days) Leaf flush params following Botta.
        pc : float
            (1/days) Leaf flush params following Botta.
        previous_ncd : int
            Previous years number of chilling days
        """
        self.pa = pa 
        self.pb = pb
        self.pc = pc
        self.Tbase = 5.0 # degC
        self.last_yrs_accumulated_ncd = previous_ncd
        self.accumulated_ncd = 0.0
        self.accum_gdd = 0.0
        self.project_day = 0    
        self.leaf_on = 0.0
        self.leaf_off = 0.0
        self.leaf_on_found = False
        self.leaf_off_found = False
        self.drop_leaves = False
        self.fluxes = fluxes
        self.state = state
        self.control = control
        self.store_transfer_len = store_transfer_len
        self.growing_seas_len = None
        
    def calculate_phenology_flows(self, daylen, met_data, yr_days, 
                                  project_day):
        self.project_day = project_day
        self.calculate_leafon_off(daylen, met_data, yr_days)
        self.calculate_days_left_in_growing_season(yr_days)
        self.calculate_growing_season_fluxes()
        
    def calc_gdd(self, Tavg):
        """ calculate the number of growing degree days, hypothesis is that
        leaves appear after a threshold has been passed.
        """
        return max(0.0, Tavg - self.Tbase)
  
    def gdd_chill_thresh(self, ncd):
        """ Leaf out has a chilling requirement, num chill days reduces the GDD
        requirement 
        """
        return self.pa + self.pb * exp(self.pc * ncd)
    
    def calc_ncd(self, Tmean):
        """ Calculate the number of chilling days from fixed dates (1 Nov), 
        following Murray et al. 1989, same as Botta does. """
        if Tmean < 5.0:
            return 1.0
        else:
            return 0.0

    def leaf_drop(self, daylen, Tsoil, Tsoil_next_3days):
        """ Thresholds to drop leaves come from White et al.
        Note 655 minutes = 10.916 hrs.
        
        References:
        -----------
        White, M. A. et al. (1997) GBC, 11, 217*234.
        """
        if (daylen <= 10.9166667 and Tsoil <= 11.15) or Tsoil_next_3days < 2.0:
            return True
        else:
            return False
        
    def ini_phen_calcs(self):
        self.accumulated_ncd = 0.0
        self.accum_gdd = 0.0
        self.leaf_on = 0.0
        self.leaf_off = 0.0
        self.leaf_on_found = False
        self.leaf_off_found = False
        self.drop_leaves = False
        
    def calculate_leafon_off(self, daylen, met_data, yr_days):
        self.ini_phen_calcs()
        gdd_thresh = self.gdd_chill_thresh(self.last_yrs_accumulated_ncd)
        
        if yr_days == 366:
            nov_doy = 306
        else:
            nov_doy = 305
        for d in xrange(yr_days):
            
            Tmean = met_data['tair'][self.project_day]
            Tsoil = met_data['tsoil'][self.project_day]
            
            if self.project_day < 362:
                Tsoil_next_3days = ((met_data['tsoil'][self.project_day] +
                                     met_data['tsoil'][self.project_day+1] + 
                                     met_data['tsoil'][self.project_day+2])/3.0)
            else:
                # i.e. end of year, didn't find this so have no effect
                Tsoil_next_3days = 999.9 
            
            # Sum the daily mean air temperature above 5degC starting on Jan 1
            # (july 1 for south hemisphere) - note code needs to be changed for
            # SH!
            self.accum_gdd += self.calc_gdd(Tmean)
            
            if self.leaf_on_found == False and self.accum_gdd >= gdd_thresh:
                self.leaf_on = d
                self.leaf_on_found = True
                
            if self.leaf_off_found == False and (self.accum_gdd >= gdd_thresh):
                # I am prescribing that no leaves can fall off before doy=180
                # Had issue with KSCO simulations where the photoperiod was
                # less than the threshold very soon after leaf out.
                if d > 180:
                    self.drop_leaves = self.leaf_drop(daylen[d], Tsoil, 
                                                      Tsoil_next_3days)   
                    if self.drop_leaves:
                        self.leaf_off_found = True
                        self.leaf_off = d
                    
            # Calculated NCD from fixed date following Murray et al 1989.
            if d+1 >= nov_doy:
                self.accumulated_ncd += self.calc_ncd(Tmean) 
                
            self.project_day += 1
    
        self.last_yrs_accumulated_ncd = self.accumulated_ncd
        
        # Length of time taken for new growth from storage to be allocated.
        # This is either some site-specific calibration or the midpoint of the
        # length of the growing season. The litterfall takes place over an 
        # identical period. Dividing by a larger number would increase the
        # rate the C&N is allocated.
        self.growing_seas_len = self.leaf_off - self.leaf_on
        if self.store_transfer_len == None:
            self.len_groloss = floor((self.growing_seas_len) / 2.0)
        else:
            self.len_groloss = self.store_transfer_len
        
        
        
        
    def calculate_days_left_in_growing_season(self, yr_days):
        """ Calculate 2 arrays to signify the days left of growing period
        and days left before all the leaves fall off. In both cases these will
        be 2 lists, with 0.0 outside of the growing period and a series of 
        numbers e.g. day 46 to 0 for growing_days. 0.5 is subtracted from the
        doy to get round the issue of approximating an integral with discrete
        time steps -> trapezoidal type solution
        """
        
        self.state.remaining_days = [] 
        self.state.growing_days = [] 
        self.state.leaf_out_days= [] 
        for doy in xrange(1, yr_days+1):
            if doy > self.leaf_off - self.len_groloss and doy <= self.leaf_off:
                self.state.remaining_days.append((doy - 0.5) - 
                                                  self.leaf_off + 
                                                  self.len_groloss)
            else:
                self.state.remaining_days.append(0.0)
            
            if doy > self.leaf_on and doy <= self.len_groloss+self.leaf_on:
                self.state.growing_days.append(self.len_groloss + 
                                               self.leaf_on - (doy - 0.5))             
            else:
                self.state.growing_days.append(0.0)
            
            if doy > self.leaf_on and doy < self.leaf_off:
                self.state.leaf_out_days.append(1.0)             
            else:
                self.state.leaf_out_days.append(0.0)
        
    def calculate_growing_season_fluxes(self):
        
        # C allocation rates across growing season
        self.fluxes.lrate = (2.0 * self.state.c_to_alloc_shoot / 
                             self.len_groloss**2)
        self.fluxes.wrate = (2.0 * self.state.c_to_alloc_stem / 
                             self.len_groloss**2)   
        self.fluxes.brate = (2.0 * self.state.c_to_alloc_branch / 
                             self.len_groloss**2)   
        
        
        # N allocation rates across growing season
        self.fluxes.lnrate = (2.0 * self.state.n_to_alloc_shoot / 
                              self.len_groloss**2)
        self.fluxes.bnrate = (2.0 * self.state.n_to_alloc_branch / 
                              self.len_groloss**2)
        self.fluxes.wnimrate = (2.0 * self.state.n_to_alloc_stemimm / 
                                self.len_groloss**2)
        self.fluxes.wnmobrate = (2.0 * self.state.n_to_alloc_stemmob / 
                                 self.len_groloss**2)
       
        