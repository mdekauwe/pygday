import math 
from utilities import day_length


class Phenology(object):
    """ Prognostic phenology model calibrated from EO data
    
    References:
    -----------
    Botta, A. et al. (2000) GCB, 6, 709-725.
    """
    
    def __init__(self, fluxes, state, previous_ncd, pa=-68., pb=638., 
                 pc=0.01):
        """
        Parameters:
        ----------
        pa : float
            (days) Leaf flush params following Botta.
        pb : float
            (days) Leaf flush params following Botta.
        pc : float
            (1/days) Leaf flush params following Botta.
        """
        self.pa = pa 
        self.pb = pb
        self.pc = pc
        self.Tbase = 5.0 # degC
        self.last_yrs_accumulated_ncd = previous_ncd
        self.accumulated_ncd = 0.0
        self.accumulated_gdd = 0.0
        self.project_day = 0    
        self.leaf_on = 0.0
        self.leaf_off = 0.0
        self.leaf_on_found = False
        self.leaf_off_found = False
        self.drop_leaves = False
        self.fluxes = fluxes
        self.state = state
        
    def calculate_phenology_flows(self, daylen, met_data, yr_days, 
                                  project_day):
        self.project_day = project_day
        self.calculate_leafon_off(daylen, met_data, yr_days)
        self.calculate_days_left_in_growing_season(yr_days)
        self.calculate_growing_season_fluxes()
        
    def calc_gdd(self, Tavg):
        return max(0.0, Tavg - self.Tbase)
  
    def gdd_chill_thresh(self, ncd):
        """ Leaf out has a chilling requirement, num chill days reduces the GDD
        requirement 
        """
        return self.pa + self.pb * math.exp(self.pc * ncd)
    
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
        if (daylen < 10.9166667 and Tsoil < 11.15) or Tsoil_next_3days < 2.0:
            return True
        else:
            return False
    
    def ini_phen_calcs(self):
        self.accumulated_ncd = 0.0
        self.accumulated_gdd = 0.0
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
                Tsoil_next_3days = met_data['tsoil'][self.project_day] + \
                                    met_data['tsoil'][self.project_day+1] + \
                                    met_data['tsoil'][self.project_day+2] / 3.0
            else:
                # i.e. end of year, didn't find this so have no effect
                Tsoil_next_3days = 999.9 
            self.accumulated_gdd += self.calc_gdd(Tmean)
            
            if self.leaf_on_found == False and self.accumulated_gdd >= gdd_thresh:
                self.leaf_on = d
                self.leaf_on_found = True
            if self.leaf_off_found == False and (self.accumulated_gdd >= gdd_thresh):
                
                self.drop_leaves = self.leaf_drop(daylen[d], Tsoil, Tsoil_next_3days)   
                if self.drop_leaves:
                    self.leaf_off_found = True
                    self.leaf_off = d
                
            # Calculated NCD from fixed date following Murray et al 1989.
            if d+1 >= nov_doy:
                self.accumulated_ncd += self.calc_ncd(Tmean) 
                
            self.project_day += 1
    
        self.last_yrs_accumulated_ncd = self.accumulated_ncd
        self.mid_point = math.floor((self.leaf_off - self.leaf_on) / 4.0)
        
    def calculate_days_left_in_growing_season(self, yr_days):
        self.state.remaining_days = [] 
        self.state.growing_days = [] 
        self.state.leaf_out_days= [] 
        for day in xrange(1, yr_days+1):
            if day>=self.leaf_off-self.mid_point and day<=self.leaf_off:
                self.state.remaining_days.append(day-self.leaf_off+self.mid_point)
            else:
                self.state.remaining_days.append(0.0)
            
            if day>=self.leaf_on and day<self.mid_point+self.leaf_on:
                self.state.growing_days.append(self.mid_point+self.leaf_on-day)             
            elif day==self.mid_point+self.leaf_on:
                self.state.growing_days.append(-999.9)
            else:
                self.state.growing_days.append(0.0)
            
            if day>self.leaf_on and day<self.leaf_off:
                self.state.leaf_out_days.append(1.0)             
            else:
                self.state.leaf_out_days.append(0.0)
        
    def calculate_growing_season_fluxes(self):
        
        # C allocation rates across growing season
        self.fluxes.lrate = 2.0 * self.state.c_to_alloc_shoot / self.mid_point**2
        self.fluxes.wrate = 2.0 * self.state.c_to_alloc_stem / self.mid_point**2   
            
        # N allocation rates across growing season
        self.fluxes.lnrate = 2.0 * self.state.n_to_alloc_shoot / self.mid_point**2 
        
        len_gs = len([i for i in self.state.growing_days if i > 0])
        totcf = sum([self.fluxes.lrate*i for i in xrange(len_gs+1)])
        totnf = sum([self.fluxes.lnrate*i for i in xrange(len_gs+1)])
        
        excess_alloc = totcf - self.state.c_to_alloc_shoot
        excess_allocn = totnf - self.state.n_to_alloc_shoot
        
        self.fluxes.lrate = 2.0 * (self.state.c_to_alloc_shoot-excess_alloc) / self.mid_point**2
        self.fluxes.lnrate = 2.0 * (self.state.n_to_alloc_shoot-excess_allocn) / self.mid_point**2 
        totcf = sum([self.fluxes.lrate*i for i in xrange(len_gs+1)])
        totnf = sum([self.fluxes.lnrate*i for i in xrange(len_gs+1)])
        self.state.left_over_c = self.state.c_to_alloc_shoot-totcf
        self.state.left_over_n = self.state.n_to_alloc_shoot-totnf
        #print self.state.c_to_alloc_shoot,totcf, self.state.left_over_c,totcf+self.state.left_over_c
        #import sys; sys.exit()
        
      