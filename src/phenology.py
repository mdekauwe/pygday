from math import floor, exp
from utilities import day_length


class Phenology(object):
    """ 
    There are two phenology schemes currently implemented, one which should 
    generally be applicable for deciduous broadleaf forests and one for 
    short-lived grasses.
    
    The tree scheme is based on Botta et al. who calibrated their model against
    EO data. The grasses scheme is really a combination of approaches, I've 
    tried to keep it is as simple as possible...user beware :)
    
    There are some key issues:
    - leaf drop for trees won't work where the daylength isn't applicable, i.e
      the tropics would be my guess.
    - the grass phenology drop takes no account of available water. It of course
      should, but this would need a complete re-write of the logic. Perhaps 
      someone else can do that.
      
    One potential thing to look at if you ever get bored is:
      
    - Look at Caldararu et al 2014, Phenology as a strategy for carbon 
      optimality: a global model. Seems interesting at a quick glance.
      
    The distribution of C&N is pre-calculated here using a ramping function
    based on leaf out/off dates.
    
    Finally, no account has been taken for the southern hemisphere! This won't
    work there.
    
    References:
    -----------
    * Botta, A. et al. (2000) GCB, 6, 709-725.
    * Foley, J. et al. (1996) GBC, 10, 603-628.
    * Krinner, G. et al. (2005) GBC, 19, GB1015
    * White, M. A. et al. (1997) GBC, 11, 217-234.
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
        
        - Dependance on daylength means that this is only valid outside of the
          tropics. 
        
        References:
        -----------
        White, M. A. et al. (1997) GBC, 11, 217-234.
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
        
    def calc_ini_grass_pheno_stuff(self, met_data, yr_days):
        """ Series of constraints based on temp and precip need to be 
        pre-calculated for grasses to determine leaf on/off
        """
        
        # Save this as we need to loop over the data once to pre-calculate 
        # everything 
        project_day_save = self.project_day 
        
        
        tmax_ann = 0.0
        tmin_ann = 70.0
        tavg_ann = 0.0
        ppt_sum = 0.0
        for d in xrange(yr_days):
            tair = met_data['tair'][self.project_day]
            tam = met_data['tam'][self.project_day]
            tpm = met_data['tpm'][self.project_day]
            ppt_sum += met_data['rain'][self.project_day]
            
            if tair > tmax_ann:
                tmax_ann = tair
            
            if tair < tmin_ann:
                tmin_ann = tair
            
            tavg_ann += tair
        
            self.project_day += 1
        
        # reset date index
        self.project_day = project_day_save
    
        Trange = tmax_ann - tmin_ann
        tavg_ann /= yr_days
        
        # Cool or warm grassland? Definitions are from Botta, Table 1, pg 712.
        # But thresholds are from Foley et al.
        if Trange > 20.0 or tmin_ann < 5.0:
            grass_temp_threshold = 0.0 # cool
        elif Trange <= 20.0 or tmin_ann >= 5.0:
            grass_temp_threshold = 5.0  # warm
       
        # 92% of tmax_ann is the threshold used in grass offset below 
        # Note this has to be done below the range calcs as they use the tmax
        tmax_ann *= 0.92
        
        # Based on White et al. 1997 this is a threshold for grasses so they 
        # have enough accumulated rain. It is essentially a fudge for soil
        # moisture availability and the 15% is somewhat arbitary
        ppt_sum_crit = ppt_sum * 0.15
        
        return (grass_temp_threshold, tmax_ann, ppt_sum_crit)
        
    def calculate_leafon_off(self, daylen, met_data, yr_days):
        self.ini_phen_calcs()
        
        # Krinner et al. 2005, page 26, alternatively Foley et al. 1996 suggests
        # the same value = 100 for both pathways
        if self.control.alloc_model == "GRASSES":
            if self.control.ps_pathway == "C3":
                gdd_thresh = 185
            elif self.control.ps_pathway == "C4":
                gdd_thresh = 400
            (grass_temp_threshold, 
             tmax_ann, 
             ppt_sum_crit) = self.calc_ini_grass_pheno_stuff(met_data, yr_days)
        else:
            gdd_thresh = self.gdd_chill_thresh(self.last_yrs_accumulated_ncd)
        
        if yr_days == 366:
            nov_doy = 306
        else:
            nov_doy = 305
        
        ppt_sum = 0.0
        for d in xrange(yr_days):
            
            Tmean = met_data['tair'][self.project_day]
            Tsoil = met_data['tsoil'][self.project_day]
            ppt_sum += met_data['rain'][self.project_day]
            
            # Calculate ppt total from the next 7 days
            if self.project_day < 358:
                st = self.project_day+1
                en = self.project_day+8
                ppt_sum_next = sum(met_data['rain'][st:en])
            else:
                # i.e. end of year, didn't find this so have no effect
                ppt_sum_next = 0.0
            
            
            # Calculate ppt total from the previous 30 days
            st = self.project_day-30
            en = self.project_day
            ppt_sum_prev = sum(met_data['rain'][st:en])
            
            if d < 362:
                Tsoil_next_3days = ((met_data['tsoil'][self.project_day] +
                                     met_data['tsoil'][self.project_day+1] + 
                                     met_data['tsoil'][self.project_day+2])/3.0)
                                     
                Tair_next_3days = ((met_data['tair'][self.project_day] +
                                    met_data['tair'][self.project_day+1] + 
                                    met_data['tair'][self.project_day+2])/3.0)
            else:
                # i.e. end of year, didn't find this so have no effect
                Tsoil_next_3days = 999.9 
                Tair_next_3days = 999.9

            # Sum the daily mean air temperature above 5degC starting on Jan 1
            self.accum_gdd += self.calc_gdd(Tmean)
            
            
            #
            ## Calculate leaf on
            #
            if self.control.alloc_model == "GRASSES":
                if (self.leaf_on_found == False and 
                    self.accum_gdd >= gdd_thresh and 
                    ppt_sum >= ppt_sum_crit):
                
                    self.leaf_on = d
                    self.leaf_on_found = True
            else:
                if self.leaf_on_found == False and self.accum_gdd >= gdd_thresh:
                    self.leaf_on = d
                    self.leaf_on_found = True
            
            #
            ## Calculate leaf off
            #
            if self.control.alloc_model == "GRASSES":
                if self.leaf_off_found == False:
                    # test for hot and dry conditions 
                    # Based on white et al. 1997
                    if (ppt_sum_prev < 11.4 and
                        ppt_sum_next < 9.7 and 
                        Tmean > tmax_ann): 
                        self.leaf_off_found = True
                        self.leaf_off = d
                        
                    # test for cold offset condition
                    # Leaf drop constraint is based on Foley et al. 1996 as
                    # we dont have access to the Tmin for the constraint from
                    # White et al. This is likely more straightforward anyway
                    elif d > 182 and Tair_next_3days < grass_temp_threshold:
                        self.leaf_off_found = True
                        self.leaf_off = d
            else:
                if (self.leaf_off_found == False and 
                    self.accum_gdd >= gdd_thresh):
                    # I am prescribing that no leaves can fall off before doy=180
                    # Had issue with KSCO simulations where the photoperiod was
                    # less than the threshold very soon after leaf out.
                    if d > 182:
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
        
        if self.leaf_on_found == False or self.leaf_off_found == False:
             raise RuntimeError, "Problem in phenology leaf on/off not found" 
    
        
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
        self.fluxes.crate = (2.0 * self.state.c_to_alloc_croot / 
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
        self.fluxes.cnrate = (2.0 * self.state.n_to_alloc_croot / 
                             self.len_groloss**2)   
        