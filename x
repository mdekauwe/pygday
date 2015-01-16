diff --cc src/soil_cn_model.py
index 72155db,da6eebc..0000000
--- a/src/soil_cn_model.py
+++ b/src/soil_cn_model.py
@@@ -1041,19 -1043,29 +1041,17 @@@ class NitrogenSoilFlows(object)
                                 self.fluxes.nlittrelease)
          
      def calc_root_exudation_uptake_of_N(self):
 -        """ When N mineralisation is large enough to allow a small amount of N
 -        immobilisation, the amount of N which enters the active pool is 
 -        calculated according to REXC divided by the CN of the active pool. When
 -        exudation enters the active pool, the CN ratio of the exudates drops 
 -        from REXC/REXN to the CN of the active pool. Which is consistent with 
 -        the CENTURY framework, where C flows between pools lead to either
 -        mineralisation (N gain) or immobilisation (N loss) due to differences
 -        in the CN ratio of the outgoing and incoming pools.
 -        
 -        The amount of N added to the active pool is independent of the CUE of 
 -        the microbial pool in response to root exudation (REXCUE).
 -        
 -        """
          active_CN_ratio = self.state.activesoil / self.state.activesoiln
--        
-         C_to_active_pool = self.fluxes.root_exc * (1.0 - self.fluxes.rexc_cue)
-         N_to_active_pool = C_to_active_pool / active_CN_ratio
 -        # Need to account for the increase in available N
 -        N_miss = (max(0.0, self.fluxes.root_exc / active_CN_ratio - 
 -                           self.fluxes.root_exn))
 -        
++        N_to_active_pool = self.fluxes.root_exc / active_CN_ratio
 +    
 +        # N immobilisation (loss) due to REXN sequestration in the active pool
-         N_miss = (max(0.0, C_to_active_pool / active_CN_ratio) - 
++        N_miss = (max(0.0, self.fluxes.root_exc / active_CN_ratio) - 
 +                      self.fluxes.root_exn)
 +    
 +        # N added to the active pool is independent of the CUE of the microbial
 +        # pool in response to root exudation
          if N_miss < self.fluxes.nmineralisation:
 -            
 -            N_to_active_pool = self.fluxes.root_exc / active_CN_ratio
 -            
 +        
              # update active pool
              self.state.activesoiln += N_to_active_pool
          
