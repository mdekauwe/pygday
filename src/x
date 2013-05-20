64,65c64
<     def calc_day_growth(self, project_day, fdecay, rdecay, daylen, doy, 
<                         days_in_yr):
---
>     def calc_day_growth(self, project_day, fdecay, rdecay, daylen, doy, days_in_yr):
96,97d94
<         self.update_plant_state(fdecay, rdecay, project_day, doy)
<         self.precision_control()
98a96,97
>         self.update_plant_state(fdecay, rdecay, project_day, doy)
>          
196c195
<              self.state.wtfac_root) = self.sm.calculate_soil_water_fac()
---
>                 self.state.wtfac_root) = self.sm.calculate_soil_water_fac()
253c252
<         #print self.state.alleaf, self.state.alroot, self.state.albranch, self.state.alstem
---
>         
257a257,259
>         # allocation fractions varying per yr
>         
>         """
259c261,270
<         self.state.c_to_alloc_root = self.state.alroot * self.state.cstore
---
>         self.state.c_to_alloc_root = self.state.alroot[0] * self.state.cstore
>         self.state.c_to_alloc_branch = self.state.albranch * self.state.cstore
>         self.state.alstem = 1.0 - self.state.alleaf - self.state.alroot[0]
>         self.state.c_to_alloc_stem = self.state.alstem * self.state.cstore
>         #self.state.c_to_alloc_rootexudate = (self.state.alroot_exudate *
>         #                                        self.state.cstore)
>     
>         # annual available N for allocation to leaf
>         self.state.n_to_alloc_shoot = (self.state.c_to_alloc_shoot *
>                                         self.state.shootnc_yr)
260a272,274
>         """
>         self.state.c_to_alloc_shoot = self.state.alleaf * self.state.cstore
>         self.state.c_to_alloc_root = self.state.alroot * self.state.cstore
265a280
>         
267,268c282,284
<         #self.state.n_to_alloc_branch = self.state.albranch * ntot
<         #self.state.n_to_alloc_stem = self.state.alstem * ntot
---
>         self.state.n_to_alloc_branch = self.state.albranch * self.state.nstore
>         
>         self.state.n_to_alloc_stem = self.state.alstem * self.state.nstore
274c290
<         #                                self.params.ncrfac))
---
>         #                                 self.params.ncrfac))
276,279d291
<         
<         
<         self.state.n_to_alloc_branch = self.state.albranch * self.state.nstore
<         self.state.n_to_alloc_stem = self.state.alstem * self.state.nstore
287a300,313
>         
>         """
>         self.state.c_to_alloc_shoot = self.state.alleaf * self.state.cstore
>         self.state.c_to_alloc_root = self.state.alroot[i] * self.state.cstore
>         self.state.c_to_alloc_branch = self.state.albranch * self.state.cstore
>         self.state.alstem = 1.0 - self.state.alleaf - self.state.alroot[i]
>         sumx = self.state.alleaf + self.state.alstem+ self.state.alroot[i]
>         #print self.state.alleaf, self.state.alstem, self.state.alroot[i]
>         self.state.c_to_alloc_stem = self.state.alstem * self.state.cstore
>         self.state.n_to_alloc_root = (min(self.state.nstore,
>                                           self.state.c_to_alloc_root *
>                                           self.state.rootnc))
>         
>         """
292a319,321
>         # allocate remainder to stem
>         self.state.alstem = (1.0 - self.state.alleaf - self.state.alroot - 
>                              self.state.albranch - self.state.alroot_exudate)
295,296c324,326
<         #self.state.n_to_alloc_branch = self.state.albranch * ntot
<         #self.state.n_to_alloc_stem = self.state.alstem * ntot
---
>         self.state.n_to_alloc_branch = self.state.albranch * self.state.nstore
>         
>         self.state.n_to_alloc_stem = self.state.alstem * self.state.nstore
302c332
<         #                                self.params.ncrfac))
---
>         #                                 self.params.ncrfac))
304,307d333
<         
<         
<         self.state.n_to_alloc_branch = self.state.albranch * self.state.nstore
<         self.state.n_to_alloc_stem = self.state.alstem * self.state.nstore
310a337,344
>        
>         # if we want to put back a floating N:C then we need to have
>         # self.state.c_to_alloc_shoot + self.state.c_to_alloc_stem * some factor
> 
>         # annual available N for allocation to leaf
>         #self.state.n_to_alloc_shoot = (self.state.c_to_alloc_shoot *
>         #                                self.state.shootnc_yr)
> 
337c371
<         self.fluxes.retrans = self.nitrogen_retrans(fdecay, rdecay, doy)
---
>         self.fluxes.retrans = self.nitrogen_retrans(fdecay, rdecay)
394d427
<             
400a434
>             self.fluxes.npbranch = self.fluxes.bnrate * self.state.growing_days[doy]
402d435
<             
404,405d436
<                                   self.state.growing_days[doy])
<             self.fluxes.npbranch = (self.fluxes.bnrate * 
407d437
<             
447c477
<     def nitrogen_retrans(self, fdecay, rdecay, doy):
---
>     def nitrogen_retrans(self, fdecay, rdecay):
465,476c495,508
<             leafretransn = (self.params.fretrans * self.fluxes.lnrate * 
<                             self.state.remaining_days[doy])
<         else:
<             leafretransn = self.params.fretrans * fdecay * self.state.shootn
<         
<         arg1 = (leafretransn +
<                 self.params.rretrans * rdecay * self.state.rootn +
<                 self.params.bretrans * self.params.bdecay *
<                 self.state.branchn)
<         arg2 = (self.params.wretrans * self.params.wdecay *
<                 self.state.stemnmob + self.params.retransmob *
<                 self.state.stemnmob)
---
>             arg1 = (self.fluxes.leafretransn  +
>                     self.params.rretrans * rdecay * self.state.rootn +
>                     self.fluxes.branchretransn)
>             arg2 = (self.params.wretrans * self.params.wdecay *
>                     self.state.stemnmob + self.params.retransmob *
>                     self.state.stemnmob)
>         else:
>             arg1 = (self.params.fretrans * fdecay * self.state.shootn +
>                     self.params.rretrans * rdecay * self.state.rootn +
>                     self.params.bretrans * self.params.bdecay *
>                     self.state.branchn)
>             arg2 = (self.params.wretrans * self.params.wdecay *
>                     self.state.stemnmob + self.params.retransmob *
>                     self.state.stemnmob)
568c600
<             if float_eq(self.state.shoot, 0.0):
---
>             if self.state.shoot == 0.0:
582,613c614,619
<                               (self.state.sla * const.M2_AS_HA / 
<                               (const.KG_AS_TONNES * self.params.cfracts)) -
<                               (self.fluxes.deadleaves + 
<                                self.fluxes.ceaten) *
<                                self.state.lai / self.state.shoot)
<    
<     def precision_control(self, tolerance=1E-08):
<         """ Detect very low values in state variables and force to zero to 
<         avoid rounding and overflow errors """       
<         
<         # C & N state variables 
<         if self.state.shoot < tolerance:
<             self.fluxes.deadleaves += self.state.shoot
<             self.fluxes.deadleafn += self.state.shootn
<             self.fluxes.deadbranch += self.state.branch
<             self.fluxes.deadbranchn += self.state.branchn
<             self.state.shoot = 0.0 
<             self.state.shootn = 0.0 
<             self.state.branch = 0.0
<             self.state.branchn = 0.0
<         
<         if self.state.root < tolerance:
<             self.fluxes.deadrootn += self.state.rootn
<             self.fluxes.deadroot += self.state.root
<             self.state.root = 0.0
<             self.state.rootn = 0.0
<     
<         if self.state.stem < tolerance:     
<             self.fluxes.deadstemn += self.state.stem
<             self.state.stem = 0.0
<             self.state.stemnimm = 0.0
<             self.state.stemnmob = 0.0
---
>                                   (self.state.sla * const.M2_AS_HA / 
>                                   (const.KG_AS_TONNES * self.params.cfracts)) -
>                                   (self.fluxes.deadleaves + 
>                                    self.fluxes.ceaten) *
>                                    self.state.lai / self.state.shoot)
>             
615,616d620
<         # should add check for soil pools - excess goes where?
<         
629c633
<                              self.fluxes.ceaten)
---
>                                 self.fluxes.ceaten)
634c638,644
<         if self.control.deciduous_model:       
---
>         if self.control.deciduous_model:
>             #self.state.shootn += (self.fluxes.npleaf - 
>             #                      self.fluxes.deadleafn - 
>             #                      self.fluxes.neaten )
>             #self.state.branchn += (self.fluxes.npbranch - 
>             #                        self.fluxes.deadbranchn)
>                                     
638c648,659
<                                   self.fluxes.neaten)                        
---
>                                   self.fluxes.neaten)
>             self.state.branchn += (self.fluxes.npbranch - 
>                                    self.fluxes.bnrate * 
>                                    self.state.remaining_days[doy])                         
>             
>             # potential for floating errors to be an issue, so effectively 
>             # zero if pool becomes v.v.small.
>             self.state.shoot = max(0.0, self.state.shoot)                       
>             self.state.shootn = max(0.0, self.state.shootn)
>             self.state.branch = max(0.0, self.state.branch)                       
>             self.state.branchn = max(0.0, self.state.branchn)                       
>                                     
643,645c664,668
<                                 
<         self.state.branchn += (self.fluxes.npbranch - self.params.bdecay *
<                                self.state.branchn)
---
>             self.state.branchn += (self.fluxes.npbranch - self.params.bdecay *
>                                    self.state.branchn)                     
>         
>         
>         
659,663c682,688
<             
<         # This doesn't make sense for the deciduous model because of the ramp
<         # function. The way the deciduous logic works we now before we start
<         # how much N we have to allocate so it is impossible to allocate in 
<         # excess. Therefore this is only relevant for evergreen model.
---
>         
>         
>         # This doesn't make sense for the deciduous model. This is because of 
>         # the ramp function. So there will be a period where we will be above
>         # the ncmaxf, so we will end up just cutting things back. This logic
>         # essentially doesn't fit with the ramp function which allocates 
>         # everything we have to allocate, i.e. I don't think this is necessary
666,669c691,693
<             # if foliage N:C ratio exceeds its max, then nitrogen uptake is cut 
<             # back n.b. new ring n/c max is already set because it is related to 
<             # leaf n:c
<             
---
>             # if foliage N:C ratio exceeds its max, then nitrogen uptake is cut back 
>             # n.b. new ring n/c max is already set because it is related to leaf n:c
>     
722c746,751
<     def calculate_cn_store(self):        
---
>     def calculate_cn_store(self, tolerance=1.0E-05):        
>         cgrowth = (self.fluxes.cpleaf + self.fluxes.cproot + 
>                    self.fluxes.cpbranch + self.fluxes.cpstem)
>         ngrowth = (self.fluxes.npleaf + self.fluxes.nproot + 
>                    self.fluxes.npbranch + self.fluxes.npstemimm + 
>                    self.fluxes.npstemmob)
725,726c754,757
<         self.state.cstore += self.fluxes.npp
<         self.state.nstore += self.fluxes.nuptake + self.fluxes.retrans 
---
>         self.state.cstore += self.fluxes.npp - cgrowth
>         self.state.nstore += self.fluxes.nuptake + self.fluxes.retrans - ngrowth
>         
>         
728c759
<   
---
>         
