14c14
< from litter_production import Litter
---
> from litter_production import LitterProduction
19a20
> 
81,84c82
<         
<         # build list of variables to print
<         (self.print_state, self.print_fluxes) = self.pr.get_vars_to_print()
<         
---
> 
93a92,103
>         # set initial lai -> m2/m2
>         self.state.lai = (self.params.slainit * const.M2_AS_HA /
>                           const.KG_AS_TONNES / self.params.cfracts *
>                           self.state.shoot)
> 
>         # Specific leaf area (m2 onesided/kg DW)
>         self.state.sla = self.params.slainit
> 
>         self.time_constants = ['rateuptake', 'rateloss', 'retransmob',
>                                 'fdecay', 'fdecaydry', 'rdecay', 'rdecaydry',
>                                 'bdecay', 'wdecay', 'kdec1', 'kdec2', 'kdec3',
>                                 'kdec4', 'kdec5', 'kdec6', 'kdec7', 'nuptakez']
101c111,112
<         self.lf = Litter(self.control, self.params, self.state, self.fluxes)
---
>         self.lf = LitterProduction(self.control, self.params, self.state,
>                                    self.fluxes)
109c120
<         
---
> 
111c122
<             self.pg.initialise_deciduous_model()
---
>             self.initialise_deciduous_model()
121,124d131
<         self.state.sla = self.params.slainit # Specific leaf area (m2/kg DW)
<         self.state.lai = (self.params.slainit * const.M2_AS_HA /
<                           const.KG_AS_TONNES / self.params.cfracts *
<                           self.state.shoot)
126,127d132
<         
<         
156,157c161
<     
<     #@profile
---
> 
161,168c165,168
<         
<         # figure out the number of years for simulation and the number of
<         # days in each year
<         years = uniq(self.met_data["year"])
<         days_in_year = [self.met_data["year"].count(yr) for yr in years]
<         
<         for i, yr in enumerate(years):
<             daylen = calculate_daylength(days_in_year[i], self.params.latitude)
---
>         project_day2 = 0
>         for yr in uniq(self.met_data["year"]):
>             days_in_year = len([x for x in self.met_data["year"] if x == yr])
>             daylen = calculate_daylength(days_in_year, self.params.latitude)
172c172
<                                             days_in_year[i], project_day)
---
>                                             days_in_year, project_day)
175c175
<             for doy in xrange(days_in_year[i]):
---
>             for doy in xrange(days_in_year):
178c178
<                 (fdecay, rdecay) = self.lf.calculate_litter(doy)
---
>                 (fdecay, rdecay) = self.lf.calculate_litter_flows(doy)
182,183c182
<                                         daylen[doy], doy, 
<                                         float(days_in_year[i]))
---
>                                         daylen[doy], doy, float(days_in_year))
194c193,194
<                 self.day_end_calculations(project_day, days_in_year[i])
---
>                 self.day_end_calculations(project_day, days_in_year)
> 
196d195
<                 
201c200
<                 print self.fluxes.gpp * 100
---
>                 
215c214
<                 self.pg.allocate_stored_c_and_n()
---
>                 self.allocate_stored_c_and_n()
241a241,242
> 
> 
250,255d250
<         time_constants = ['rateuptake', 'rateloss', 'retransmob',
<                           'fdecay', 'fdecaydry', 'rdecay', 'rdecaydry',
<                           'bdecay', 'wdecay', 'kdec1', 'kdec2', 'kdec3',
<                           'kdec4', 'kdec5', 'kdec6', 'kdec7', 'nuptakez']
<         conv = const.NDAYS_IN_YR
<         
257,258c252,254
<             for i in time_constants:
<                 setattr(self.params, i, getattr(self.params, i) / conv)
---
>             for i in self.time_constants:
>                 setattr(self.params, i, getattr(self.params, i) /
>                         const.NDAYS_IN_YR)
260,261c256,258
<             for i in time_constants:
<                 setattr(self.params, i, getattr(self.params, i) * conv)
---
>             for i in self.time_constants:
>                 setattr(self.params, i, getattr(self.params, i) *
>                         const.NDAYS_IN_YR)
332a330,385
>     def initialise_deciduous_model(self):
>         # Divide up NPP based on annual allocation fractions
> 
>         self.state.c_to_alloc_shoot = (self.state.alleaf *
>                                         self.state.cstore)
> 
> 
>         self.state.c_to_alloc_root = (self.state.alroot *
>                                         self.state.cstore)
> 
>         self.state.c_to_alloc_branch = (self.state.albranch *
>                                         self.state.cstore)
>         self.state.c_to_alloc_stem = (self.state.alstem *
>                                         self.state.cstore)
>         #self.state.c_to_alloc_rootexudate = (self.state.alroot_exudate *
>         #                                        self.state.cstore)
> 
>         # annual available N for allocation to leaf
>         self.state.n_to_alloc_shoot = (self.state.c_to_alloc_shoot *
>                                         self.state.shootnc_yr)
> 
> 
> 
>     def allocate_stored_c_and_n(self):
>         """
>         At the end of the year allocate everything for the coming year
>         based on stores from the previous year avaliable N for allocation
>         """
>         self.state.c_to_alloc_shoot = (self.state.alleaf *
>                                         self.state.cstore)
>         self.state.c_to_alloc_root = (self.state.alroot *
>                                         self.state.cstore)
>         self.state.c_to_alloc_branch = (self.state.albranch *
>                                         self.state.cstore)
>         self.state.c_to_alloc_stem = (self.state.alstem *
>                                         self.state.cstore)
> 
>         self.state.n_to_alloc_root = (min(self.state.nstore,
>                                           self.state.c_to_alloc_root *
>                                           self.state.rootnc))
> 
>         # constant N:C of foliage during the growing season(kgN kg-1C)
>         self.state.shootnc_yr = ((self.state.nstore -
>                                   self.state.n_to_alloc_root) /
>                                   (self.state.c_to_alloc_shoot))
>         # if we want to put back a floating N:C then we need to have
>         # self.state.c_to_alloc_shoot + self.state.c_to_alloc_stem * some factor
> 
>         # annual available N for allocation to leaf
>         self.state.n_to_alloc_shoot = (self.state.c_to_alloc_shoot *
>                                         self.state.shootnc_yr)
> 
>         total = (self.state.c_to_alloc_shoot + self.state.c_to_alloc_root + 
>                  self.state.c_to_alloc_stem)
>         
> 
345,348c398,408
<         for var in self.print_state:
<             output.append(getattr(self.state, var))
<         for var in self.print_fluxes:
<             output.append(getattr(self.fluxes, var))
---
>         for var in self.print_opts:
>             try:
>                 if hasattr(self.state, var):
>                     value = getattr(self.state, var)
>                     output.append(value)
>                 else:
>                     value = getattr(self.fluxes, var)
>                     output.append(value)
>             except AttributeError:
>                 err_msg = "Error accessing var to print: %s" % var
>                 raise AttributeError, err_msg
