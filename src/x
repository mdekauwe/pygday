9c9
< from math import fabs
---
> import math
16c16
< from update_pools import CarbonSoilPools, NitrogenSoilPools
---
> from update_pools import CarbonPools, NitrogenPools
113c113
<         self.cpl = CarbonSoilPools(self.control, self.params, self.state, 
---
>         self.cpl = CarbonPools(self.control, self.params, self.state, 
115c115
<         self.npl = NitrogenSoilPools(self.control, self.params, self.state, 
---
>         self.npl = NitrogenPools(self.control, self.params, self.state, 
144,146c144,146
<         while (fabs(prev_plantc - self.state.plantc) > tolerance and 
<                fabs(prev_soilc - self.state.soilc) > tolerance and 
<                fabs(prev_litterc - self.state.litterc) > tolerance): 
---
>         while (math.fabs(prev_plantc - self.state.plantc) > tolerance and 
>                math.fabs(prev_soilc - self.state.soilc) > tolerance and 
>                math.fabs(prev_litterc - self.state.litterc) > tolerance): 
177,178c177,178
<                 self.pg.calc_day_growth(project_day, fdecay, rdecay, 
<                                         daylen[doy], doy, float(days_in_year))
---
>                 self.pg.grow(project_day, fdecay, rdecay, daylen[doy], doy, 
>                         float(days_in_year))
184c184
<                 # Update model soil C&N pools
---
>                 # soil model - update pools
186c186
<                 self.npl.calculate_npools(cact, cslo, cpas)
---
>                 self.npl.calculate_npools(cact, cslo, cpas, project_day)
191,192c191,193
<                 #if self.spin_up == False:
<                 #    print self.fluxes.gpp * 100, self.state.lai
---
>                 if self.spin_up == False:
>                     print self.fluxes.gpp * 100, self.state.lai, self.fluxes.nrootexudate
>                 
196,197c197,198
<                 #if self.control.print_options == 0:
<                 #    self.save_daily_outputs(yr, doy+1)
---
>                 if self.control.print_options == 0:
>                     self.save_daily_outputs(yr, doy+1)
411,412c412,413
<     fname = "/Users/mdekauwe/research/NCEAS_face/GDAY_duke_simulation/params/NCEAS_dk_youngforest.cfg"
<     #fname = "test.cfg"
---
>     #fname = "/Users/mdekauwe/research/NCEAS_face/GDAY_duke_simulation/params/NCEAS_dk_youngforest.cfg"
>     fname = "test.cfg"
438,439c439,440
<     #main()
<     profile_main()
---
>     main()
>     #profile_main()
