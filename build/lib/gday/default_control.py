"""
G'DAY default control flags

Read into the model unless the user changes these at runtime with definitions
in the .INI file

"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (05.09.2011)"
__email__   = "mdekauwe@gmail.com"

startday = 1           # simulation start day
startyear = 1998       # simulation start year
startmonth = 1         # simulation start month
water_model = 1        # Water balance model
assim_model = 7        # 0-4=McM, 5+6=Bewdy, 7=mate
wue_model = 3          #
nuptake_model = 1      # 0=constant uptake, 1=func of N inorgn, 2=depends on rate of soil N availability
trans_model = 1        # 0=trans from WUE, 1=Penman-Monteith, 2=Priestley-Taylor
fixroot = 0            #
fixleafnc = 0          # fixed leaf N C ?
leuning_func = 0       #
sel_noprod1 = 0        # switch off 'Resid' flow
sel_noprod2 = 0        # switch off 'Resid' flow
sel_noprod3 = 0        # switch off 'Resid' flow
sel_noprod4 = 0        # switch off 'Resid' flow
passiveconst = 0       # hold passive pool at passivesoil
print_options = 1      # 0=every timestep, 1=end of run
grazing = 0            # Is foliage grazed 0=No, 1=Yes
use_eff_nc = 0         # use constant leaf n:c for  metfrac s
strfloat = 0           # Structural pool input N:C floats
use_leuning = 0
co2_conc = 0           # 0=ambient, 1=elevated
fixed_stem_nc = 1      # 1=fixed stem N:C, 0=vary stem N:C with foliage
