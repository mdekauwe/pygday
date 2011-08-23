"""
G'DAY default control flags

Read into the model unless the user changes these at runtime with definitions 
in the .INI file

"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (14.02.2011)"
__email__   = "mdekauwe@gmail.com"


#set control parameters
startday = 1           # simulation start day
startyear = 1998       # simulation start year
startmonth = 1         # simulation start month
water_model = 1        # Water balance model
model_number = 7       # assimilation model
wue_model = 3          #  
constuptake = 0        # Hold N uptake constant at Nuptakez ? BM
fixroot = 0            #
fixleafnc = 0          # fixed leaf N C ?
leuning_func = 0       #
sel_noprod1 = 0        # switch off 'Resid' flow
sel_noprod2 = 0        # switch off 'Resid' flow
sel_noprod3 = 0        # switch off 'Resid' flow
sel_noprod4 = 0        # switch off 'Resid' flow
passiveconst = 0       # hold passive pool at passivesoil
print_options = 1      # 0=every timestep, 1=end of run
grazing = 0            # Is foliage grazed (Y or N)
use_eff_nc = 0          # use constant leaf n c for  metfrac s
strfloat = 0           # Structural pool input N:C floats 
consnuptake = 0        # const N uptake (t/ha/yr) if Constuptake
use_leuning = 0
scale_nup_with_availb = 0
trans_model = 1        # 0=trans from WUE, 1=Penman-Monteith, 2=Priestley-Taylor
co2_conc = 0