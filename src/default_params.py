"""
G'DAY default model parameters

Read into the model unless the user changes these at runtime with definitions
in the .INI file

"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (14.02.2011)"
__email__   = "mdekauwe@gmail.com"

# N-Cycle off
prescribed_leaf_NC = 0.03

# Age parameters - if equal to each other this is turned off
ageold           = 1000.0         #Plant age when max leaf N C ratio is lowest
ageyoung         = 0.0            #Plant age when max leaf N C ratio is highest

# miscellaneous 
latitude         = 39.11999 #latitude (degrees, negative for south)
albedo           = 0.18     #albedo
liteffnc         = 0.0
passivesoilz = 1.0          # constant vals
passivesoilnz = 1.0         # constant vals
d0 = 0.0
d1 = 0.0
a1 = 0.0

# phenology
previous_ncd = 17 # In the first year we don't have last years data, so I have precalculated the average of all the november-jan chilling values
store_transfer_len = None

# set environment parameters
finesoil          = 0.5            # clay+silt fraction

# light interception parameters
direct_frac       = 0.5            # direct beam fraction of incident radiation
kext              = 0.5            # extinction coefficient
lai_cover         = 0.5            # LAI when max cover fraction is reached (m2 (leaf) m-2 (ground) ~ 2.5
sla               = 3.9            # specific leaf area (m2 one-sided/kg DW)
slazero           = 3.9            # (if equal slamax=no effect) specific leaf area new fol at zero leaf N/C (m2 one-sided/kg DW)
slamax            = 3.9            # (if equal slazero=no effect) specific leaf area new fol at max leaf N/C (m2 one-sided/kg DW)

# set photosynthetic parameters
alpha_j           = 0.26           # initial slope of rate of electron transport, used in calculation of quantum yield. Value calculated by Belinda
alpha_c4          = 0.06           # quantium efficiency for C4 plants has no Ci and temp dependancy, so if a fixed constant.
cfracts           = 0.5            # carbon fraction of dry biomass
cue               = 0.5            # carbon use efficiency, or the ratio of NPP to GPP
delsj             = 644.4338       # Deactivation energy for electron transport (J mol-1 k-1)
eac               = 79430.0        # Activation energy for carboxylation [J mol-1]
eao               = 36380.0        # Activation energy for oxygenation [J mol-1]
eag               = 37830.0        # Activation energy at CO2 compensation point [J mol-1]
eaj               = 43790.0        # Activation energy for electron transport (J mol-1)
eav               = 51560.0        # Activation energy for Rubisco (J mol-1)
edj               = 2e+05          # Deactivation energy for electron transport (J mol-1)
gamstar25         = 42.75          # Base rate of CO2 compensation point at 25 deg C [umol mol-1]
growth_efficiency = 0.7            # growth efficiency (yg) - used only in Bewdy
jmaxna            = 40.462         # slope of the reln btween jmax and leaf N content (g N m-2) - (umol/g n/s) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010.
jmaxnb            = 13.691         # intercept of jmax vs n (umol/g n/s) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010.
jmax              = -999.9         # maximum rate of electron transport (umol m-2 s-1)
kc25              = 404.9          # Base rate for carboxylation by Rubisco at 25degC [mmol mol-1]
ko25              = 278400.0       # Base rate for oxygenation by Rubisco at 25degC [umol mol-1]. Note value in Bernacchie 2001 is in mmol!!
kq10              = 0.08           # exponential coefficient for Rm vs T
measurement_temp  = 25.0           # temperature Vcmax/Jmax are measured at, typical 25.0 (celsius) 
oi                = 205000.0       # intercellular concentration of O2 [umol mol-1]
theta             = 0.7            # curvature of photosynthetic light response curve
vcmaxna           = 20.497         # slope of the reln btween vcmax and leaf N content (g N m-2) - (umol/g n/s) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010.
vcmaxnb           = 8.403          # intercept of vcmax vs n (umol/g n/s) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010.
vcmax             = -999.9         # maximum rate of carboxylation (umol m-2 s-1) 

# water model parameters
dz0v_dh           = 0.123    # Rate of change of vegetation roughness length for momentum with height. Value from Jarvis? for conifer 0.075
displace_ratio    = 0.67     # Value for coniferous forest (0.78) from Jarvis et al 1976, taken from Jones 1992 pg 67. More standard assumption is 2/3
z0h_z0m           = 0.1      # Assume z0m = z0h, probably a big assumption [as z0h often < z0m.], see comment in code!! But 0.1 might be a better assumption
wcapac_root       = 240.0    # Max plant avail soil water -root zone, i.e. total (mm) (smc_sat-smc_wilt) * root_depth (750mm) = [mm (water) / m (soil depth)]
wcapac_topsoil    = 100.0    # Max plant avail soil water -top soil (mm)
rooting_depth     = 20000.0  # Rooting depth (mm)
topsoil_depth     = 450.0    # Topsoil depth (mm)
topsoil_type      = None
rootsoil_type     = None
ctheta_topsoil    = None     # Fitted parameter based on Landsberg and Waring
ntheta_topsoil    = None     # Fitted parameter based on Landsberg and Waring
ctheta_root       = None     # Fitted parameter based on Landsberg and Waring
ntheta_root       = None     # Fitted parameter based on Landsberg and Waring
fractup_soil      = 0.5      #fraction of uptake from top soil layer
wetloss           = 0.5      #daily rainfall lost per lai (mm/day)
rfmult            = 1.0
intercep_frac     = 0.15     # Maximum intercepted fraction, values in Oishi et al 2008, AFM, 148, 1719-1732 ~13.9% +/- 4.1, so going to assume 15% following Landsberg and Sands 2011, pg. 193.
max_intercep_lai  = 3.0      # canopy LAI at which interception is maximised.
qs                = 1.0      # exponent in water stress modifier, =1.0 JULES type representation, the smaller the values the more curved the depletion. 
g1                = 4.8      # stomatal conductance parameter: Slope of reln btw gs and assimilation (fitted by species/pft).
b_topsoil         = None
b_root            = None
psi_sat_root      = None     # MPa
psi_sat_topsoil   = None     # MPa
theta_sat_root     = None
theta_sat_topsoil  = None

# set carbon allocation parameters & allometric parameters
c_alloc_fmax = 0.25    # allocation to leaves at leaf n_crit. If using allometric model this is the max alloc to leaves
c_alloc_fmin = 0.25    # allocation to leaves at zero leaf n/c. If using allometric model this is the min alloc to leaves
c_alloc_rmax = 0.05    # allocation to roots at root n_crit. If using allometric model this is the max alloc to fine roots
c_alloc_rmin = 0.05    # allocation to roots at zero root n/c. If using allometric model this is the min alloc to fine roots
c_alloc_bmax = 0.2     # allocation to branches at branch n_crit. If using allometric model this is the max alloc to branches
c_alloc_bmin = 0.2     # allocation to branches at zero branch n/c. If using allometric model this is the min alloc to branches
c_alloc_cmax = 0.2     # allocation to coarse roots at n_crit. If using allometric model this is the max alloc to coarse roots
heighto      = 4.826   # constant in avg tree height (m) - stem (t C/ha) reln
htpower      = 0.35    # Exponent in avg tree height (m) - stem (t C/ha) reln
height0      = 5.0     # Height when leaf:sap area ratio = leafsap0 (trees)
height1      = 30.0    # Height when leaf:sap area ratio = leafsap1 (trees)
leafsap0     = 7500.0  # leaf area  to sapwood cross sectional area ratio when Height = Height0 (mm^2/mm^2)
leafsap1     = 2700.0  # leaf to sap area ratio when Height = Height1 (mm^2/mm^2)
branch0      = 5.61    # constant in branch-stem allometry (trees)
branch1      = 0.346   # exponent in branch-stem allometry
croot0       = 0.34     # constant in coarse_root-stem allometry (trees)
croot1       = 0.84     # exponent in coarse_root-stem allometry
targ_sens    = 0.5     # sensitivity of allocation (leaf/branch) to track the target, higher values = less responsive.
density      = 420.0   # sapwood density kg DM m-3 (trees)
nf_min       = 0.005   # leaf N:C minimum N concentration which allows productivity
nf_crit      = 0.015   # leaf N:C below which N availability limits productivity 

#set nitrogen allocation parameters
ncwnewz          = 0.003          #New stem ring N:C at zero leaf N:C (mobile)
ncwnew           = 0.003          #New stem ring N:C at critical leaf N:C (mob)
ncwimmz          = 0.003          #Immobile stem N C at zero leaf N C
ncwimm           = 0.003          #Immobile stem N C at critical leaf N C
ncbnewz          = 0.003          #new branch N C at zero leaf N C
ncbnew           = 0.003          #new branch N C at critical leaf N C
nccnewz          = 0.003          #new coarse root N C at zero leaf N C
nccnew           = 0.003          #new coarse root N C at critical leaf N C

ncrfac           = 0.8            #N:C of fine root prodn / N:C c of leaf prodn
ncmaxfyoung      = 0.04           #max N:C ratio of foliage in young stand, if the same as old=no effect
ncmaxfold        = 0.04           #max N:C ratio of foliage in old stand, if the same as young=no effect
ncmaxr           = 0.03           #max N:C ratio of roots
retransmob       = 0.0            #Fraction stem mobile N retranscd (/yr)
fhw              = 0.8            # n:c ratio of stemwood - immobile pool and new ring
nmin             = 0.95           # (bewdy) minimum leaf n for +ve p/s (g/m2)

# grazing parameters
fracteaten       = 0.5   #Fractn of leaf prodn eaten by grazers
fracfaeces       = 0.3   #Fractn of grazd C that ends up in faeces (0..1)
ligfaeces        = 0.25  #Faeces lignin as fractn of biomass
faecescn         = 25.0  #Faeces C:N ratio
fractosoil       = 0.85  #Fractn of grazed N recycled to soil:faeces+urine

#nitrogen cycling parameters
nuptakez         = 0.0    # constant N uptake per year (1/yr)
rateuptake       = 5.7    # rate of N uptake from mineral N pool (/yr) from here? http://face.ornl.gov/Finzi-PNAS.pdf
rateloss         = 0.5    # Rate of N loss from mineral N pool (/yr)
fretrans         = 0.5    # foliage n retranslocation fraction - 46-57% in young E. globulus trees - see Corbeels et al 2005 ecological modelling 187, pg 463. Roughly 50% from review Aerts '96
rretrans         = 0.0    # root n retranslocation fraction
cretrans         = 0.0    # coarse root n retranslocation fraction
bretrans         = 0.0    # branch n retranslocation fraction
wretrans         = 0.0    # mobile wood N retranslocation fraction
kr               = 0.5    # N uptake coefficent (0.05 kg C m-2 to 0.5 tonnes/ha) see Silvia's PhD, Dewar and McM, 96.
nmax             = 0.24   # Maximum rate of N uptake by the vegetation (tonnes/hectare/year)
knl              = 0.01   # Concentration of N available at which N uptake proceeds at one half its maximum rate (tonnes/hectare).
ac               = 0.5    # relative amount of effort allocation to C vs. N uptake[0,1].
adapt            = 0.012  # rate of adaption of the vegetation to changing nutrient abundance (yr -1)

# set litter parameters
fdecay           = 0.5            #foliage turnover rate (1/yr)
fdecaydry        = 0.5            #Foliage turnover rate - dry soil (1/yr)
rdecay           = 0.5            #root turnover rate (1/yr)
rdecaydry        = 0.5            #root turnover rate - dry soil (1/yr)
bdecay           = 0.03           #branch and large root turnover rate (1/yr)
wdecay           = 0.02           #wood turnover rate (1/yr)
crdecay          = 0.02           #coarse roots turnover rate (1/yr)
watdecaydry      = 0.0            #water fractn for dry litterfall rates
watdecaywet      = 0.1            #water fractn for wet litterfall rates
sapturnover      = 0.1            #Sapwood turnover rate: conversion of sapwood to heartwood (1/yr)
ligshoot         = 0.25           #lignin-to-biomass ratio in leaf litter - Value in smith et al. 2013 = 0.2, note subtly difference in eqn C9.
ligroot          = 0.25           #lignin-to-biomass ratio in root litter - Value in Smith et al. 2013 = 0.16, note subtly difference in eqn C9.
structcn         = 150.0          #C:N ratio of structural bit of litter input
structrat        = 0.0            #structural input n:c as fraction of metab

#set decomposition parameters - converted from yr to day in model!
kdec1            = 3.965571       #surface structural decay rate (1/yr)
kdec2            = 14.61          #surface metabolic decay rate (1/yr)
kdec3            = 4.904786       #soil structural decay rate (1/yr)
kdec4            = 18.262499      #soil metabolic decay rate(1/yr)
kdec5            = 7.305          #active pool decay rate (1/yr)
kdec6            = 0.198279       #slow pool decay rate (1/yr)
kdec7            = 0.006783       #passive pool decay rate (1/yr)

# Set N:C ratios of soil pools [units: g/m2]
actncmax         = 0.333333  # Active pool (=1/3) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
actncmin         = 0.066667  # Active pool (=1/15) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
slowncmax        = 0.066667  # Slow pool (=1/15) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
slowncmin        = 0.025     # Slow pool (=1/40) N:C of new SOM - when Nmin=Nmin0" [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
passncmax        = 0.142857  # Passive pool (=1/7) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
passncmin        = 0.1       # Passive pool (=1/10) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
nmincrit         = 2.0       # Critical mineral N pool at max soil N:C (g/m2) (Parton et al 1993, McMurtrie et al 2001).
nmin0            = 0.0       # Mineral N pool corresponding to Actnc0,etc (g/m2)

# root model stuff
d0x    = 0.35   # Length scale for exponential decline of Umax(z)
r0     = 0.1325 # root C at half-maximum N uptake (kg C/m3)

# Disturbance
return_interval = 10 # yrs
disturbance_doy = 1
burn_specific_yr = None
hurricane_doy = None
hurricane_yr = None

# priming stuff
root_exu_CUE = None
prime_y = 0.0
prime_z = 0.0


#============== Not publicly accessible to the user ==========================#
# decay rates
decayrate = [None] * 7

# metabolic pool C fractions
fmfaeces = 0.0
fmleaf = 0.0
fmroot = 0.0
faecesn = 0.0
#==============================================================================#
