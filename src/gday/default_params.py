"""
G'DAY default model parameters

Read into the model unless the user changes these at runtime with definitions 
in the .INI file

"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (14.02.2011)"
__email__   = "mdekauwe@gmail.com"


#set environment parameters
finesoil         = 0.5            #clay+silt fraction 
ambient_co2      = 350.0          #ambient co2 concentration (ppm)
co2_unlim_c_prod = 1.0            #unlimited C prodn (t/ha/yr)
co2_effect_on_wue = 0.705          #CO2 effect on WUE (btwn 0 and 1# try )

#set photosynthetic parameters
n_crit            = 0.04           #critical leaf n/c ratio for no prod n loss, Kirschbaum etal 1994, pg 1086, =0.016/0.45
ncpower           = 1.0            #n/c power for n-limitation (try 1)
ci_ca_ratio       = 0.666667       #ci/ca ratio
kext              = 0.5            #extinction coefficient
slainit           = 3.9            #specific leaf area (m2/kg_biomass)
slazero           = 3.9            #specific leaf area new fol at zero leaf N/C (m2/kg)
slamax            = 3.9            #specofic leaf area new fol at max leaf N/C (m2/kg)
lai_cover         = 0.5            #Onesided LAI correspdg to complete ground cover
cfracts           = 0.5            #carbon fraction of dry biomass 
nmin              = 0.95           #minimum leaf n for +ve p/s (g/m2)
jmaxn             = 56.0           #slope of jmax vs n (umol/g n/s)
vcmaxn            = 28.0           #slope of vcmax vs n (umol/g n/s)
growth_efficiency = 0.7            #growth efficiency (yg)
alpha_j           = 0.3            # initial slope of rate of electron transport
alpha             = 0.05           # quantum yield (mol mol-1) used in mate
direct_frac       = 0.5            #direct beam fraction of incident radiation
kq10              = 0.08           #exponential coefficient for Rm vs T
epsilon           = 1.0
eav               = 51560.0        # Activation energy for Rubisco (J mol-1)
eaj               = 43790.0        # Activation energy for electron transport (J mol-1)
edj               = 2e+05          # Deactivation energy fro electron transport (J mol-1)
delsj             = 644.4338       # J mol-1 k-1
theta             = 0.7            # curvature of photosynthetic light response curve
cue               = 0.5            # carbon use efficiency, or the ratio of NPP to GPP
g1                = 4.8            # fitted param, proportional to sqrt(gamma_star * marginal water cost of carbon)




#set carbon allocation & grazing parameters
callocf_crit     = 0.25  #allocation to leaves at leaf  n_crit 
callocf          = 0.25  #allocation to leaves at zero leaf n/c
callocr_crit     = 0.05 #allocation to roots at leaf  n_crit 
callocr          = 0.05 #allocation to roots at zero leaf n/c
callocb_crit     = 0.2  #allocation to branches at leaf  n_crit 
callocb          = 0.2  #allocation to branches at zero leaf n/c
fracteaten       = 0.5  #Fractn of leaf prodn eaten by grazers
fracfaeces       = 0.3  #Fractn of grazd C that ends up in faeces (0..1)
ligfaeces        = 0.25 #Faeces lignin as fractn of biomass
faecescn         = 25.0 #Faeces C:N ratio
fractosoil       = 0.85 #Fractn of grazed N recycled to soil:faeces+urine

#nitrogen cycling params
#ninflow          = 0.016         #Atmospheric N input (t/ha/yr) (/yr)
rateuptake       = 5.7            #rate of N uptake from mineral N pool (/yr) from here? http://face.ornl.gov/Finzi-PNAS.pdf
rateloss         = 0.5            #Rate of N loss from mineral N pool (/yr)
fretrans         = 0.5            #foliage n retranslocation fraction
rretrans         = 0.0            #root n retranslocation fraction
bretrans         = 0.0            #branch n retranslocation fraction
wretrans         = 0.0            #mobile wood N retranslocation fraction
uo = 2.737850787E-4 # Supply rate of available N (0.01 kg N m-2 yr-1 to t/ha/day)
kr = 0.5 # N uptake coefficent (0.05 kg C m-2 to 0.5 tonnes/ha)

#set nitrogen allocation parameters
ncwnew           = 0.003          #New stem ring N:C at zero leaf N:C (mobile)
ncwnew_crit      = 0.003          #New stem ring N:C at critical leaf N:C (mob)
ncwimm           = 0.003          #Immobile stem N C at zero leaf N C
ncwimm_crit      = 0.003          #Immobile stem N C at critical leaf N C
ncbnew           = 0.003          #new branch N C at zero leaf N C
ncbnew_crit      = 0.003          #new branch N C at critical leaf N C
ncrfac           = 0.8            #N:C of root prodn / N:C c of leaf prodn
ageold           = 1000.0         #Plant age when max leaf N C ratio is lowest 
ageyoung         = 0.0            #Plant age when max leaf N C ratio is highest
ncmaxfyoung      = 0.04           #max N:C ratio of foliage in young stand
ncmaxfold        = 0.04           #max N:C ratio of foliage in old stand
ncmaxr           = 0.03           #max N:C ratio of roots
retransmob       = 0.0            #Fraction stem mobile N retranscd (/yr)

#set litter parameters
fdecay           = 0.5     #foliage decay rate (1/yr)
fdecaydry        = 0.5      #Foliage decay rate - dry soil (1/yr)
rdecay           = 0.5      #root decay rate (1/yr)
rdecaydry        = 0.5    #root decay rate - dry soil (1/yr)
bdecay           = 0.03       #branch and large root decay rate (1/yr)
wdecay           = 0.02     #wood decay rate (1/yr)
watdecaydry      = 0.0            #water fractn for dry litterfall rates
watdecaywet      = 0.1            #water fractn for wet litterfall rates
ligshoot         = 0.25           #shoot litter lignin as fraction of c
ligroot          = 0.25           #root litter lignin as fraction of c
brabove          = 0.5           #above-ground fraction of branch pool litter
structcn         = 150.0          #C:N ratio of structural bit of litter input
metfrac0         = 0.85           #litter metabolic fraction
metfrac1         = -0.018         #litter metabolic fraction
structrat        = 0.0            #structural input n:c as fraction of metab

#set decomposition parameters
kdec1            = 3.965571    #surface structural decay rate (1/yr)
kdec2            = 14.61       #surface metabolic decay rate (1/yr)
kdec3            = 4.904786   #soil structural decay rate (1/yr)
kdec4            = 18.262499      #soil metabolic decay rate(1/yr)
kdec5            = 7.305        #active pool decay rate (1/yr)
kdec6            = 0.198279    #slow pool decay rate (1/yr)
kdec7            = 0.006783    #passive pool decay rate (1/yr)

# Set N:C ratios of soil pools [units: g/m2]
actncmax         = 0.333333  #Active pool N:C ratio of new SOM - maximum [units: g/m2]
actnc0           = 0.066667  #Active pool N:C of new SOM - when Nmin=Nmin0 [units: g/m2]
slownc0          = 0.025  #Slow pool N:C of new SOM - when Nmin=Nmin0" [units: g/m2]
slowncmax        = 0.066666  #Slow pool N:C ratio of new SOM - maximum [units: g/m2]
passncmax        = 0.142857  #Passive pool N:C ratio of new SOM - maximum [units: g/m2]
passnc0          = 0.1  #Passive pool N:C of new SOM - when Nmin=Nmin0 [units: g/m2]
nmincrit         = 2.0  # Critical mineral N pool at max soil N:C (g/m2)
nmin0            = 0.0  # Mineral N pool corresponding to Actnc0,etc (g/m2)

#set water model parameters
wcapac_root       = 240.0 #Max plant avail soil water -root zone, i.e. total (mm) (smc_sat-smc_wilt) * root_depth (750mm) = [mm (water) / m (soil depth)]
wcapac_topsoil    = 100.0 #Max plant avail soil water -top soil (mm)
fwpmax            = 0.52 #Fractional water content at saturation point (max production)
fwpmin            = 0.2  #Fractional water content at wilting point (no production)
fractup_soil      = 0.5            #fraction of uptake from top soil layer
extraction        = 0.007          #water extractn by unit root mass(ha/tC/d)
wue0              = 3.0            #WUE if VPD=1kPa, CO2=350ppm (gC*kPa/kgH2O)
wetloss           = 0.5  #daily rainfall lost per lai (mm/day)
rfmult            = 1.0

# misc
latitude         = 39.11999 #latitude (degrees, negative for south)
albedo           = 0.18            #albedo
liteffnc         = 0.0
canht            = 17.0   # Canopy height increased from 16m in 2001 to 18m in 2004 at Duke
dz0v_dh          = 0.075  # Rate of change of vegetation roughness length for momentum with height. 
displace_ratio   = 0.78   # Value for coniferous forest from Jarvis et al 1976, taken from Jones 1992 pg 67.  



vpd = 0.5
age = 0.0       #Initial stand age (years)
nuptakez = 0.0  # (1/yr)
passivesoilz = 1.0  # constant vals
passivesoilnz = 1.0 # constant vals
d0 = 0.0
d1 = 0.0
a1 = 0.0




#============== Not publicly accessible to the user ==========================#

# decay rates
decayrate = [0.0] * 7 # Decay rates

# metabolic pool C fractions
fmfaeces = 0.0
fmleaf = 0.0
fmroot = 0.0
faecesn = 0.0


#==============================================================================#