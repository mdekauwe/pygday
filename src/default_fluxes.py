"""
G'DAY default fluxes

Read into the model unless the user changes these at runtime with definitions 
in the .INI file

"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (14.02.2011)"
__email__   = "mdekauwe@gmail.com"

# Carbon fluxes
gpp_gCm2 = None
npp_gCm2 = None
gpp_am_pm = [0.0, 0.0]
gpp = None
npp = None
nep = None
auto_resp = None
hetero_resp = None
retrans = None

# Nitrogen fluxes
nuptake = None
nloss = None
npassive = None # n passive -> active
ngross = None   # N gross mineralisation
nimmob = None   # N immobilisation in SOM
nlittrelease = None # N rel litter = struct + metab
activelossf = None # frac of active C -> CO2
nmineralisation = None

# water fluxes
wue = 0.0
et = 0.0
soil_evap = 0.0
transpiration = 0.0
erain = 0.0
interception = 0.0
runoff = 0.0
gs_mol_m2_sec = 0.0
ga_mol_m2_sec = 0.0
omega = 0.0

# daily C production
cpleaf = None
cproot = None
cpcroot = None
cpbranch = None
cpstem = None

# daily N production
npleaf = None
nproot = None
npcroot = None
npbranch = None
npstemimm = None
npstemmob = None

# dying stuff
deadleaves = None   # Leaf litter C production (t/ha/yr)
deadroots = None    # Root litter C production (t/ha/yr)
deadcroots = None   # Coarse Root litter C production (t/ha/yr)
deadbranch = None   # Branch litter C production (t/ha/yr)
deadstems = None    # Stem litter C production (t/ha/yr)
deadleafn = None    # Leaf litter N production (t/ha/yr)
deadrootn = None    # Root litter N production (t/ha/yr)
deadcrootn = None   # Root litter N production (t/ha/yr)
deadbranchn = None  # Branch litter N production (t/ha/yr)
deadstemn = None    # Stem litter N production (t/ha/yr)

# grazing stuff
ceaten = None       # C consumed by grazers (t C/ha/y)
neaten = None       # N consumed by grazers (t C/ha/y)
faecesc = None      # Flux determined by faeces C:N
nurine = None       # Rate of N input to soil in urine (t/ha/y)

leafretransn = None

# C&N Surface litter
surf_struct_litter = None
surf_metab_litter = None
n_surf_struct_litter = None
n_surf_metab_litter = None

# C&N Root Litter
soil_struct_litter = None
soil_metab_litter = None
n_soil_struct_litter = None
n_soil_metab_litter = None

# C&N litter fluxes to slow pool
surf_struct_to_slow = None
soil_struct_to_slow = None
n_surf_struct_to_slow = None
n_soil_struct_to_slow = None

# C&N litter fluxes to active pool
surf_struct_to_active = None
soil_struct_to_active = None
n_surf_struct_to_active = None
n_soil_struct_to_active = None

# Metabolic fluxes to active pool
surf_metab_to_active = None
soil_metab_to_active = None
n_surf_metab_to_active = None
n_surf_metab_to_active = None

# C fluxes out of active pool
active_to_slow = None
active_to_passive = None
n_active_to_slow = None
n_active_to_passive = None

# C&N fluxes from slow to active pool
slow_to_active = None
slow_to_passive = None
n_slow_to_active = None
n_slow_to_passive = None

# C&N fluxes from passive to active pool
passive_to_active = None
n_passive_to_active = None

# C & N source fluxes from the active, slow and passive pools
c_into_active = None  
c_into_slow = None    
c_into_passive = None 

# Microbial respiration -> CO2
co2_to_air = [None] * 7 


# priming fluxes
root_exc = 0.0
root_exn = 0.0
co2_released_exud = 0.0
factive = 0.0
rtslow = 0.0
rexc_cue = 0.0

# C allocated fracs - NB these are at the annual timestep for the deciduous model
alleaf = 0.0
alroot = 0.0
albranch = 0.0
alstem = 0.0

# Misc stuff
cica_avg = None # used in water balance, only when running mate model
apar = None
rabove = 0.0
microbial_resp = 0.0
tfac_soil_decomp = None 
co2_rel_from_surf_struct_litter = None
co2_rel_from_soil_struct_litter = None
co2_rel_from_surf_metab_litter = None
co2_rel_from_soil_metab_litter = None
co2_rel_from_active_pool = None
co2_rel_from_slow_pool = None
co2_rel_from_passive_pool = None
