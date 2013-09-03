"""
G'DAY default module: Initial State

Read into the model unless the user changes these at runtime with definitions
in the .INI file

"""
__author__  = "Martin De Kauwe"
__version__ = "1.0 (14.02.2011)"
__email__   = "mdekauwe@gmail.com"



# Carbon state variables (t/ha)
shoot           = 3.38042             # shoot c
root            = 0.614058            # root c
branch          = 0.0                 # branch c
stem            = 403.635040          # stem c
structsurf      = 10.275649           # surface structural c
metabsurf       = 0.0                 # surface metabolic c
structsoil      = 0.686331            # soil structural c
metabsoil       = 0.0                 # soil metabolic c
activesoil      = 3.607960            # active c
slowsoil        = 61.052441           # slow c
passivesoil     = 29.694              # passive c
cstore          = 0.0                    # C store for deciduous model

# Nitrogen state variables (t/ha)
shootn       = 0.0635                 # shoot n
rootn        = 0.006951               # root n
branchn      = 0.0                    # branch n
sapwood      = 0.0                    # Sapwood: Stem=sapwood+heartwood
stemnimm     = 0.807270               # Immobile stem N (t/ha)
stemnmob     = 0.0                    # Stem N mobile pool
structsurfn  = 0.068504               # surface structural n
metabsurfn   = 0.0                    # surface metabolic n
structsoiln  = 0.004576               # soil structural n
metabsoiln   = 0.0                    # soil metabolic n
activesoiln  = 0.667656               # active n
slowsoiln    = 2.583660               # slow n
passivesoiln = 2.9694                 # passive n
inorgn       = 0.007740               # Inorganic soil N pool - dynamic (t/ha)
stemn        = stemnimm + stemnmob    # Stem N (t/ha)
nstore       = 0.0                    # N store for deciduous model

# Misc state variables
pawater_root  = 200.0              # plant available water - root zone (mm)
pawater_tsoil = 50.0               # plant available water - top soil(mm)
wtfac_root = None
wtfac_tsoil = None
lai = None
sla = None
light_interception = None
ncontent = None
age = 0.0                               #Current stand age (years)

# C allocated fracs - NB these are at the annual timestep for the deciduous model
alleaf = 0.0
alroot = 0.0
albranch = 0.0
alstem = 0.0
delta_sw_store = 0.0

# decid model
c_to_alloc_shoot = 0.0
c_to_alloc_root = 0.0
c_to_alloc_stem = 0.0
n_to_alloc_shoot = 0.0
c_to_alloc_root = 0.0
c_to_alloc_branch = 0.0
n_to_alloc_root = 0.0
n_to_alloc_shoot = 0.0
n_to_alloc_stem = 0.0
n_to_alloc_branch = 0.0
n_to_alloc_stemimm = 0.0
n_to_alloc_stemmob = 0.0
n_to_alloc_branch =0.0

# annual NPP
anpp = 0.0

remaining_days = None
growing_days = None
leaf_out_days = None

root_depth = -9999.9  # rooting depth, Dmax (m)

#============== Not publicly accessible to the user ==========================#
# total plant, soil, litter and system carbon
soilc = 0.0
littercag = 0.0
littercbg = 0.0
litterc = 0.0
plantc =  0.0
totalc = 0.0

# total plant, soil & litter nitrogen
soiln = 0.0
litternag = 0.0
litternbg = 0.0
littern = 0.0
plantn = 0.0
totaln = 0.0

# N:C ratios
rootnc = None
shootnc = None
stemnc = None
branchnc = None
actncslope = None
slowncslope = None
passncslope = None

#==============================================================================#
