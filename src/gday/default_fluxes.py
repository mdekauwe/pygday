"""
G'DAY default fluxes

Read into the model unless the user changes these at runtime with definitions 
in the .INI file

"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (14.02.2011)"
__email__   = "mdekauwe@gmail.com"

gpp_gCm2 = None
npp_gCm2 = None
gpp = None
npp = None
nep = None
auto_resp = None
hetero_resp = None

# n
nuptake = None
nloss = None
npassive = None # n passive -> active
ngross = None   # N gross mineralisation
nimmob = None   # N immobilisation in SOM
nlittrelease = None # N rel litter = struct + metab
activelossf = None # frac of active C -> CO2

# water fluxes
wue = 0.0
et = 0.0
soil_evap = 0.0
transpiration = 0.0
erain = 0.0
interception = 0.0
runoff = 0.0

# daily C production
cpleaf = None
cproot = None
cpbranch = None
cpstem = None

# daily N production
npleaf = None
nproot = None
npbranch = None
npstemimm = None
npstemmob = None

# dying stuff
deadleaves = None   # Leaf litter C production (t/ha/yr)
deadroots = None    # Root litter C production (t/ha/yr)
deadbranch = None   # Branch litter C production (t/ha/yr)
deadstems = None    # Stem litter C production (t/ha/yr)
deadleafn = None    # Leaf litter N production (t/ha/yr)
deadrootn = None    # Root litter N production (t/ha/yr)
deadbranchn = None  # Branch litter N production (t/ha/yr)
deadstemn = None    # Stem litter N production (t/ha/yr)

# grazing stuff
ceaten = None       # C consumed by grazers (t C/ha/y)
neaten = None       # N consumed by grazers (t C/ha/y)
faecesc = None      # Flux determined by faeces C:N
nurine = None       # Rate of N input to soil in urine (t/ha/y)

#C N root/shoot to struct and metab pools
cresid = [None] * 4  # Cshoot -> surf struct, Croot -> soil sturct, Cshoot -> surf metab, Croot ->surf
nresid = [None] * 4  # Nshoot -> surf struct, Nroot -> soil sturct, Nshoot -> surf metab, Nroot ->surf metab

cstruct = [None] * 4  # Csurf struct -> slow, Csurf struct -> active, Csoil struct -> slow, Csoil struct -> active
nstruct = [None] * 4  # Nsurf struct -> slow, Nsurf struct -> active, Nsoil struct -> slow, Nsoil struct -> active


# C flows to the air
co2_to_air = [None] * 7 

cactive = [None] * 2 # C active -> slow/passive
nactive = [None] * 2 # N active -> slow/passive
passive = None

nslow = [None] * 2 # Nslow -> active and -> passive
cslow = [None] * 2  # Cslow -> active and -> passive
nmetab = [None] * 2 # N surf metab and N soil metab -> active
cmetab = [None] * 2 # C surf metab and C soil metab -> active

cica_avg = None # used in water balance, only when running mate model
nmineralisation = None
apar = None