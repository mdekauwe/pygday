"""
A series of everyday constants, well.

Module defines a series of constant, e.g.
import gday_constants as const
>>>print const.radius_of_earth
>>>6.37122e+6

Refs
====
* McCree, K. J., 1972, Test of current definitions of photosynthetically
active radiation against leaf photosynthesis data. A
gricultural Meteorology, 10, 442-453.

"""
M2_AS_HA = 1E-4
HA_AS_M2 = 1.0 / 1E-4
G_AS_TONNES = 1E-6
TONNES_AS_G = 1.0 / G_AS_TONNES
KG_AS_TONNES = 1E-3
TONNES_AS_KG = 1.0 / KG_AS_TONNES
KG_AS_G = 1E+3
DEG_TO_KELVIN = 273.15
G_TO_KG = 0.001
JOULES_TO_MJ = 1E-6
MJ_TO_JOULES = 1E6
WATT_HR_TO_MJ = 0.0036
HRS_TO_SECS = 3600.0
SECS_TO_DAY = 1 / (60.0 * 60.0 * 24.)
MOL_C_TO_GRAMS_C = 12.0
GRAMS_C_TO_MOL_C = 1.0 / 12.0
MOL_WATER_TO_GRAMS_WATER = 18.0
RGAS = 8.314 # universal gas constant (J mol-1 K-1)
RAD_TO_PAR = 0.48 # McCree, 1972
UMOL_TO_MOL = 1E-6
MOL_TO_UMOL = 1E6
HRS_TO_SECS = 3600.0
DAY_TO_SECS = 60.0 * 60.0 * 24.0
MM_TO_M = 0.001
TONNES_PER_HA_TO_G_M2 = 100.0
NDAYS_IN_YR = 365.25
MOL_TO_MILLIMOLES = 1000.0
KPA_2_PA = 1000.0
PPM_VOL_2_MOL_MOL = 1E-6
MJ_TO_MOL = 4.6
RATIO_DIFF_H2O_TO_CO2 = 1.6



# Converts conductance from units of mol m-2 s-1 to m s-1 at 25 degC
# See Jones Appendix 3 or Diaz et al 2007.
CONV_CONDUCT = 0.0245
