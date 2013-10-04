#!/usr/bin/env python

""" Translate GDAY output file 

Match the NCEAS format and while we are at it carry out unit conversion so that 
we matched required standard. Data should be comma-delimited
"""

import os
import numpy as np
import csv
import sys
import pandas as pd
import datetime as dt
from cStringIO import StringIO
import shutil

__author__  = "Martin De Kauwe"
__version__ = "1.0 (06.04.2011)"
__email__   = "mdekauwe@gmail.com"

def date_converter(*args): 
    return dt.datetime.strptime(str(int(float(args[0]))) + " " +\
                                str(int(float(args[1]))), '%Y %j')

def translate_output(infname, met_fname):
    outdir = "runs"
    UNDEF = -9999.
    units = setup_units()
    variable, variable_names = setup_varnames()

    # load met stuff, i.e. the stuff needed for NCEAS output that G'day 
    # does not output
    envir = load_met_input_data(met_fname)
    
    # load the rest of the g'day output
    gday = load_gday_output(infname)
    
    # merge dictionaries to ease output
    data_dict = dict(envir, **gday)
    
    #f = open("D1GDAYDKAMB_MGDK.txt", 'w')
    ofname = os.path.join(outdir, "temp.nceas")
    f = open(ofname, "w")
    # write output in csv format
    writer = csv.writer(f, dialect=csv.excel)
    
    
    
    # write header for csv file
    writer.writerow(variable)
    writer.writerow(units)
    writer.writerow(variable_names)
    for i in xrange(len(gday['DOY'])):
        writer.writerow([("%.8f" % (float(data_dict[k][i])) \
                         if data_dict.has_key(k) else UNDEF)
                         for k in variable_names])
    
    # Need to replace the temp file with the infname which is actually
    # the filename we want to use
    shutil.move(ofname, infname)
       
    
    
    
def remove_comments_from_header(fname):
    """ I have made files with comments which means the headings can't be 
    parsed to get dictionary headers for pandas! Solution is to remove these
    comments first """
    s = StringIO()
    with open(fname) as f:
        for line in f:
            if '#' in line:
                line = line.replace("#", "").lstrip(' ')
            s.write(line)
    s.seek(0) # "rewind" to the beginning of the StringIO object
    
    return s


def load_met_input_data(fname):
    MJ_TO_MOL = 4.6
    SW_TO_PAR = 0.48
    DAYS_TO_HRS = 24.0
    UMOL_TO_MOL = 1E-6
    tonnes_per_ha_to_g_m2 = 0.01
    
    s = remove_comments_from_header(fname)
    met_data = pd.read_csv(s, parse_dates=[[0,1]], skiprows=4, index_col=0, 
                           sep=",", keep_date_col=True, 
                           date_parser=date_converter)
    
    precip = met_data["rain"]
    #par = met_data[:,1] * MJ_TO_MOL * SW_TO_PAR
    par = met_data["par"] * UMOL_TO_MOL
    air_temp = met_data["tair"]
    soil_temp = met_data["tsoil"]
    vpd = met_data["vpd_avg"] 
    co2 = met_data["co2"]
    ndep = met_data["ndep"] / tonnes_per_ha_to_g_m2
    
    return {'CO2': co2, 'PPT':precip, 'PAR':par, 'AT':air_temp, 'ST':soil_temp, 
            'VPD':vpd, 'NDEP':ndep}
    
    
    
def load_gday_output(fname):
    TOPSOIL=42.6712851218
    SW_RAD_TO_PAR = 2.0
    UNDEF = -9999.
    tonnes_per_ha_to_g_m2 = 0.01
    yr_to_day = 365.25
    
    s = remove_comments_from_header(fname)
    out = pd.read_csv(s, parse_dates=[[0,1]], skiprows=0, index_col=0, 
                      sep=",", keep_date_col=True, date_parser=date_converter)
    
    
    year = out["year"]
    doy = out["doy"]
    
    
    # state outputs
    
    pawater_topsoil = (out["pawater_topsoil"] / TOPSOIL) * 100.0
    shoot = out["shoot"] / tonnes_per_ha_to_g_m2
    stem = out["stem"] / tonnes_per_ha_to_g_m2
    branch = out["branch"] / tonnes_per_ha_to_g_m2
    fine_root = out["root"] / tonnes_per_ha_to_g_m2
    litterc = out["litterc"] / tonnes_per_ha_to_g_m2
    littercag = out["littercag"] / tonnes_per_ha_to_g_m2
    littercbg = out["littercbg"] / tonnes_per_ha_to_g_m2
    soilc = out["soilc"] / tonnes_per_ha_to_g_m2
    lai = out["lai"]
    shootn = out["shootn"] / tonnes_per_ha_to_g_m2
    stemn = out["stemn"] / tonnes_per_ha_to_g_m2
    branchn = out["branchn"] / tonnes_per_ha_to_g_m2
    rootn = out["rootn"] / tonnes_per_ha_to_g_m2
    litternag = out["litternag"] / tonnes_per_ha_to_g_m2
    litternbg = out["litternbg"] / tonnes_per_ha_to_g_m2
    nsoil = out["soiln"] / tonnes_per_ha_to_g_m2
    inorgn = out["inorgn"] / tonnes_per_ha_to_g_m2
    tnc = out["cstore"] / tonnes_per_ha_to_g_m2
    nstorage = out["nstore"] / tonnes_per_ha_to_g_m2

    
    # fluxes outputs
    nep = out["nep"] / tonnes_per_ha_to_g_m2
    gpp = out["gpp"] / tonnes_per_ha_to_g_m2
    npp = out["npp"] / tonnes_per_ha_to_g_m2
    rh = out["hetero_resp"] / tonnes_per_ha_to_g_m2
    ra = out["auto_resp"] / tonnes_per_ha_to_g_m2
    et = out["et"] # mm of water' are same value as kg/m2
    trans = out["transpiration"] # mm of water' are same value as kg/m2
    soil_evap = out["soil_evap"] # mm of water' are same value as kg/m2
    can_evap = out["interception"] # mm of water' are same value as kg/m2
    runoff = out["runoff"] # mm of water' are same value as kg/m2
    gl = out["cpleaf"] / tonnes_per_ha_to_g_m2
    # gw summed from cpstem and cpbranch below
    cpstem = out["cpstem"] / tonnes_per_ha_to_g_m2
    cpbranch = out["cpbranch"] / tonnes_per_ha_to_g_m2
    gr = out["cproot"] / tonnes_per_ha_to_g_m2
    deadleaves = out["deadleaves"] / tonnes_per_ha_to_g_m2
    deadroots = out["deadroots"] / tonnes_per_ha_to_g_m2
    deadbranch = out["deadbranch"] / tonnes_per_ha_to_g_m2
    deadstems = out["deadstems"] / tonnes_per_ha_to_g_m2
    deadleafn = out["deadleafn"] / tonnes_per_ha_to_g_m2
    deadbranchn = out["deadbranchn"] / tonnes_per_ha_to_g_m2
    deadstemn = out["deadstemn"] / tonnes_per_ha_to_g_m2
    deadrootn = out["deadrootn"] / tonnes_per_ha_to_g_m2
    nup = out["nuptake"] / tonnes_per_ha_to_g_m2
    ngross = out["ngross"] / tonnes_per_ha_to_g_m2
    nmin = out["nmineralisation"] / tonnes_per_ha_to_g_m2
    npleaf = out["npleaf"] / tonnes_per_ha_to_g_m2
    nproot = out["nproot"] / tonnes_per_ha_to_g_m2
    npstemimm = out["npstemimm"] / tonnes_per_ha_to_g_m2
    npstemmob = out["npstemmob"] / tonnes_per_ha_to_g_m2
    npbranch = out["npbranch"] / tonnes_per_ha_to_g_m2
    apar = out["apar"] / SW_RAD_TO_PAR
    gcd = out["gs_mol_m2_sec"]
    ga = out["ga_mol_m2_sec"]
    nleach = out["nloss"] / tonnes_per_ha_to_g_m2

    
    # Misc stuff we don't output
    
    drainage = [UNDEF] * len(doy)
    rleaf = [UNDEF] * len(doy)
    rwood = [UNDEF] * len(doy)
    rgrow = [UNDEF] * len(doy)
    rsoil = [UNDEF] * len(doy)
    cex = [UNDEF] * len(doy)
    cvoc = [UNDEF] * len(doy)
    lh = [UNDEF] * len(doy)
    sh = [UNDEF] * len(doy)
    coarse_root = [UNDEF] * len(doy)
    ccoarse_lit = [UNDEF] * len(doy)
    gcr = [UNDEF] * len(doy)
    coarse_rootn = [UNDEF] * len(doy)
    ndw = [UNDEF] * len(doy)
    npoolo = [UNDEF] * len(doy)
    nfix = [UNDEF] * len(doy)
    nvol = [UNDEF] * len(doy)
    ngcr = [UNDEF] * len(doy)
    
    # Misc calcs from fluxes/state
    lma = shoot / lai
    ncon = shootn / shoot 
    recosys = rh + ra
    cw = stem + branch
    gw = cpstem + cpbranch
    cwn = stemn + branchn
    cwin = deadstems + deadbranch
    ndeadwood = deadbranchn + deadstemn
    nwood_growth = npstemimm + npstemmob + npbranch
    
    
    
    return {'YEAR':year, 'DOY':doy, 'SW':pawater_topsoil, 
            'NEP':nep, 'GPP':gpp, 'NPP':npp, 'CEX':cex, 'CVOC':cvoc, 
            'RECO':recosys, 'RAUTO':ra, 'RLEAF':rleaf, 'RWOOD':rwood, 
            'RGROW':rgrow, 'RHET':rh, 'RSOIL':rsoil, 'ET':et, 'T':trans,
            'ES':soil_evap, 'EC':can_evap, 'RO':runoff, 'DRAIN':drainage,
            'LE':lh, 'SH':sh, 'CL':shoot, 'CW':cw, 'CCR':coarse_root,
            'CFR':fine_root, 'TNC':tnc, 'CFLIT':litterc, 'CFLITA':littercag,
            'CFLITB':littercbg, 'CCLITB':ccoarse_lit, 'CSOIL':soilc,
            'GL':gl, 'GW':gw, 'GCR':gcr, 'GR':gr, 'CLLFALL':deadleaves, 
            'CRLIN':deadroots, 'CWIN':cwin, 'LAI':lai, 'LMA':lma, 'NCON':ncon,
            'NCAN':shootn, 'NWOOD':cwn, 'NCR':coarse_rootn, 'NFR':rootn,
            'NSTOR':nstorage, 'NLIT':litternag, 'NRLIT':litternbg, 'NDW':ndw, 
            'NSOIL':nsoil, 'NPOOLM':inorgn, 'NPOOLO':npoolo, 'NFIX':nfix,
            'NLITIN':deadleafn, 'NWLIN':ndeadwood, 'NRLIN':deadrootn, 'NUP':nup,
            'NGMIN':ngross, 'NMIN':nmin, 'NVOL': nvol, 'NLEACH':nleach, 
            'NGL':npleaf, 'NGW':nwood_growth, 'NGCR':ngcr, 'NGR':nproot, 
            'APARd':apar, 'GCd':gcd, 'Ga':ga}
  
def setup_units():    
    units = ['--','--','Mean ppm', 'PPT', 'mol m-2', 'Mean DegC', 'Mean DegC', 
                'kPa h', 'Mean %', 'gN m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'kgH2O m-2 d-1', 
                'kgH2O m-2 d-1', 'kgH2O m-2 d-1', 'kgH2O m-2 d-1', 
                'kgH2O m-2 d-1', 'kgH2O m-2 d-1', 'MJ m-2', 'MJ m-2', 
                'gC m-2', 'gC m-2', 'gC m-2', 'gC m-2', 'gC m-2', 'gC m-2', 
                'gC m-2', 'gC m-2', 'gC m-2', 'gC m-2 0 to 30 cm', 
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'm2 m-2', 'gC m-2',
                'gN gd.m.-1', 'gN m-2', 'gN m-2', 'gN m-2', 'gN m-2', 'gN m-2', 
                'gN m-2', 'gN m-2', 'gN m-2', 'gN m-2 0 to 30 cm', 
                'gN m-2 0 to 30 cm', 'gN m-2 0 to 30 cm', 'gN m-2 d-1', 
                'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1', 
                'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1', 
                'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1', 
                'MJ m-2 d-1', 'mol H2O m-2 s-1', 'mol H2O m-2 s-1']
    return units 
    
def setup_varnames():
    variable = ['Year', 'Day of the year', 'CO2', 'Precipitation', 'PAR', 
                    'Air temp canopy', 'Soil temp 10 cm', 'Vapour Pres Def',
                    'Surface soil water', 'N deposition', 'Net Eco Prod', 
                    'Gross Prim Prod', 'Net Prim Prod', 'C exudation', 
                    'C VOC Flux', 'Resp ecosystem', 'Resp autotrophic',
                    'Resp leaves (maint)', 'Resp Wood (maint)', 
                    'Resp Fine Root (maint)', 'Resp growth', 
                    'Resp heterotrophic', 'Resp from soil', 
                    'Evapotranspiration', 'Transpiration', 'Soil Evaporation',
                    'Canopy evaporation', 'Runoff', 'Drainage', 'Latent Energy',
                    'Sensible Heat', 'C Leaf Mass', 'C Wood Mass', 
                    'C Coarse Root mass', 'C Fine Root mass', 
                    'C Storage as TNC', 'C Fine Litter Total', 
                    'C Fine Litter above', 'C Fine Litter below', 
                    'C Coarse Litter', 'C Soil', 'C Leaf growth', 
                    'C Wood growth', 'C Coarse Root growth', 
                    'C Fine Root growth', 'C Leaf Litterfall', 
                    'C Root litter inputs', 'C Wood/branch inputs', 
                    'LAI projected', 'Leaf gC/leaf area', 'N Conc Leaves',
                    'N Mass Leaves', 'N Mass Wood', 'N Mass Coarse Roots',
                    'N Mass Fine Roots', 'N storage', 'N litter aboveground',
                    'N litter belowground', 'N Dead wood', 'N Soil Total',
                    'N in Mineral form', 'N in Organic form', 'N fixation',
                    'N Leaf Litterfall', 'N Wood/brch litterfall', 
                    'N Root litter input', 'N Biomass Uptake', 
                    'N Gross Mineralization', 'N Net mineralization', 
                    'N Volatilization', 'N Leaching', 'N Leaf growth',
                    'N Wood growth', 'N CR growth', 'N Fine Root growth',
                    'Aborbed PAR', 'Average daytime canopy conductance',
                    'Average daytime aerodynamic conductance']
    
    variable_names = ['YEAR', 'DOY', 'CO2', 'PPT', 'PAR', 'AT', 'ST', 'VPD', 
                        'SW', 'NDEP', 'NEP', 'GPP', 'NPP', 'CEX', 'CVOC', 
                        'RECO', 'RAUTO', 'RLEAF', 'RWOOD', 'RROOT', 'RGROW', 
                        'RHET', 'RSOIL',
                        'ET', 'T', 'ES', 'EC', 'RO', 'DRAIN', 'LE', 'SH',
                        'CL', 'CW', 'CCR', 'CFR', 'TNC', 'CFLIT', 'CFLITA',
                        'CFLITB', 'CCLITB', 'CSOIL', 'GL', 'GW', 'GCR', 'GR',
                        'CLLFALL', 'CRLIN', 'CWIN', 'LAI', 'LMA', 'NCON', 
                        'NCAN', 'NWOOD', 'NCR', 'NFR', 'NSTOR', 'NLIT', 
                        'NRLIT', 'NDW', 'NSOIL', 'NPOOLM', 'NPOOLO', 'NFIX',
                        'NLITIN', 'NWLIN', 'NRLIN', 'NUP', 'NGMIN', 'NMIN',
                        'NVOL', 'NLEACH', 'NGL', 'NGW', 'NGCR', 'NGR', 'APARd',
                        'GCd', 'Ga']
    
    return variable, variable_names
    
        
if __name__ == "__main__":
    fname = "dk_fixco2_fixndep_forest_equilib.out"
    met_fname = "duke_equilibrium_metdata_fixndep_0.004_fixco2_270.gin"
    translate_output(fname, met_fname)
    
    
    
    
    
    
    
    
    
    
    
    