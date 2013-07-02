#!/usr/bin/env python

""" Oak Ridge Ambient Simulation for NCEAS FACE experiment 

-> Run ORNL CO2 experiment: 1998-2008
"""

import os
import shutil
import sys
import subprocess
from gday import gday as model
from gday import adjust_gday_param_file as ad

__author__  = "Martin De Kauwe"
__version__ = "1.0 (31.07.2012)"
__email__   = "mdekauwe@gmail.com"

def main(experiment_id, GROW_FOREST=True, RUN_SIM=True):
    
    # dir names
    base_dir = "/Users/mdekauwe/research/NCEAS_face/GDAY_ornl_simulation"
    param_dir = os.path.join(base_dir, "params")
    met_dir = os.path.join(base_dir, "forcing")
    run_dir = os.path.join(base_dir, "runs")
    
    if GROW_FOREST == True:
        # copy spin-upbase files to make two new experiment files
        shutil.copy(os.path.join(param_dir, experiment_id + "_or_model_indust.cfg"),
                    os.path.join(param_dir, experiment_id + "_or_model_indust_adj.cfg"))
    
        # Grow forest from 1984 to start of FACE experiment
        # we are swapping grass params for forest params now
        itag = experiment_id + "_or_model_indust_adj"
        otag = experiment_id + "_or_youngforest_amb"
        mtag = "met_expgrw.gin"
        out_fn = itag + "_youngforest.out"
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        
        replace_dict = { 
                        # files
                         "out_param_fname": '"%s"' % (out_param_fname),
                         "cfg_fname": '"%s"' % (cfg_fname),
                         "met_fname": '"%s"' % (met_fname),
                         "out_fname": '"%s"' % (out_fname),
                         
                         # parameters
                         "canht": "12.4", # Canopy height increased from 16m in 2001 to 18m in 2004 at Duke
                         "displace_ratio": "0.67",
                         "z0h_z0m": "1.0",         # Assume z0m = z0h, probably a big assumption [as z0h often < z0m.], see comment in code!! 
                         "dz0v_dh": "0.1",       # However I have used value from Jarvis, quoted in Jones 1992, pg. 67. Produces a value within the bounds of 3.5-1.1 mol m-2 s-1 Drake, 2010, GCB for canht=17
                         "ncfmin": "0.0151",     # Minimum of leaf litter values from each year, divided by an average of 46.3% C in leaf litter, calculated from Silvia's worksheet.
                         "sla": "11.3",
                         "slainit": "11.3",
                         "slamax": "11.3",
                         "slazero": "11.3",
                         "wdecay": "0.011",
                         "wdecaydry": "0.011",
                         "fdecay": "1.277",         # using same as roots, roughly 9 months
                         "fdecaydry": "1.277",      # using same as roots, roughly 9 months
                         "rdecay": "1.277",         # from silvia's spreadsheet
                         "rdecaydry": "1.277",      # from silvia's spreadsheet
                         "bdecay": "0.011",
                         "brabove": "0.5",
                         "finesoil": "0.5",
                         "fretrans": "0.5",
                         "wcapac_root": "332.0",       # [mm] (FC (m3/m-3)-WP (m3/m-3)) * rooting_depth (mm) using derived values and depth from protocol, 2000 mm
                         "wcapac_topsoil": "33.2",      # [mm] (FC (m3/m-3)-WP (m3/m-3)) * rooting_depth (mm) using derived values and depth from protocol, assuming 200 mm top soil following Corbeels 2005a.
                         "fwpmax_tsoil": "0.366",      # protocol
                         "fwpmin_tsoil": "0.2",        # protocol
                         "fwpmax_root": "0.366",       # protocol
                         "fwpmin_root": "0.2",         # protocol
                         "topsoil_type": '"silty_clay_loam"',
                         "rootsoil_type": '"silty_clay_loam"',
                         "lai_cover": "0.5",
                         "n_crit": "0.03556",
                         "ncrfac": "0.8",
                         "ncwnewz": "0.003",          #New stem ring N:C at zero leaf N:C (mobile)
                         "ncwnew": "0.003",           #New stem ring N:C at critical leaf N:C (mob)
                         "ncwimmz": "0.003",          #Immobile stem N C at zero leaf N C
                         "ncwimm": "0.003",           #Immobile stem N C at critical leaf N C
                         "ncbnewz": "0.003",          #new branch N C at zero leaf N C
                         "ncbnew": "0.003",           #new branch N C at critical leaf N C
                         "ncmaxfyoung": "0.04",       #max N:C ratio of foliage in young stand, if the same as old=no effect
                         "ncmaxfold": "0.04",         #max N:C ratio of foliage in old stand, if the same as young=no effect
                         "ncmaxr": "0.03",            #max N:C ratio of roots
                         "eav": "58520.0",    # silva excel spreadsheet
                         "eaj": "38670.0",    # silva excel spreadsheet
                         "edj": "2e+05",      # silva excel spreadsheet
                         "delsj": "638.1",    # silva excel spreadsheet
                         "theta": "0.95",     # silva excel spreadsheet
                         "latitude": "35.9",  # silva excel spreadsheet
                         "g1": "4.4",         # used to be"7.5"
                         #"rateuptake": "5.7",
                         #"rateloss": "0.8",
                         "rateuptake": "2.23727164",
                         "rateloss": "0.5",
                         "ligshoot": "0.279",
                         "ligroot": "0.523",
                         "cfracts": "0.467",
                         "ligshoot": "0.279",
                         "jmaxna": "40.462",
                         "jmaxnb": "13.691",
                         "vcmaxna": "20.497",
                         "vcmaxnb": "8.403",
                         "d0x": "0.35",
                         "r0": "0.1325",
                         "top_soil_depth": "0.3",
                         "previous_ncd": "19.0",
                         
                         "callocf": "0.26",           #allocation to leaves at leaf n_crit, number from silva xls spreadsheet.
                         "callocfz": "0.26",          #allocation to leaves at zero leaf n/c, number from silva xls spreadsheet.
                         "callocr": "0.11",           #allocation to roots at root n_crit, number from silva xls spreadsheet.
                         "callocrz": "0.11",          #allocation to roots at root leaf n/c, number from silva xls spreadsheet.
                         "callocb": "0.06",           #allocation to branches at branch  n_crit, number from silva xls spreadsheet. Fraction of wood in branches is 0.0941 * fraction of C in wood (0.63)
                         "callocbz": "0.06",          #allocation to branches at zero branch n/c, number from silva xls spreadsheet.
                         
                         
                         
                         # state
                         "age": "0.0",
                         "cstore": "0.01",
                         "nstore": "0.01",
                         "branch": "0.01",
                         "branchn": "0.001",
                         "shoot": "0.01",
                         "shootn": "0.001",
                         "root": "0.01",
                         "rootn": "0.001",
                         "stem": "0.01",
                         "stemn": "2E-03",
                         "stemnimm": "0.001",
                         "stemnmob": "0.001",
                         "metabsurf": "0.0",
                         "metabsurfn": "0.0",
                         
                         
                         # control
                         "assim_model": '"mate"',
                         "calc_sw_params": "0",   #0 uses fwp values, 1= derive them
                         "deciduous_model": '"true"',
                         "fixed_stem_nc": "1",
                         "fixleafnc": '"false"',
                         "grazing": '"false"',
                         "model_optroot": '"false"',
                         "modeljm": '"true"',
                         "nuptake_model": "1",
                         "passiveconst": '"false"',
                         "print_options": '"end"',
                         "strfloat": "0",
                         "trans_model": "1",
                         "use_eff_nc": "0",
                         "use_leuning": "0",
                         "water_stress": '"true"'
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        G = model.Gday(cfg_fname)
        G.run_sim()
        
    if RUN_SIM == True:
        # Run ORNL simulation 
        itag = experiment_id + "_or_youngforest_amb"
        otag = experiment_id + "_or_simulation"
        mtag = "met_exprun_AMB.gin"
        out_fn = "D1GDAYORAMB.csv"
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        
        
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = { 
                         # files
                         "out_param_fname": '"%s"' % (out_param_fname),
                         "cfg_fname": '"%s"' % (cfg_fname),
                         "met_fname": '"%s"' % (met_fname),
                         "out_fname": '"%s"' % (out_fname),
                         
                         # state 
                         "age": "14.0",
                          
                             
                         # parameters
                         "previous_ncd": "35.0",
                         
                         # control
                         "assim_model": '"mate"',
                         "calc_sw_params": "0",   #0 uses fwp values, 1= derive them
                         "deciduous_model": '"true"',
                         "fixed_stem_nc": "1",
                         "fixleafnc": '"false"',
                         "grazing": '"false"',
                         "model_optroot": '"false"',
                         "modeljm": '"true"',
                         "nuptake_model": "1",
                         "passiveconst": '"false"',
                         "print_options": '"daily"',
                         "strfloat": "0",
                         "trans_model": "1",
                         "use_eff_nc": "0",
                         "use_leuning": "0",
                         "water_stress": '"true"',
                         
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        #print cfg_fname
        G = model.Gday(cfg_fname)
        G.run_sim()
        
        # translate output to NCEAS style output
        
        # add this directory to python search path so we can find the scripts!
        sys.path.append(os.path.join(base_dir, "scripts"))
        import translate_GDAYoutput_to_NCEAS_format as tr
        tr.translate_output(out_fname, met_fname, "ornl_ambient_simulation")
        

        
    
if __name__ == "__main__":
    
    experiment_id = "NCEAS"
    main(experiment_id)
    
    
    
    
    
    
    
    
    
    
    
    
    