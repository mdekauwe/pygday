#!/usr/bin/env python

""" 
Example script for code testing, the result is not necessarily sensible...this
is meant more for checking changes don't drastically affect code speed.
"""

import os
import shutil
import sys
import subprocess
import timeit
from gday import gday as model
from gday import adjust_gday_param_file as ad

__author__  = "Martin De Kauwe"
__version__ = "1.0 (12.02.2012)"
__email__   = "mdekauwe@gmail.com"

def main(experiment_id, GROW_FOREST=False, RUN_SIM=True):
    
    # dir names
    base_dir = "/Users/mdekauwe/src/python/pygday/example/"
    param_dir = os.path.join(base_dir, "params")
    met_dir = os.path.join(base_dir, "forcing")
    run_dir = os.path.join(base_dir, "runs")
    
    if GROW_FOREST == True:
        # Grow forest from 1984 to start of FACE experiment
        # we are swapping grass params for forest params now
        
        
        itag = experiment_id + "_example"
        otag = experiment_id + "_example_sim"
        mtag = "met_data.gin"
        out_fn = "EXAMPLE.csv"
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        # copy base files to make two new experiment files
        shutil.copy(os.path.join(cfg_fname),                
                    os.path.join(param_dir, itag + "_run" +  ".cfg"))
        cfg_fname = os.path.join(param_dir, itag + "_run" +  ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = { 
                         "out_param_fname": '"%s"' % (out_param_fname),
                         "cfg_fname": '"%s"' % (cfg_fname),
                         "met_fname": '"%s"' % (met_fname),
                         "out_fname": '"%s"' % (out_fname),
                         "age": "0.0",
                         "ageold": "1000",
                         "print_options": "end",
                         "albedo": "0.123", # modis site avg
                         "finesoil": "0.5", 
                         "n_crit": "0.04",
                         "slamax": "4.6",  #Drake, 2010, PCE, 33, 1756-1766, fig5
                         "slainit": "4.6", #Drake, 2010, PCE, 33, 1756-1766, fig5
                         "slazero": "4.6", #Drake, 2010, PCE, 33, 1756-1766, fig5
                         "cfracts": "0.5",
                         "lai_cover": "0.5",
                         "callocf": "0.25",   
                         "callocfz": "0.25",
                         "callocrx": "0.0",  
                         "callocr": "0.05", 
                         "callocrz": "0.05", 
                         "callocb": "0.2", 
                         "callocbz": "0.2",
                         "fretrans": "0.5",  
                         "rretrans": "0.0",  
                         "bretrans": "0.0",  
                         "wretrans": "0.0",
                         "ncwnewz": "0.003",    
                         "ncwnew": "0.003",     
                         "ncwimmz": "0.003",    
                         "ncwimm": "0.003",     
                         "ncbnewz": "0.003",    
                         "ncbnew": "0.003",     
                         "ncrfac": "0.8",     
                         "ncmaxfyoung": "0.04",
                         "ncmaxfold": "0.04",  
                         "ncmaxr": "0.03",     
                         "retransmob": "0.0", 
                         "fdecay": "0.63",   # 19 mth turnover * 1/30, McCarthy, 2007, GCB, 13, 2479-2497  
                         "fdecaydry": "0.63", # 19 mth turnover * 1/30, McCarthy, 2007, GCB, 13, 2479-2497 
                         "rdecay": "0.63",      
                         "rdecaydry": "0.63",   
                         "bdecay": "0.02",      
                         "wdecay": "0.02",      
                         "watdecaydry": "0.0", 
                         "watdecaywet": "0.1", 
                         "ligshoot": "0.25",    
                         "ligroot": "0.25",     
                         "brabove": "0.5",     
                         "branch": "0.01",
                         "shoot": "0.01",
                         "root": "0.01",
                         "stem": "0.01",
                         "metabsurf": "0.0",
                         #"metabsoil": "0.0",
                         "metabsurfn": "0.0",
                         #"metabsoiln": "0.0",
                         "branchn": "0.001",
                         "shootn": "0.001",
                         "rootn": "0.001",
                         "stemn": "2E-03",
                         "stemnimm": "0.001",
                         "stemnmob": "0.001",
                         "wcapac_root": "67.387",      # [mm] (FC-WP)*rooting_depth using derived values and depth from Oren et al 1998 
                         "wcapac_topsoil": "6.045",    # [mm] (FC-WP)*rooting_depth using derived values and depth from Oren et al 1998
                         "fwpmax_tsoil": "0.52",      # this is the saturation point so this is wrong, derive value
                         "fwpmin_tsoil": "0.2",
                         "fwpmax_root": "0.52",       # this is the saturation point so this is wrong, derive value
                         "fwpmin_root": "0.2",
                         "topsoil_type": '"clay_loam"',
                         "rootsoil_type": '"clay"',
                         "calc_sw_params": "1",   #0 uses fwp values, 1= derive them
                         "co2_conc": '"AMB"',
                         "trans_model": "1",
                         "rateuptake": "5.7", 
                         "rateloss": "0.5",
                         "g1": "2.74", # used to be 4.8
                         "deciduous_model": "0",
                         "canht": "17.0",          # Canopy height increased from 16m in 2001 to 18m in 2004 at Duke
                         "dz0v_dh": "0.075",       # However I have used value from Jarvis, quoted in Jones 1992, pg. 67. Produces a value within the bounds of 3.5-1.1 mol m-2 s-1 Drake, 2010, GCB for canht=17
                         "displace_ratio": "0.78",
                         "z0h_z0m": "1.0",         # Assume z0m = z0h, probably a big assumption [as z0h often < z0m.], see comment in code!! 
                         "callocrx": "0.0",
                         "vxfix": "5.0e-9", # N fixation per unit rhizodeposition 5g/1000g -> tonnes 
                         "modeljm": '"true"',
                         "jmaxna": "60.0",  # Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "jmaxnb": "0.0",   # Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "vcmaxna": "30.61",# Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "vcmaxnb": "0.0",  # Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "model_optroot": '"false"',
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        #print cfg_fname
        G = model.Gday(cfg_fname)
        G.run_sim()
    
    if RUN_SIM == True:
        # Run Duke simulation 
        itag = experiment_id + "_example"
        otag = experiment_id + "_example_sim"
        mtag = "met_data.gin"
        out_fn = "EXAMPLE.csv"
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = { 
                         "out_param_fname": '"%s"' % (out_param_fname),
                         "cfg_fname": '"%s"' % (cfg_fname),
                         "met_fname": '"%s"' % (met_fname),
                         "out_fname": '"%s"' % (out_fname),
                         "age": "12.0",
                         "print_options": '"daily"',
                         "co2_conc": '"AMB"',
                         "water_stress": '"true"',
                         "model_optroot": '"false"',
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        
        
        G = model.Gday(cfg_fname)
        G.run_sim()

def profile_main():
    """ profile code """
    import cProfile, pstats
    prof = cProfile.Profile()
    prof = prof.runctx("main()", globals(), locals())
    print "<pre>"
    stats = pstats.Stats(prof)
    stats.sort_stats("cumulative")  # Or cumulative
    stats.print_stats(500)  # 80 = how many to print
    # The rest is optional.
    # stats.print_callees()
    # stats.print_callers()
    print "</pre>"        

if __name__ == "__main__":
    
    import time
    import filecmp
    
    experiment_id = 'TEST'
    N = 1
    times = []
    for i in xrange(N):
        start_time = time.time()
        main(experiment_id, GROW_FOREST=False, RUN_SIM=True) 
        end_time = time.time()
        times.append(end_time - start_time)
    avg_time = sum(times) / float(len(times))
    
    print "\nAverage time (n=%d) : %2.5f seconds\n\n" % (N, avg_time)
    
    
    
    file_dif = filecmp.cmp('runs/EXAMPLE.csv', 'runs/PREVIOUS.csv') 
    if file_dif == False:
        print "Houston we have a problem!!"