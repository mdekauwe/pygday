#!/usr/bin/env python

""" Duke Ambient Simulation for NCEAS FACE experiment 

-> Grow young forest from 1984, varying NDEP/CO2
-> Run Duke experiment (1996-2007), varying NDEP/CO2
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

def main(experiment_id):
    
    # dir names
    base_dir = "/Users/mdekauwe/src/python/pygday/example/"
    param_dir = os.path.join(base_dir, "params")
    met_dir = os.path.join(base_dir, "forcing")
    run_dir = os.path.join(base_dir, "runs")
      
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
                     "wcapac_root": "75.0",      # Oren et al 1998
                     "wcapac_topsoil": "35.0",   # Oren et al 1998
                     "age": "12.0",
                     "print_options": "0",
                     "co2_conc": "0"
                    }
    ad.adjust_param_file(cfg_fname, replace_dict)
    G = model.Gday(cfg_fname)
    G.run_sim()
    
if __name__ == "__main__":
    
    experiment_id = "TEST"
    main(experiment_id)
    
    
    
    
    
    
    
    
    
    
    
    
    