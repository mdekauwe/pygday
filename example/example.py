#!/usr/bin/env python

""" 
Example script of how I would run the model, the result is not necessarily 
sensible...

* Note you don't need to change the parameter values this way, though I think
it is preferable.

"""

import os
import sys

# How to import G'day
from gday import gday as model
from gday import adjust_gday_param_file as ad

__author__  = "Martin De Kauwe"
__version__ = "1.0 (10.05.2013)"
__email__   = "mdekauwe@gmail.com"

def main(experiment_id):
    
    # --- FILE PATHS, DIR NAMES ETC --- #
    base_dir = ""
    param_dir = os.path.join(base_dir, "params")
    met_dir = os.path.join(base_dir, "forcing")
    run_dir = os.path.join(base_dir, "runs")
    
    # --- CHANGE PARAM VALUES ON THE FLY --- #
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
                     "canht": "17.0",   # m
                    }
    ad.adjust_param_file(cfg_fname, replace_dict)
    
    # --- RUN THE MODEL --- #
    G = model.Gday(cfg_fname)
    G.run_sim()


if __name__ == "__main__":
    
    experiment_id = 'TEST'
    main(experiment_id) 
   