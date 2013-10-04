#!/usr/bin/env python

""" 
Example script of how I would run the model, the result is not necessarily 
sensible...but essentially the Duke experiment.

* Note you don't need to change the parameter values this way, though I think
it is preferable.

"""

import os
import sys

# How to import G'day
from gday import gday as model
from gday import adjust_gday_param_file as ad

__author__  = "Martin De Kauwe"
__version__ = "1.0 (05.10.2013)"
__email__   = "mdekauwe@gmail.com"

def main(experiment_id):
    
    # --- FILE PATHS, DIR NAMES ETC --- #
    base_dir = os.getcwd()
    param_dir = os.path.join(base_dir, "params")
    met_dir = os.path.join(base_dir, "forcing")
    run_dir = os.path.join(base_dir, "runs")
    
    # --- CHANGE PARAM VALUES ON THE FLY --- #
    itag = experiment_id + "_dk_youngforest"
    otag = experiment_id + "_dk_simulation"
    mtag = "nceas_met_AMB.csv"
    out_fn = "D1GDAYDKAMB.csv"
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
                         "age": "12.0",
                         
                         
                         # control
                         "alloc_model": '"fixed"',
                         "assim_model": '"mate"',
                         "calc_sw_params": "1",   #0 uses fwp values, 1= derive them
                         "deciduous_model": '"false"',
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
                         "sw_stress_model": "1",     # Landsberg
                    }
    ad.adjust_param_file(cfg_fname, replace_dict)
    
    # --- RUN THE MODEL --- #
    G = model.Gday(cfg_fname)
    G.run_sim()


if __name__ == "__main__":
    
    experiment_id = 'NCEAS'
    main(experiment_id) 
   