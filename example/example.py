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
__version__ = "1.0 (11.02.2014)"
__email__   = "mdekauwe@gmail.com"

def main(experiment_id, site, treatment):
    
    # --- FILE PATHS, DIR NAMES ETC --- #
    base_dir = os.getcwd()
    param_dir = os.path.join(base_dir, "params")
    met_dir = os.path.join(base_dir, "met_data")
    run_dir = os.path.join(base_dir, "runs")
    
    # --- CHANGE PARAM VALUES ON THE FLY --- #
    itag = "%s_%s_model_youngforest_%s" % (experiment_id, site, treatment)
    otag = "%s_%s_model_simulation_%s" % (experiment_id, site, treatment)
    mtag = "%s_met_data_%s.csv" % (site, treatment)
    out_fn = "D1GDAY%s%s.csv" % (site, treatment.upper())
    
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
                         "ps_pathway": '"c3"',
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
    
    # translate output to NCEAS style output
    
    # add this directory to python search path so we can find the scripts!
    sys.path.append(os.path.join(base_dir, "scripts"))
    import translate_GDAY_output_to_NCEAS_format as tr
    tr.translate_output(out_fname, met_fname)
    

if __name__ == "__main__":
    
    
    # Ambient
    experiment_id = "NCEAS"
    site = "DUKE"
    treatment="amb"
    main(experiment_id, site, treatment)