#!/usr/bin/env python

""" 
Spin up example for Jianyang at the Duke site.
-> Spinup with forest params, fixed NDEP (0.004 t/ha/yr), fixed CO2 (270 ppm).

Site History:
-------------
* Potentially forest until 1700
* grassland until 1983, mowed (+/- annually)
* burnt prior to planting
* At start of experiment aboveground biomass 5.5-11 kg C m-2

Spin-up the model to a steady state. Recycling the met data in batches of a 
50 years, until the SOM, plant and litter C pools cease to change.
"""

import os
import shutil
import sys
import numpy as np
from gday import gday as model
from gday import adjust_gday_param_file as ad

__author__  = "Martin De Kauwe"
__version__ = "1.0 (25.06.2014)"
__email__   = "mdekauwe@gmail.com"


def main(experiment_id, site, SPIN_UP=None):
    
    # dir names
    base_param_name = "base_start"
    base_dir = os.getcwd()
    param_dir = os.path.join(base_dir, "params")
    met_dir = os.path.join(base_dir, "met_data")
    run_dir = os.path.join(base_dir, "outputs")
    
    
    if SPIN_UP == True:
    
        # copy base files to make two new experiment files
        shutil.copy(os.path.join(param_dir, base_param_name + ".cfg"),                
                    os.path.join(param_dir, "%s_%s_model_spinup.cfg" % \
                    (experiment_id, site)))
        
        # Run model to equilibrium assuming forest, growing C pools from 
        # effectively zero
        itag = "%s_%s_model_spinup" % (experiment_id, site)
        otag = "%s_%s_model_spunup" % (experiment_id, site)
        mtag = "%s_met_data_equilibrium_50_yrs.csv" % (site)
        out_fn = itag + "_equilib.out" 
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
                         
                         # control - using Fixed allocation coeffs
                         "alloc_model": '"fixed"',
                         "assim_model": '"mate"',
                         "calc_sw_params": '"true"',   #false=use fwp values, true=derive them
                         "deciduous_model": '"false"',
                         "disturbance": "0",
                         "fixed_stem_nc": '"true"',
                         "fixleafnc": '"false"',
                         "grazing": '"false"',
                         "gs_model": '"medlyn"',
                         "model_optroot": '"false"',
                         "modeljm": '"true"',
                         "nuptake_model": "1",
                         "passiveconst": '"false"',
                         "print_options": '"end"',
                         "ps_pathway": '"c3"',
                         "strfloat": "0",
                         "sw_stress_model": "1",  # Sands and Landsberg 
                         "trans_model": "1",
                         "use_eff_nc": "0",
                         "use_leuning": "0",
                         "water_stress": '"true"',
                         
                         # state - default C:N 25.
                         "age": "0.0",
                         "canht": "17.0", # Canopy height increased from 16m in 2001 to 18m in 2004 at Duke
                         "activesoil": "0.001",
                         "activesoiln": "0.00004",
                         "age": "0.0",
                         "branch": "0.001",
                         "branchn": "0.00004",
                         "cstore": "0.001",
                         "inorgn": "0.00004",
                         "metabsoil": "0.0",
                         "metabsoiln": "0.0",
                         "metabsurf": "0.0",
                         "metabsurfn": "0.0",
                         "nstore": "0.00004",
                         "passivesoil": "0.001",
                         "passivesoiln": "0.0004",
                         "prev_sma": "1.0",
                         "root": "0.001",
                         "root_depth": "-9999.9",
                         "rootn": "0.00004",
                         "sapwood": "0.001",
                         "shoot": "0.001",
                         "shootn": "0.00004",
                         "slowsoil": "0.001",
                         "slowsoiln": "0.00004",
                         "stem": "0.001",
                         "stemn": "0.00004",
                         "stemnimm": "0.00004",
                         "stemnmob": "0.0",
                         "structsoil": "0.001",
                         "structsoiln": "0.00004",
                         "structsurf": "0.001",
                         "structsurfn": "0.00004",
                         
                         
                         # parameters
                         "latitude": "35.9",
                         "intercep_frac": "0.15",
                         "max_intercep_lai": "3.0",
                         "albedo": "0.123",   # modis site avg
                         "finesoil": "0.5",  
                         "slamax": "4.6",     # Drake, 2010, PCE, 33, 1756-1766, fig5 (values are in DM and projected, i.e. one sided, what we need for GDAY). Take value for ~20 years. y=47.2-0.063*x
                         "sla": "4.6",    # Drake, 2010, PCE, 33, 1756-1766, fig5 (values are in DM and projected, i.e. one sided, what we need for GDAY). Take value for ~20 years. y=47.2-0.063*x
                         "slazero": "4.6",    # Drake, 2010, PCE, 33, 1756-1766, fig5 (values are in DM and projected, i.e. one sided, what we need for GDAY). Take value for ~20 years. y=47.2-0.063*x
                         "cfracts": "0.5",
                         "cfracts": "0.5",
                         "lai_cover": "0.5",
                         "jmaxn": "60.0",   # Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "vcmaxn": "30.61", # Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "c_alloc_fmax": "0.25",   
                         "c_alloc_fmin": "0.25",
                         "c_alloc_rmax": "0.05", 
                         "c_alloc_rmin": "0.05", 
                         "c_alloc_bmax": "0.2", 
                         "c_alloc_bmin": "0.2",
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
                         "ncmaxfyoung": "0.06",
                         "ncmaxfold": "0.06",  
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
                         "rateuptake": "3.0",           # set somewhat (very) arbitarly to get an LAI ~ 4.
                         "rateloss": "0.5",
                         "topsoil_depth": "350.0",     # Not needed as I have supplied the root zone water and topsoil water available
                         "rooting_depth": "750.0",     # Not needed as I have supplied the root zone water and topsoil water available
                         "ctheta_topsoil": "0.5",      # Derive based on soil type clay_loam
                         "ntheta_topsoil": "5.0",      # Derive based on soil type clay_loam
                         "ctheta_root": "0.4",         # Derive based on soil type clay
                         "ntheta_root": "3.0",         # Derive based on soil type clay
                         "topsoil_type": '"clay_loam"',
                         "rootsoil_type": '"clay"',
                         "measurement_temp": "25.0",
                         "dz0v_dh": "0.075", # However I have used value from Jarvis, quoted in Jones 1992, pg. 67. Produces a value within the bounds of 3.5-1.1 mol m-2 s-1 Drake, 2010, GCB for canht=17
                         "displace_ratio": "0.78",
                         "g1": "2.74", 
                         "jmaxna": "60.0",  # Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "jmaxnb": "0.0",   # Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "vcmaxna": "30.61",# Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "vcmaxnb": "0.0",  # Original values Belinda had, plus largely match scatter plot fig 7, pg 235 Ellsworth 2011
                         "sapturnover": "0.1",
                         "heighto": "4.826",  
                         "htpower": "0.35",   
                         "height0": "5.0",    
                         "height1": "20.0",   
                         "leafsap0": "8000.0",
                         "leafsap1": "3060.0",  # Duke protocol
                         "branch0": "5.61",   
                         "branch1": "0.346",  
                         "targ_sens": "0.5",  
                         "density": "420.0",  
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        G = model.Gday(cfg_fname, spin_up=True)
        G.spin_up_pools()
        
    
    
if __name__ == "__main__":
    
    experiment_id = "NCEAS"
    site = "DUKE"
    main(experiment_id, site, SPIN_UP=True)
    
    
    
    
    
    
    
    
    
    
    
    
    