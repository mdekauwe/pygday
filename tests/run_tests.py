#!/usr/bin/env python

""" 
Series of unit-tests for GDAY code to pass at installation time, or
when new code has been changed. 
"""

import os
import sys
import numpy as np
import unittest
from math import exp, sqrt, sin, pi
from gday.mate import MateC3
from gday.file_parser import read_met_forcing
import gday.default_control as control
import gday.default_files as files
import gday.default_params as params
import gday.default_fluxes as fluxes
import gday.default_state as state

__author__  = "Martin De Kauwe"
__version__ = "1.0 (09.012.2014)"
__email__   = "mdekauwe@gmail.com"


    

def testMate():
    
    # set up parameters
    params.alpha_j           = 0.26    
    params.cfracts           = 0.5     
    params.cue               = 0.5     
    params.delsj             = 644.4338
    params.eac               = 79430.0 
    params.eao               = 36380.0 
    params.eag               = 37830.0 
    params.eaj               = 43790.0 
    params.eav               = 51560.0 
    params.edj               = 2e+05   
    params.gamstar25         = 42.75   
    params.jmaxna            = 40.462  
    params.jmaxnb            = 13.691  
    params.kc25              = 404.9   
    params.ko25              = 278400.0
    params.measurement_temp  = 25.0    
    params.oi                = 205000.0
    params.theta             = 0.7     
    params.vcmaxna           = 20.497  
    params.vcmaxnb           = 8.403   
    params.g1                = 4.8 
    params.ctheta_root       = 0.4
    params.ntheta_root       = 3.0
    
    met_fname = "../example/met_data/DUKE_met_data_amb_co2.csv"
    #met_data = read_met_forcing(met_fname, met_header=4)
    
    day = 1
    daylen = 12.0
    
    met_data = {}
    met_data['tam'] = {} 
    met_data['tam'][day] = {} 
    met_data['tpm'] = {} 
    met_data['tpm'][day] = {} 
    met_data['vpd_am'] = {} 
    met_data['vpd_am'][day] = {} 
    met_data['vpd_pm'] = {} 
    met_data['vpd_pm'][day] = {} 
    met_data['co2'] = {} 
    met_data['co2'][day] = {} 
    met_data['par'] = {} 
    met_data['par'][day] = {} 
    
    met_data['tam'][day] = 18.0
    met_data['tpm'][day] = 25.0
    met_data['vpd_am'][day] = 0.8
    met_data['vpd_pm'][day] = 1.5
    met_data["co2"][day] = 380.0
    met_data['par'][day] = 4674600.0
    
    M = MateC3(control, params, state, fluxes, met_data)
    
    state.ncontent = 4.5
    state.lai = 3.0
    state.fipar = (1.0 - exp(-params.kext * state.lai))
    
    M.calculate_photosynthesis(day, daylen)


class GdayTests(unittest.TestCase):
    testMate()
    print "Testing MATE"
    
    def test_total_gpp(self):
        # Values pre-calculated
        gpp_gCm2 = 1.91973402377
        self.assertAlmostEqual(gpp_gCm2, fluxes.gpp_gCm2)
    
    def test_total_npp(self):
        # Values pre-calculated
        npp_gCm2 = 0.959867011885
        self.assertAlmostEqual(npp_gCm2, fluxes.npp_gCm2)
    
    def test_gpp_am(self):    
        # Values pre-calculated
        gpp_am = 1.01899406127
        self.assertAlmostEqual(gpp_am, fluxes.gpp_am)
    
    def test_gpp_pm(self):       
        # Values pre-calculated
        gpp_pm = 0.900739962502
        self.assertAlmostEqual(gpp_pm, fluxes.gpp_pm)
        
        
        
        
if __name__ == "__main__":
    
    unittest.main()
    #testMate()