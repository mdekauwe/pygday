#!/usr/bin/env python

""" run GDAY, plot LAI, NPP, transpiration """

import os
import sys
import numpy as np
import pandas as pd
import datetime as dt
import unittest
import matplotlib.pyplot as plt

__author__  = "Martin De Kauwe"
__version__ = "1.0 (11.02.2014)"
__email__   = "mdekauwe@gmail.com"

# Run example simulation
site = "DUKE"
treatment="amb"
os.chdir("../example/")
os.system("example.py")
os.chdir("../tests/")

def date_converter(*args): 
    return dt.datetime.strptime(str(int(float(args[0]))) + " " +\
                                str(int(float(args[1]))), '%Y %j')

def read_data(fname):
    df = pd.read_csv(fname, parse_dates=[[0,1]], index_col=0, sep=",", 
                     keep_date_col=True, date_parser=date_converter, 
                     na_values=["-9999"], skiprows=2)
    return df
    
# load data
df = read_data("../example/runs/D1GDAY%s%s.csv" % (site, treatment.upper()))

def evap_trans():
    print "ET = T+ES+EC"
    ET = df.groupby("YEAR").ET.sum()
    T = df.groupby("YEAR").T.sum()
    ES = df.groupby("YEAR").ES.sum()
    EC = df.groupby("YEAR").EC.sum()

    X = ET
    Y = T+ES+EC
    
    return X.values, Y.values
       

def soil_water():
    print "deltaSW = PPT - ET - RO - DRAIN"
    yrs = np.unique(df["YEAR"])
    
    deltaSW = np.zeros(0)
    for yr in yrs:
        yrs_data = df[df["YEAR"] == yr]
        change = (yrs_data["SW"][-1] - 
                  yrs_data["SW"][0])
        deltaSW = np.append(deltaSW, change)
    
    PPT = df.groupby("YEAR").PPT.sum()
    ET = df.groupby("YEAR").ET.sum()
    RO = df.groupby("YEAR").RO.sum()
    DRAIN = df.groupby("YEAR").DRAIN.sum()
    
    X = deltaSW
    
    Y = PPT - ET
    if np.all(np.isnan(RO)) == False:
        Y -= RO
    if np.all(np.isnan(DRAIN)) == False:
        Y -= DRAIN
    
    return X, Y.values

def NPP_resp():
    print "NPP = GPP - Rauto"
    yrs = np.unique(df["YEAR"])
            
    NPP = df.groupby("YEAR").NPP.sum()
    GPP = df.groupby("YEAR").GPP.sum()
    RAUTO = df.groupby("YEAR").RAUTO.sum()
    
    X = NPP
    Y = GPP-RAUTO    
    
    return X.values, Y.values

def NPP_growth():
    print "NPP = GL + GW + GCR + GR + GREPR + change in TNC + CVOC"
    yrs = np.unique(df["YEAR"])
    deltaTNC = np.zeros(0)
    for yr in yrs:
        yrs_data = df[df["YEAR"] == yr]
        change = (yrs_data["TNC"][-1] - 
                  yrs_data["TNC"][0])
        deltaTNC = np.append(deltaTNC, change)
    
    NPP = df.groupby("YEAR").NPP.sum()
    GL = df.groupby("YEAR").GL.sum()
    GW = df.groupby("YEAR").GW.sum()
    GCR = df.groupby("YEAR").GCR.sum()
    GR = df.groupby("YEAR").GR.sum()
    GREPR = df.groupby("YEAR").GREPR.sum()
    CVOC = df.groupby("YEAR").CVOC.sum()
            
    X = NPP
    
    Y = GL + GW
    if np.all(np.isnan(GCR)) == False:
        Y += GCR
    if np.all(np.isnan(GR)) == False:
        Y += GR
    if np.all(np.isnan(deltaTNC)) == False:
        Y += deltaTNC
    if np.all(np.isnan(CVOC)) == False:
        Y += CVOC
    if np.all(np.isnan(GREPR)) == False:
        Y += GREPR
    
    return X.values, Y.values


def NEP():
    print "NEP = GPP - Reco"
    yrs = np.unique(df["YEAR"])
            
    NEP = df.groupby("YEAR").NEP.sum()
    GPP = df.groupby("YEAR").GPP.sum()
    RECO = df.groupby("YEAR").RECO.sum()

    X = NEP
    Y = GPP - RECO
    
    return X.values, Y.values        
 
 
 
def cf_stock():
    print "deltaCL = GL - CLLFALL"
    yrs = np.unique(df["YEAR"])
    deltaCL = np.zeros(0)
    for yr in yrs:
        yrs_data = df[df["YEAR"] == yr]
        change = (yrs_data["CL"][-1] - 
                  yrs_data["CL"][0])
        deltaCL = np.append(deltaCL, change)
            
    GL = df.groupby("YEAR").GL.sum()
    CLLFALL = df.groupby("YEAR").CLLFALL.sum()
            
    X = deltaCL
    Y = GL - CLLFALL
    return X, Y.values   
          
#
# - Series of Tests - #
#

class GdayTests(unittest.TestCase):
    
    def test_water_fluxes(self):
        (X, Y) = evap_trans()
        np.testing.assert_array_almost_equal(X, Y)
        (X, Y) = soil_water()
        np.testing.assert_array_almost_equal(X, Y, decimal=1)
        
    def test_carbon_fluxes(self):
        (X, Y) = NPP_resp()
        np.testing.assert_array_almost_equal(X, Y)
        (X, Y) = NPP_growth()
        np.testing.assert_array_almost_equal(X, Y)
        (X, Y) = NEP()
        np.testing.assert_array_almost_equal(X, Y)
    
    def test_carbon_stocks(self):  
        (X, Y) = cf_stock()
        np.testing.assert_array_almost_equal(X, Y, decimal=1)
    
         
if __name__ == "__main__":
    
    unittest.main()
