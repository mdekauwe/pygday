#!/usr/bin/env python

"""
Create a G'DAY met forcing file for the Duke site, used in NCEAS simulation.

Note we are outputting "daytime" averages not whole day averages...

And with the equilibrium simulation I am selecting random years from the 
available met-data

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (25.01.2011)"
__email__ = "mdekauwe@gmail.com"

import sys
import os
import csv
import math
import numpy as np
from datetime import date 
import datetime
import calendar
import scikits.timeseries as ts

def main(fname=None, ofname=None, start_sim=None, end_sim=None, equil=False, 
            num_sequences=None):
    
    
    ovar_names = ['prjday', 'sw_rad', 'tair', 'rain', 'tsoil', 
                    'tam', 'tpm', 'vpd_am', 'vpd_pm', 'vpd_avg', 'co2', 'ndep',
                    'wind', 'atmos_press']
    ounits = ['--', 'mj/m2/day', 'c', 'mm', 'c', 'c', 'c', 'kPA', 'kPa', 'kPa',
                'ppm', 't/ha/year', 'm/s', 'kPa']
    
    # rU -> universal newline support, as this was done in excel
    csv_data = csv.reader(open(fname, 'rU'), delimiter=' ')
    
    # Read the variable names and units
    var_names = csv_data.next()
    units = csv_data.next()
    
    # build empty dictionary
    data = {}
    for name in var_names:
        data[name] = []
    
    # read data into dictionary 
    for row in csv_data:
        values = [float(i) for i in row]
        for name, value in zip(var_names, values):
            data[name].append(value)
    
    # estimate DAYTIME tmean, tsoil.
    (tmean, tsoil, 
        tam, tpm) = estimate_temp_stuff(data['tmin'], data['tmax'])
    
    # convert PAR from mol to MJ
    sw_rad = convert_par_units(data['par'], rtn_par=False)
    
    co2 = [350.0] * len(data['year'])
    ndep = [0.016/365.25] * len(data['year'])
    
    # read wind and pressure..
    data_wind_pres = np.loadtxt("/Users/mdekauwe/research/NCEAS_face/met_data/day_avg_wind_sp_atmos_pressure.asc", skiprows=1)
    wind = data_wind_pres[:,2]
    pres = data_wind_pres[:,3]
    
    
    
    # calculate vpd stuff used in MATE, should we be checking consistency with
    # "observed" vpd?
    (vpd_am, vpd_pm, vpd_avg) = calculate_vpd_stuff(data['tmin'], data['tmax'])
    
    try:
        f = open(ofname, 'w')
        print >>f, "# Duke daily met forcing for NCEAS"
        print >>f, "# Data from %s-%s" % (start_sim, end_sim) 
        print >>f, "# Created by Martin De Kauwe, %s" % date.today()
        
        # print units out
        print >>f, '# Units:' ,
        for i in ounits:
            print >> f, i ,
        print >> f    
        
        # print variable names
        print >>f, '#' ,
        for i in ovar_names:
            print >> f, i ,
        print >> f 
        
    except IOError:
        raise IOError('Could not read met file: "%s"' % fname) 
        
    if equil == False:
        prj_day = 1
        for i in xrange(len(data['year'])):
            if data['year'][i] >= start_sim and data['year'][i] <= end_sim:
                
                print >> f, prj_day, sw_rad[i], tmean[i], data['rain'][i], \
                                        tsoil[i], tam[i], tpm[i], \
                                        vpd_am[i], vpd_pm[i], vpd_avg[i], \
                                        co2[i], ndep[i], wind[i], pres[i]
                prj_day += 1
    else:
        # For the Equilibirum weather we want random years repeated...
        prj_day = 1
        
        
        # how many repeats?
        num_seq = np.ones(num_sequences)
        yrs = np.arange(start_sim, end_sim + 1)
        num_years = len(yrs) * len(num_seq)
        
        # shuffle the years so we get "random" weather for our equilibrium
        # simulations
        shuff_years = (yrs * num_seq[:,None]).reshape(num_years)
        np.random.shuffle(shuff_years)
        year_doy = np.array([data['year'],data['doy']])
        for y in shuff_years:
            # get the index of a given shuffled year
            yrs_index = np.asarray(np.where(year_doy[0,:] == y))
            for i in yrs_index[0,:]:
                print >> f, prj_day, sw_rad[i], tmean[i], data['rain'][i], \
                                        tsoil[i], tam[i], tpm[i], \
                                        vpd_am[i], vpd_pm[i], vpd_avg[i]
                prj_day += 1
            
            
       
    f.close()

    


def estimate_temp_stuff(tmin, tmax):
    # G'DAY needs daytime mean temp/mean soil temp, not daily mean temp, 
    # so generate from tmin, tmax.
    # Routinue comes from original G'DAY code.
    #
    # Additionally temp routinues from MATE
    tmean, tsoil, temp_am, temp_pm = [], [], [], []
    for i in xrange(len(tmin)):
        
        # from gday
        avg = (tmax[i] + tmin[i]) / 2.0 + (tmax[i] - tmin[i]) / (3.0 * math.pi)
        soil = (tmax[i] + tmin[i]) / 2.0
        tmean.append(avg)   # original gday stuff
        tsoil.append(soil)  # original gday stuff 
        
        # followin stuff is from mate
        am = ((tmin[i] + tmax[i]) / 2.0 - (tmax[i] - tmin[i]) / 2.0 *
                math.sqrt(2.0) / 1.5 / math.pi)
        pm = ((tmin[i] + tmax[i]) / 2.0 + (tmax[i] - tmin[i]) / 2.0 * 
                (4.0 + 2.0 * math.sqrt(2.0)) / 3.0 / math.pi)
        temp_am.append(am) 
        temp_pm.append(pm)
        
    return tmean, tsoil, temp_am, temp_pm

def convert_par_units(par_in_moles, rtn_par=False):
    # return par in mj/m2/day, originally in mol/m2/day
    # photon of energy = h * c / lambda, h = planck, c = speed of light and
    # lambda = wavelength 
    # time conversion cancels out...
    
    #J_TO_MJ = 1.0E6
    #NANOM_TO_M = 1.0E-9 
    #planck_constant = 6.626068E-34 # j.s
    #speed_of_light = 299792458.0 # m/s
    #avogadro_no = 6.0221415E23
    #appx_par_wavel = 550 * NANOM_TO_M #(i.e. midpoint 400-700nm)
    # the proportion of PAR to global radiation, equal to 0.48 (McCree 1972)
    #PAR_TO_IRRADIANCE = 0.48 
    #MOLE_TO_MJ = (((avogadro_no * planck_constant * speed_of_light) / 
    #                appx_par_wavel) / J_TO_MJ)
    #sw_rad = []
    #for i in xrange(len(par_in_moles)):
    #    sw_rad.append((MOLE_TO_MJ * par_in_moles[i] / PAR_TO_IRRADIANCE))
    #    
    #return sw_rad
    
    # Note all of the above is identical to this...
    UMOLPERJ = 4.6
    
    # the proportion of PAR to global radiation, equal to 0.48 (McCree 1972)
    PAR_TO_IRRADIANCE = 1.0 / 0.48 
    
    rad = []
    for i in xrange(len(par_in_moles)):
        # return SW_radiation
        if rtn_par == False:
            rad.append(par_in_moles[i] / UMOLPERJ * PAR_TO_IRRADIANCE)
        # return par, just change units to MJ/m2/day
        else:
            rad.append(par_in_moles[i] / UMOLPERJ)
    
        #print rad[i] *  2.208, par_in_moles[i]
    return rad

def calculate_vpd_stuff(tmin, tmax):
    """ G'DAY/MATE use daytime daily vpd not whole day averages...
    
    output units are kPa
    """
    
    vpd_am, vpd_pm, vpd_avg = [], [], []
    for i in xrange(len(tmin)):
        tavg = tmin[i] + tmax[i] / 2.0
        vpd_am.append(0.61078 * (math.exp(17.269 * tavg / (237.3 + tavg)) - 
                        math.exp(17.269 * tmin[i] / (237.3 + tmin[i]))) * 
                        (1.0 - math.sqrt(2.0) / 1.5 / math.pi))

        vpd_pm.append(0.5 * 0.61078 * (math.exp(17.269 * tmax[i] / 
                        (237.3 + tmax[i])) - math.exp(17.269 * tmin[i] / 
                        (237.3 + tmin[i]))) * (1.0 + (4.0 + 2.0 * 
                        math.sqrt(2.0)) / 3.0 / math.pi))
        
        
        if vpd_am[i] < 0.05:
            vpd_am[i] = 0.05
        
        if vpd_pm[i] < 0.05:
            vpd_pm[i] = 0.05
        
        vpd_avg.append((vpd_am[i] + vpd_pm[i]) / 2.0)
        
    return vpd_am, vpd_pm, vpd_avg
    

    
if __name__ == "__main__":
    
    fname = "duke_daily_met.dat"
    ofname = "duke_equilibrium_metdata.gin"
    num_seq = 11 # 9 yrs of met data repeated 11 times, i.e. 99 yrs of met data
    main(fname=fname, ofname=ofname, start_sim=1993.0, end_sim=2007.0, 
            equil=True, num_sequences=num_seq)

    fname = "duke_daily_met.dat"
    ofname = "duke_metdata.gin"
    num_seq = 1 # 9 yrs of met data repeated 11 times, i.e. 99 yrs of met data
    main(fname=fname, ofname=ofname, start_sim=1998.0, end_sim=2007.0, 
            equil=False, num_sequences=num_seq)