#!/usr/bin/env python
""" Load all the model initialisation data, see docstring below"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (22.02.2011)"
__email__   = "mdekauwe@gmail.com"


import re
import os
import sys
import keyword
import default_params as p
import default_control as c
import default_state as s
import default_files as fi
import default_fluxes
import ConfigParser
from utilities import str2boolean

def initialise_model_data(fname, met_header, DUMP=True):
    """ Load default model data, met forcing and return
    If there are user supplied input files initialise model with these instead

    Parameters:
    ----------
    fname : string
        filename of input options, parameters. Filename should include path!
    met_header : int
            row number of met file header with variable names
    DUMP : logical
        dump a the default parameters to a file

    Returns:
    --------
    control : integers, object
        model control flags
    params: floats, object
        model parameters
    state: floats, object
        model state
    fluxes : floats, object
        model fluxes
    met_data : floats, dictionary
        meteorological forcing data

    """
    # if this code is run in a monte carlo fashion python doesn't reimport the
    # modules with each instances! Force it to
    reload(default_fluxes)
    reload(s)
    reload(c)
    reload(fi)
    reload(p)
    
    R = ReadConfigFile(fname)
    config_dict = R.load_files()
    (user_control, user_params, user_state,
        user_files, user_fluxes, user_print) = R.get_config_dicts(config_dict)
    
    # get driving data
    forcing_data = read_met_forcing(user_files['met_fname'], met_header)

    # read in default modules and then adjust these
    if DUMP == False:
        params = adjust_object_attributes(user_params, p)
        state = adjust_object_attributes(user_state, s)
        control = adjust_object_attributes(user_control, c)
        files = adjust_object_attributes(user_files, fi)
    else:
        params = p
        state = s
        control = c
        files = fi
   
    return (control, params, state, files, default_fluxes, forcing_data,
            user_print)

class ReadConfigFile(object):
    """ Read supplied config file (.cfg/.ini).

    Return various dictionaries based on defined sections. 
    """
    def __init__(self, fname):

        """
        Parameters:
        ----------
        fname : string
            filename of parameter (CFG) file [including path]

        """
        self.config_file = fname
        self.Config = ConfigParser.ConfigParser()
        self.Config.optionxform = str # Respect case
    
    
    def load_files(self):
        """ load config file, return a dictionary

        Returns:
        --------
        config : object
            user defined parameter file as an object

        """
        try:
            config = self.Config.read(self.config_file)
        except (ConfigObjError, IOError), e: 
            raise IOError('%s' % e)
        
        return config
        
    def get_config_dicts(self, config_dict):
        """ reak config dictionary into small dictionaries based on sections.

        Parameters:
        -----------
        config_dict : dictionary
            User supplied config dict

        Returns:
        --------
        control : integers, object
            model control flags
        params: floats, object
            model parameters
        state: floats, object
            model state
        fluxes : floats, object
            model fluxes

        """
        user_files = self.ConfigSectionMap("files")
        user_params = self.ConfigSectionMap("params")
        user_control = self.ConfigSectionMap("control")
        user_state = self.ConfigSectionMap("state")
        user_print_opts = self.ConfigSectionMap("print")
        
        # add default cfg fname, dir incase user wants to dump the defaults
        fi.cfg_fname = self.config_file

        return (user_control, user_params, user_state, user_files,
                default_fluxes, user_print_opts)
    
    def ConfigSectionMap(self, section):
        flags = ['model_optroot', "deciduous_model", "modeljm", \
                 "water_stress", "fixleafnc", "passiveconst", "calc_sw_params",\
                 'fixed_stem_nc','exudation','adjust_rtslow', 'ncycle', 'grazing']
        flags_up = ["assim_model", "print_options", "alloc_model", "ps_pathway",\
                    "gs_model"]
        
        
        
        dict1 = {}
        options = self.Config.options(section)
        for option in options:
            try:
                value = self.Config.get(section, option)
                if section == "params" or section == "state":
                    
                    if value.replace('_','').replace('"','').isalpha() and value != "None":
                        dict1[option] = value
                    elif value.replace('_','').replace('"','').isalpha() and value == "None":
                        dict1[option] = None
                    else:
                        dict1[option] = float(value)
                elif section == "control":
                    if option in flags:
                        dict1[option] = str2boolean(value)
                    elif option in flags_up:
                        dict1[option] = value.replace('"', '').upper()
                    elif value.replace('_','').replace('"','').isalpha() and value != "None":
                        dict1[option] = value
                    elif value.replace('_','').replace('"','').isalpha() and value == "None":
                        dict1[option] = None
                    else:
                        dict1[option] = int(value)
                    
                elif section == "print":
                    dict1[option] = value.replace('"', '')
                elif section == "files":
                    dict1[option] = value.replace('"', '')
                else:
                    dict1[option] = self.Config.get(section, option)
                if dict1[option] == -1:
                    DebugPrint("skip: %s" % option)
            except:
                print("exception on %s!" % option)
                dict1[option] = None
                
        return dict1
    
    
#def read_met_forcing(fname, met_header, comment='#'):
#    """ Read the driving data into a dictionary, assumes user has provided
#    location of variables names.
#
#    Parameters:
#    -----------
#    fname : string
#        filename of parameter (CFG) file
#    met_header : int
#            row number of met file header with variable names
#    comment : string, optional
#        character defining a comment
#
#    Returns:
#    --------
#    data : dictionary
#        met forcing data
#
#    """
#    f = open(fname, 'rb')
#    # Skip crap
#    for i in xrange(met_header):
#        junk = csv.reader(f).next() 
#    names = ','.join(csv.reader(f).next()) 
#    col_names = re.sub(comment, ' ', names).lstrip().rstrip().split(",")
#    rows = csv.reader(f, delimiter=',') 
#    cols = map(list, zip(*rows)) # transpose the data
#    data = dict(zip(col_names, [map(float, c) for c in cols]))
#
#    return data
#

def read_met_forcing(fname, met_header, comment='#'):
    """ Read the driving data into a dictionary
    method searches for the hash tag followed by prjday in order to build the 
    named dictionary

    Parameters:
    -----------
    fname : string
        filename of parameter (CFG) file
    met_header : int
            row number of met file header with variable names
    comment : string, optional
        character defining a comment

    Returns:
    --------
    data : dictionary
        met forcing data

    """
    try:
        data = {}
        f = open(fname.replace('"', ''), 'r')
       
        for line_number, line in enumerate(f):            
            if line_number == met_header:
                # remove comment tag
                var_names = re.sub(r'#', ' ', line).lstrip().rstrip().split(",")
            elif not line.lstrip().startswith("#"):
                values = [float(i) for i in line.split(",")]
                for name, value in zip(var_names, values):
                    data.setdefault(name, []).append(value) 
        f.close()
    except IOError:
        raise IOError('Could not read met file: "%s"' % fname)

    return data
    
def adjust_object_attributes(user_dict, obj):
    """Loop through the user supplied dict and change relevant attributes

    Parameters:
    -----------
    user_dict : dictionary
        dictionary that contains values to change
    obj : object
        default model parameters

    Returns:
    --------
    obj : object
        adjusted parameters object

    """
    # check user hasn't specified a parameter we are not expecting...
    # make sure parameters is not named a reserved python word
    ignore = ['cfg_fname']
    special = ['topsoil_type', 'rootsoil_type']
    bad_words = keyword.kwlist
    bad_vars = [method for method in dir(str) if method[:2]=='__']
    for key, value in user_dict.iteritems():
        #print key, value
        if key in bad_words:
            err_msg = "You cant name your parameter anything from:\n\n %s" \
                            % bad_words
            raise RuntimeError, err_msg
        elif key in bad_vars:
            err_msg = "You cant name your parameter anything from:\n\n %s" \
                           % bad_words
            raise RuntimeError, err_msg
        elif key in special and not ignore:
            setattr(obj, key, value)
        elif key in ignore and not special:
            # we need to keep the default location names for the output dump
            setattr(obj, key, value)
        elif hasattr(obj, key):
            setattr(obj, key, value)
        else:
            err_msg = ".cfg file contains variable not in the model: %s" % key
            raise RuntimeError, err_msg
    return obj



if __name__ == "__main__":

    # pylint: disable=C0103

    fname = "../example/params/NCEAS_DUKE_model_youngforest_amb.cfg"
    met_header = 4
    
    # read in user defined variables (stored in dictionaries)
    (control, params, state, 
     files, fluxes, met_data, 
     print_opts) = initialise_model_data(fname, met_header, DUMP=False)
   
    print state.shootn
    
    par = met_data['par']
    sys.exit()
    
    import matplotlib.pyplot as plt
    plt.plot(par)
    plt.show()
