#!/usr/bin/env python
""" Load all the model initialisation data, see docstring below"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (22.02.2011)"
__email__   = "mdekauwe@gmail.com"

import os
import sys
import keyword
import default_params as p
import default_control as c
import default_state as s
import default_files as fi
import default_fluxes 
from read_config_file import LoadConfigFile

def initialise_model_data(fname, DUMP=True): 
    """ Load default model data, met forcing and return
    If there are user supplied input files initialise model with these instead
    
    Parameters:
    ----------
    fname : string
        filename of input options, parameters. Filename should include path!
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
    
    R = ReadUserConfigFile(fname) 
    config_dict = R.load_files()
    (user_control, user_params, user_state, 
        user_files, user_fluxes, user_print) = R.get_config_dicts(config_dict)
    
    # get driving data
    forcing_data = read_met_forcing(fname=user_files['met_fname'])
    
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


class ReadUserConfigFile(LoadConfigFile):
    """ Read user supplied config file.
    
    Return various dictionaries based on defined sections. Note this is a 
    sub-class of the generic LoadConfigFile class.
    """
    def get_config_dicts(self, config_dict):
        """ split the dictionary into chuncks
        
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
        user_files = config_dict['files']
        user_params = config_dict['params'] 
        user_control = config_dict['control'] 
        user_state = config_dict['state'] 
        user_print_opts = config_dict['print']
        
        # add default cfg fname, dir incase user wants to dump the defaults
        fi.cfg_fname = self.config_file
    
        return (user_control, user_params, user_state, user_files, 
                default_fluxes, user_print_opts)


def read_met_forcing(fname, comment='#'):        
    """ Read the driving data.
    Expects to find the data:
    prjday, sw_rad, tmean, rain, vpd, tsoil
    
    it searches for the hash tag followed by prjday in order to build the named
    dictionary
    
    Parameters:
    -----------
    fname : string
        filename of parameter (CFG) file
    comment : string, optional
        character defining a comment
        
    Returns:
    --------
    data : dictionary
        met forcing data
        
    """
    import sys
    try:
        f = open(fname, 'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            if 'prjday' in line:
                var_names = line[1:].split()
                data = {}
                for name in var_names:
                    data[name] = []
            elif not line.lstrip().startswith("#"):    
                values = [float(i) for i in line.split()]
                for name, value in zip(var_names, values):
                    data[name].append(value)
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
    bad_words = keyword.kwlist
    bad_vars = [method for method in dir(str) if method[:2]=='__']
    for key, value in user_dict.iteritems():
        if key in bad_words:
            err_msg = "You cant name your parameter anything from:\n\n %s" \
                            % bad_words
            raise RuntimeError, err_msg
        elif key in bad_vars:
            err_msg = "You cant name your parameter anything from:\n\n %s" \
                           % bad_words
            raise RuntimeError, err_msg
        elif hasattr(obj, key):
            setattr(obj, key, value)
        elif key in ignore:
            # we need to keep the default location names for the output dump
            setattr(obj, key, value)
        else:
            err_msg = ".cfg file contains variable not in the model: %s" % key
            raise RuntimeError, err_msg
    return obj


        
if __name__ == "__main__":
    
    # pylint: disable=C0103
    
    fname = 'gday'
    fdir = "/Users/mdekauwe/src/python/GDAY_model/params"
    
    # read in user defined variables (stored in dictionaries)
    (control, params, state, files, fluxes, 
        met_data) = initialise_model_data(fname, default_dir=fdir, DUMP=False) 
    
    print state.shootn 
    
    par = met_data['par']
    
    import matplotlib.pyplot as plt
    plt.plot(par)
    plt.show()