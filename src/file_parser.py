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
        user_files = self.buid_dict_from_ini_file("files")
        user_params = self.buid_dict_from_ini_file("params")
        user_control = self.buid_dict_from_ini_file("control")
        user_state = self.buid_dict_from_ini_file("state")
        user_print_opts = self.buid_dict_from_ini_file("print")

        # add default cfg fname, dir incase user wants to dump the defaults
        fi.cfg_fname = self.config_file

        return (user_control, user_params, user_state, user_files,
                default_fluxes, user_print_opts)

    def buid_dict_from_ini_file(self, section):
        """
        Return the .cfg file as a series of dictionaries depending on which section is called.

        the configparser package reads everything as a string, so we have to cast it
        ourselves. This isn't entriely straightforward, as we need to (i) catch None's,
        (ii) catch underscores as isalpha() ignores these and (iii) cast float/ints

        Parameters:
        -----------
        section : string
            Identifier to grab the relevant section from the .cfg file, e.g. "params"

        Returns:
        --------
        d : dictionary
            dictionary containing stuff from the .cfg file.
        """
        flags = ['model_optroot', "deciduous_model", "modeljm", \
                 'water_stress', "fixleafnc", "passiveconst", \
                 "calc_sw_params",'fixed_stem_nc','exudation',\
                 'adjust_rtslow', 'ncycle', 'grazing','output_ascii']
        flags_up = ["assim_model", "print_options", "alloc_model", \
                    "ps_pathway","gs_model", "respiration_model"]

        d = {}
        options = self.Config.options(section)
        for option in options:
            try:
                value = self.Config.get(section, option)
                if section == "params" or section == "state":
                    if value.replace('_','').isalpha() and value != "None":
                        d[option] = value
                    elif value.replace('_','').isalpha() and value == "None":
                        d[option] = None
                    else:
                        d[option] = float(value)
                elif section == "control":
                    if option in flags:
                        d[option] = str2boolean(value)
                    elif option in flags_up:
                        d[option] = value.upper()
                    elif value.replace('_','').isalpha() and value != "None":
                        d[option] = value
                    elif value.replace('_','').isalpha() and value == "None":
                        d[option] = None
                    else:
                        d[option] = int(value)
                elif section == "print" or section == "files":
                    d[option] = value
            except:
                print("Error reading .Cfg file into dicitonary: %s!" % option)
                d[option] = None

        return d


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
        f = open(fname, 'r')
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
