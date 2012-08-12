#!/usr/bin/env python
""" Load all the model initialisation data, see docstring below"""

__author__  = "Martin De Kauwe"
__version__ = "1.0 (22.02.2011)"
__email__   = "mdekauwe@gmail.com"


import sys
import os
from configobj import ConfigObj, ConfigObjError
from validate import Validator


class LoadConfigFile(object):
    """ Read a config file (.cfg/.ini) and return the data as a dictionary.

    This class should be sub-classed. Object expects that you are using the
    configspec file as well. I guess I should make this optional one day...
    """
    def __init__(self, fname):

        """
        Parameters:
        ----------
        fname : string
            filename of parameter (CFG) file [including path]

        """
        self.config_file = fname


    def load_files(self):
        """ load config file, return a dictionary

        Returns:
        --------
        config : object
            user defined parameter file as an object

        """
        try:
            config = ConfigObj(self.config_file, unrepr=True)
        except (ConfigObjError, IOError), e:
            raise IOError('%s' % e)
        return config

    def get_config_dicts(self, config_dict):
        """ Break config dictionary into small dictionaries based on sections.

        This is a dummy method and when subclassed should be adjusted such that
        it is code specific.
        """
        pass


class TestRead(LoadConfigFile):
    def process_config_dictionary(self, config_dict):
        print config_dict



if __name__ == "__main__":

    fname = "/Users/mdekauwe/src/python/pygday/params/duke_testing"
    T = LoadConfigFile(fname)
    config_dict = T.load_files()
    T.get_config_dicts(config_dict)
