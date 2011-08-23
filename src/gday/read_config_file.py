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
    def __init__(self, fname, default_dir=None, fname_spec=None, 
                    default_spec_dir=None, ext='.cfg'):
        
        """
        Parameters:
        ----------
        fname : string
            filename of parameter (CFG) file
        default_dir : string
            path to the parameter file
        fname_spec : string
            filename of parameter (CFG) spec file. Spec file is identical to the
            CFG file except it list the various data types.
        default_spec_dir : string
            path to the parameter spec file
        ext : string
            file extension for parameter file, default = CFG
        
        """
        self.fname = fname
        self.default_dir = default_dir
        self.fname_spec = fname_spec
        self.default_spec_dir = default_spec_dir
        self.ext = ext
        
        if self.default_dir is None:
            self.default_dir = os.getcwd()
        self.config_file = os.path.join(self.default_dir, '%s%s' 
                                                % (self.fname, self.ext))
        
        # the configspec file should be in same directory as the config file?
        if self.default_spec_dir is None:
            self.default_spec_dir = self.default_dir
        if self.fname_spec is None:
            self.fname_spec = self.fname + 'spec'
            self.config_spec_file = os.path.join(self.default_spec_dir, '%s%s' 
                                                % (self.fname_spec, self.ext))
        else:
            self.config_spec_file = os.path.join(self.default_spec_dir, '%s%s' 
                                                % (self.fname_spec, self.ext))

    def load_files(self): 
        """ load config and configspec file, return a dictionary 
        
        Returns:
        --------
        config : object
            user defined parameter file as an object
        
        """
        try:
            config = ConfigObj(self.config_file, 
                                configspec=self.config_spec_file,
                                file_error=True)
            # check default values have been filled in
            validator = Validator()
            config.validate(validator)
            
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
    
    # pylint: disable=C0103
    default_dir = "/Users/mdekauwe/src/python/GDAY_model/params"
    fname = 'gday'
    T = LoadConfigFile(fname, default_dir=default_dir)
    config_dict = T.load_files()
    T.get_config_dicts(config_dict)
    
    
    