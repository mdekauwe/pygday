import os
import csv
import constants 			as const
import git_revision_info 	as git
import pdb
import numpy as np

__author__  = "Martin De Kauwe"
__version__ = "1.0 (21.03.2011)"
__email__   = "mdekauwe@gmail.com"



class PrintOutput(object):
    """Print model self.state and fluxes.

    Potential really to print anything, but for the moment the obvious.

    ** Note we are overiding the model data class here
    """
    def __init__(self, params, state, fluxes, control, files, print_opts):
        """
        Parameters
        ----------
        control : integers, object
            model control flags
        params: floats, object
            model parameters
        state: floats, object
            model state
        fluxes : floats, object
            model fluxes
        print_opts : object
            output print options

        """
        self.params 	= params
        self.state 		= state
        self.fluxes 	= fluxes
        self.control 	= control
        self.files 		= files
        self.print_opts = print_opts
        self.git 		= git
		
        # dump the default run parameters for the user to change
        self.default_param_fname = self.files.cfg_fname

        # dump the state at the end of a run, typical if user is running to
        # equilibrium
        self.out_param_fname = self.files.out_param_fname
        
        # make output files if the don't exist
        self.mk_output_dir()
			
        # daily output filename
        try:
            self.odaily = open(self.files.out_fname, 'wb')
            self.wr = csv.writer(self.odaily, delimiter=',', 
                                 quoting=csv.QUOTE_NONE, escapechar=' ')
            self.print_fluxes = []
            self.print_state = []
            for i, var in enumerate(self.print_opts):
                #print i, var
                try:
                    if hasattr(self.state, var):
                        self.print_state.append(var)
                    else:
                        self.print_fluxes.append(var)
                except AttributeError:
                    err_msg = "Error accessing var to print: %s" % var
                    raise AttributeError, err_msg
    
            self.write_daily_output_header()
        except IOError:
            raise IOError("Can't open %s file for write" % self.odaily)
        
        self.day_output = []
        
    def mk_output_dir(self):
    	if not os.path.isfile(self.files.out_fname):
        	dirname=''
        	dirs=self.files.out_fname.split('/')
        	
        	for i in dirs[:len(dirs)-1]:
        		dirname=dirname+i+'/'
        	
        	if not os.path.isdir(dirname): os.makedirs(dirname) 
     
    def get_vars_to_print(self):
        """ return lists of variable names to print out """
        return (self.print_state, self.print_fluxes)
      
    def save_default_parameters(self):
        """ Print default model state, control and param files.

        User should be adjusting these files as the defaults may be utter
        rubbish.
        """
        try:
            oparams = open(self.default_param_fname, 'w')
        except IOError:
            raise IOError("Can't open %s file for write" %
                            self.default_param_fname)
        self.print_parameters(oparams=oparams)

        # tidy up
        oparams.close()

    def save_state(self):
        """ Save model state

        Keep the state with the intention of starting the next model run from
        this point.

        Need to think about this as really, if we print to a new file i have to
        change the whole file read method. Might be easier to dump an entire
        param file as in print_parameters_default

        """
        try:
            oparams = open(self.out_param_fname, 'w')
        except IOError:
            raise IOError("Can't open %s file for write" %
                            self.out_param_fname)
        self.print_parameters(oparams=oparams)
        
        # tidy up
        oparams.close()

    def print_parameters(self, oparams=None):
        """ print model parameters

        This is either called before running the program so that the user can
        see the default parameters. Otherwise it is called at the end of
        run time, in which case it will represent the state of the model. This
        file can then be used to restart a run. A typical usage of this is
        running the model to some "equilibrium" state

        Parameters:
        -----------
        oparams : fp
            output parameter file pointer

        """
        ignore = ['anpp','actncslope', 'slowncslope', 'passncslope', 'decayrate', \
                  'fmfaeces', 'fipar', 'wtfac_tsoil', 'delta_sw_store',\
                  'remaining_days', 'growing_days', 'leaf_out_days', \
                  'albranch', 'alleaf', 'alroot', 'alstem', \
                  'c_to_alloc_branch', 'c_to_alloc_root', 'c_to_alloc_croot', \
                  'c_to_alloc_shoot','c_to_alloc_stem', \
                  'n_to_alloc_branch', 'n_to_alloc_root', 'n_to_alloc_croot',\
                  'n_to_alloc_shoot', 'n_to_alloc_stem', 'n_to_alloc_stemimm',\
                  'n_to_alloc_stemmob', 'ncontent', 'fmleaf', 'fmroot',\
                  'branchnc', 'lai', 'litterc', 'littercag' ,'littercbg',\
                  'littern', 'litternag', 'litternbg','plantc','plantn',\
                  'rootnc','shootnc','soilc', 'soiln','totalc',\
                  'totaln','wtfac_topsoil','wtfac_root','plantnc',\
                  'grw_seas_stress','git']

        special = ['rootsoil_type', 'topsoil_type', 'assim_model', 'co2_conc',\
                   'deciduous_model', 'fixleafnc', 'model_optroot',\
                   'modeljm', 'passiveconst', 'print_options', 'water_stress',\
                   'calc_sw_params', 'alloc_model','fixed_stem_nc', \
                   'ps_pathway','gs_model','grazing','exudation',\
                   'ncycle','adjust_rtslow']
                   
        self.dump_ini_data("[git]\n", self.git, ignore, special, oparams, 
                            print_tag=False, print_files=True)                    
        self.dump_ini_data("\n[files]\n", self.files, ignore, special, oparams, 
                            print_tag=False, print_files=True)
        self.dump_ini_data("\n[params]\n", self.params, ignore, special,oparams, 
                            print_tag=False, print_files=False)
        self.dump_ini_data("\n[state]\n", self.state, ignore, special, oparams, 
                            print_tag=False, print_files=False)
        self.dump_ini_data("\n[control]\n", self.control, ignore, special, 
                            oparams, print_tag=False, print_files=False)
        self.dump_ini_data("\n[print]\n", self.print_opts, ignore, special,
                            oparams, print_tag=True, print_files=False)
        
    def dump_ini_data(self, ini_section_tag, obj, ignore, special, fp, 
                        print_tag=False,
                        print_files=False):
        """ Get user class attributes and exclude builitin attributes
        Returns a list
    
        Parameters:
        ----------
        ini_section_tag : string
            section header 
        obj : object
            clas object
        ignore : list
            variables to ignore when printing output
        fp : fp
            out file pointer
        print_tag : logical
            if true special print options
        print_files : logical
            if true special print files options
        """
        try:
            fp.write(ini_section_tag)
            data = [i for i in dir(obj) if not i.startswith('__') \
                    and i not in ignore]
            data.sort()
            
            if print_tag == False and print_files == False:
                for i in data:
                    if i in special:
                        fp.writelines('%s = "%s"\n' % (i, getattr(obj, i)))
                    else:
                        fp.writelines("%s = %s\n" % (i, getattr(obj, i)))
                                    
            elif print_tag == False and print_files == True:
                fp.writelines('%s = "%s"\n' % (i, getattr(obj, i)) 
                                for i in data)
            elif print_tag == True and print_files == False:
                fp.writelines('%s = "%s"\n' % (i, "yes") for i in obj)
        except IOError:
            raise IOError("Error writing params file")
    
    def write_daily_output_header(self):
    	 self.write_daily_output_header_git()
    	 self.write_daily_output_header_titles()
    
    def write_daily_output_header_git(self):
    	
        header = ["URL","local branch","remote branch (if different)","revision code","tag"]
        info   = [self.git.URL_Fetch,self.git.branch,self.git.remote_branch,
        	      self.git.revision_code,self.git.tag]
        
        for i in range(len(header)): self.wr.writerow(["# " + header[i]+': '+info[i]])
        
    def write_daily_output_header_titles(self):
        header = []
        header.extend(["year","doy"])
        header.extend(["%s" % (var) for var in self.print_state])
        header.extend(["%s" % (var) for var in self.print_fluxes])
        
        self.wr.writerow(header)
        
    
    def write_daily_outputs_file(self, day_outputs):
        """ Write daily outputs to a csv file """
        self.wr.writerows(day_outputs)
        
   
    def clean_up(self):
        """ close the output file that holds the daily output """
        self.odaily.close()