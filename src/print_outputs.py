import os
import csv
import constants as const
from _version import __version__ as git_revision

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
        self.params = params
        self.state = state
        self.fluxes = fluxes
        self.control = control
        self.files = files
        self.print_opts = print_opts

        # get the git version of the code and stamp this into the output file
        self.revision_code = git_revision

        # dump the default run parameters for the user to change
        self.default_param_fname = self.files.cfg_fname

        # dump the state at the end of a run, typical if user is running to
        # equilibrium
        self.out_param_fname = self.files.out_param_fname
        
        if self.control.output_ascii:  
            self.odaily = open(self.files.out_fname, 'wb')
            self.wr = csv.writer(self.odaily, delimiter=',',
                                 quoting=csv.QUOTE_NONE, escapechar=' ')
        else:
            self.odaily = open(self.files.out_fname, 'wb')
            binhdr_fname = self.files.out_fname.split(".")[0] + '.bin.hdr' 
            self.odaily_hdr_fp = open(binhdr_fname, 'wb')
        # daily output file hdr
        try:
            self.print_fluxes = []
            self.print_state = []
            for i, var in enumerate(self.print_opts):
                #print var
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
                  'c_to_alloc_shoot','c_to_alloc_stem', 'lai',\
                  'n_to_alloc_branch', 'n_to_alloc_root', 'n_to_alloc_croot',\
                  'n_to_alloc_shoot', 'n_to_alloc_stem', 'n_to_alloc_stemimm',\
                  'n_to_alloc_stemmob', 'ncontent', 'fmleaf', 'fmroot',\
                  'branchnc', 'litterc', 'littercag' ,'littercbg',\
                  'littern', 'litternag', 'litternbg','plantc','plantn',\
                  'rootnc','shootnc','soilc', 'soiln','totalc',\
                  'totaln','wtfac_topsoil','wtfac_root','plantnc',\
                  'grw_seas_stress']

        special = ['rootsoil_type', 'topsoil_type', 'assim_model', 'co2_conc',\
                   'deciduous_model', 'fixleafnc', 'model_optroot',\
                   'passiveconst', 'print_options', 'water_stress',\
                   'calc_sw_params', 'alloc_model','fixed_stem_nc', \
                   'ps_pathway','gs_model','grazing','exudation',\
                   'ncycle','adjust_rtslow', "respiration_model"]
        
        self.dump_ini_data("[git]\n", None, ignore, special, 
                            oparams, print_tag=False, print_files=False, git=True)
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
                        print_files=False, git=False):
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
            
            if print_tag == False and print_files == False and git == False:
                for i in data:
                    fp.writelines("%s = %s\n" % (i, getattr(obj, i)))
            elif print_tag == False and print_files == False and git == True:
                fp.writelines('%s = %s\n' % ("git_hash", self.revision_code))
            elif print_tag == False and print_files == True and git == False:
                fp.writelines('%s = %s\n' % (i, getattr(obj, i)) for i in data)
            elif print_tag == True and print_files == False and git == False:
                fp.writelines('%s = %s\n' % (i, "yes") for i in obj)
            
        except IOError:
            raise IOError("Error writing params file")

    def write_daily_output_header(self):
        
        if self.control.output_ascii:
            self.wr.writerow(["%s:%s" % ("#Git_revision_code", self.revision_code.replace(" ", ""))])
            header = []
            header.extend(["year","doy"])
            header.extend(["%s" % (var) for var in self.print_state])
            header.extend(["%s" % (var) for var in self.print_fluxes])
            self.wr.writerow(header)        
        else:
            self.odaily_hdr_fp.write("%s:%s\n" % \
                                    ("#Git_revision_code", self.revision_code))
            buff1 = (",".join(("%s" % (var) for var in self.print_state)))
            buff2 = (",".join(("%s" % (var) for var in self.print_fluxes)))
            self.odaily_hdr_fp.write("year,doy,"+ buff1+","+ buff2+"\n")
            self.odaily_bin_hdr = "year,doy,"+ buff1+","+ buff2
            
            
    def write_daily_outputs_file(self, day_outputs):
        """ Write daily outputs to a csv file """
        self.wr.writerows(day_outputs)

    def write_daily_outputs_file_to_binary(self, day_outputs):
        """ Write daily outputs to a binay file 
        
        Not properly implemented:
            - Need to add something to make sure hdr info is attached to output
              file.
            - Also need a flag in the control to switch this on.
        """
        from numpy import array
        a = array(day_outputs,'float32')
        a.tofile(self.odaily)
        
        #from array import array
        #for i in xrange(len(day_outputs)):
        #    float_array = array('d', day_outputs[i])
        #    float_array.tofile(self.odaily)

    def clean_up(self, nrows=None):
        """ close the output file that holds the daily output """
        self.odaily.close()
        if not self.control.output_ascii:
            self.odaily_hdr_fp.write("nrows=%d\n" % (nrows))
            self.odaily_hdr_fp.write("ncols=%d\n" % \
                                     (len(self.odaily_bin_hdr.split(","))))
            self.odaily_hdr_fp.close()
