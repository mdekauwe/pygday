import os
import sys
import math
import constants as const

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

        # dump the default run parameters for the user to change
        self.default_param_fname = self.files.cfg_fname

        # dump the state at the end of a run, typical if user is running to
        # equilibrium
        self.out_param_fname = self.files.out_param_fname

        # daily output filename
        try:
            self.odaily = open(self.files.out_fname, 'w')
        except IOError:
            raise IOError("Can't open %s file for write" % self.odaily)

       
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
        ignore = ['actncslope', 'slowncslope', 'passncslope', 'decayrate', \
                    'fmfaeces', 'light_interception', 'wtfac_tsoil', \
                    'wtfac_root']
        self.dump_ini_data("[files]\n", self.files, ignore, oparams, 
                            print_tag=False, print_files=True)
        self.dump_ini_data("\n[params]\n", self.params, ignore, oparams, 
                            print_tag=False, print_files=False)
        self.dump_ini_data("\n[state]\n", self.state, ignore, oparams, 
                            print_tag=False, print_files=False)
        self.dump_ini_data("\n[control]\n", self.control, ignore, oparams, 
                            print_tag=False, print_files=False)
        self.dump_ini_data("\n[print]\n", self.print_opts, ignore, oparams, 
                            print_tag=True, print_files=False)
        
    def dump_ini_data(self, ini_section_tag, obj, ignore, fp, print_tag=False,
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
                fp.writelines("%s = %s\n" % (i, getattr(obj, i)) 
                                    for i in data)
            elif print_tag == False and print_files == True:
                 fp.writelines('%s = "%s"\n' % (i, getattr(obj, i)) 
                                for i in data)
            elif print_tag == True and print_files == False:
                fp.writelines('%s = "%s"\n' % (i, "yes") for i in obj)
        except IOError:
            raise IOError("Error writing params file")
    
    def write_daily_file_headers(self):
        """write (with comment, #) column headings"""
        try:
            self.odaily.write("%s " % '# prjday year doy')
            self.odaily.writelines("%s " % (i) for i in self.print_opts)
            self.odaily.write("\n")
        except IOError:
            raise IOError("Error writing file headers to: %s" % self.odaily)

    def save_daily_output(self, project_day, date):
        """ print daily output to file

        Parameters:
        -----------
        project_day : integer
            simulation day
        date : datetime object
            format yr/mth/day

        """
        # day of year 1-365/366
        doy = int(date.strftime('%j'))
        year = date.year
        if project_day == 1:
            self.write_daily_file_headers()
        try:
            self.odaily.write("%s %s %s " % (project_day, year, doy))
            for var in self.print_opts:
                try:
                    if hasattr(self.state, var):
                        value = getattr(self.state, var)
                        self.odaily.write("%s " % value)
                    else:
                        value = getattr(self.fluxes, var)
                        self.odaily.write("%s " % value)
                except AttributeError:
                    err_msg = "Error accessing var to print: %s" % var
                    raise AttributeError, err_msg
            self.odaily.write("\n")
        except IOError:
            raise IOError("Error writing file: %s" % self.odaily)
    
    def tidy_up(self):
        self.odaily.close()
