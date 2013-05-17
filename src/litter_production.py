""" Litter Production object """

__author__  = "Martin De Kauwe"
__version__ = "1.0 (25.02.2011)"
__email__   = "mdekauwe@gmail.com"

from utilities import float_lt, float_gt

class Litter(object):
    """ Calculate C and N litter production

    Litter production for each pool is assumed to be proportional to biomass
    pool size.
    """
    def __init__(self, control, params, state, fluxes):
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
        """

        self.params = params
        self.fluxes = fluxes
        self.control = control
        self.state = state
    
    def calculate_litter(self, doy=None):
        """Various C and N litter production elements

        Returns:
        --------
        fdecay : float
            foliage decay rate [tonnes C/ha/day]
        rdecay : float
            fine root decay rate [tonnes C/ha/day]

        """
        # Leaf/root litter rates are higher during dry periods and therefore is 
        # dependent on soil water content
        fdecay = self.decay_in_dry_soils(self.params.fdecay,
                                         self.params.fdecaydry)
        
        rdecay = self.decay_in_dry_soils(self.params.rdecay,
                                         self.params.rdecaydry)
        
        # litter N:C ratios, roots and shoot
        ncflit = self.state.shootnc * (1.0 - self.params.fretrans)
        ncrlit = self.state.rootnc * (1.0 - self.params.rretrans)
        
        # C litter production
        self.fluxes.deadroots = rdecay * self.state.root   # ditto
        self.fluxes.deadstems = self.params.wdecay * self.state.stem
        self.fluxes.deadbranch = self.params.bdecay * self.state.branch
        
        if self.control.deciduous_model:
            self.fluxes.deadleaves = (self.fluxes.lrate * 
                                      self.state.remaining_days[doy])
            self.fluxes.deadleafn = (self.fluxes.lnrate * 
                                     self.state.remaining_days[doy]  * 
                                     (1.0 - self.params.fretrans))
        else:
            self.fluxes.deadleaves = fdecay * self.state.shoot
            self.fluxes.deadleafn = self.fluxes.deadleaves * ncflit
        
        # N Litter production - assuming fraction is retranslocated before
        # senescence, i.e. a fracion of nutrients is stored within the plant
        self.fluxes.deadrootn = self.fluxes.deadroots * ncrlit
        self.fluxes.deadbranchn = (self.params.bdecay * self.state.branchn *
                                  (1.0 - self.params.bretrans))
        
        # N in stemwood litter - only mobile n is retranslocated
        self.fluxes.deadstemn = (self.params.wdecay * 
                                (self.state.stemnimm + self.state.stemnmob *
                                (1.0 - self.params.wretrans)))

        # any animals at work?
        if self.control.grazing:
            (self.fluxes.ceaten, self.fluxes.neaten) = self.grazing_calc(fdecay)
        else:
            self.fluxes.ceaten = 0.0
            self.fluxes.neaten = 0.0

        return (fdecay, rdecay)

    def grazing_calc(self, fdecay):
        """ grass grazing...

        Parameters:
        -----------
        fdecay : float
            foliage decay rate

        Returns:
        --------
        ceaten : float
            C consumed by grazers [tonnes C/ha/day]
        neaten : float
            N consumed by grazers [tonnes C/ha/day]
        """
        arg = (1.0 - self.params.fracteaten)
        ceaten = fdecay * self.params.fracteaten / arg * self.state.shoot
        neaten = fdecay * self.params.fracteaten / arg * self.state.shootn
        return (ceaten, neaten)

    def decay_in_dry_soils(self, decay_rate, decay_rate_dry):
        """Decay rates (e.g. leaf litterfall) can increase in dry soil, adjust
        decay param

        Parameters:
        -----------
        decay_rate : float
            default model parameter decay rate [tonnes C/ha/day]
        decay_rate_dry : float
            default model parameter dry deacy rate [tonnes C/ha/day]

        Returns:
        --------
        decay_rate : float
            adjusted deacy rate if the soil is dry [tonnes C/ha/day]

        """
        # turn into fraction...
        smc_root = self.state.pawater_root / self.params.wcapac_root
        
        new_decay_rate = (decay_rate_dry - (decay_rate_dry - decay_rate) * 
                         (smc_root - self.params.watdecaydry) /
                         (self.params.watdecaywet - self.params.watdecaydry))

        if float_lt(new_decay_rate, decay_rate):
            new_decay_rate = decay_rate

        if float_gt(new_decay_rate, decay_rate_dry):
            new_decay_rate = decay_rate_dry

        return new_decay_rate



if __name__ == "__main__":

    from file_parser import ConfigFileParser
    # pylint: disable=C0103
    # pylint: disable=C0324
    ini = "/Users/mdekauwe/src/python/GDAY_model/params/gday.ini"
    inispec = "/Users/mdekauwe/src/python/GDAY_model/params/gdayspec.ini"
    met = "/Users/mdekauwe/src/python/GDAY_model/forcing/ornl_met.dat"

    # read in user defined variables (stored in dictionaries)
    pars = ConfigFileParser(ini_fname=ini, inispec_fname=inispec, met_fname=met)
    (adj_control, adj_params, adj_state, adj_fluxes, forcing) = pars.main()

    lf = LitterFlows(adj_control, adj_params, adj_state, adj_fluxes)

    for d in xrange(len(forcing)):
        lf.flows()
        print adj_fluxes.deadleaves
