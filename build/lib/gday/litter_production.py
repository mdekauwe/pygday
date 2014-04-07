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
        
        # ==================== 
        # C litter production
        # ====================
        self.fluxes.deadroots = rdecay * self.state.root   
        self.fluxes.deadcroots = self.params.crdecay * self.state.croot 
        self.fluxes.deadstems = self.params.wdecay * self.state.stem
        self.fluxes.deadbranch = self.params.bdecay * self.state.branch
        self.fluxes.deadsapwood = ((self.params.wdecay + 
                                    self.params.sapturnover) * 
                                    self.state.sapwood)
        
        if self.control.deciduous_model:
            self.fluxes.deadleaves = (self.fluxes.lrate * 
                                      self.state.remaining_days[doy])
        else:
            self.fluxes.deadleaves = fdecay * self.state.shoot
        
        
        # ==================== 
        # N litter production
        # ====================
        self.fluxes.deadleafn = self.fluxes.deadleaves * ncflit
    
        # Assuming fraction is retranslocated before senescence, i.e. a fracion 
        # of nutrients is stored within the plant
        self.fluxes.deadrootn = self.fluxes.deadroots * ncrlit
        self.fluxes.deadcrootn = (self.params.crdecay * self.state.crootn *
                                 (1.0 - self.params.cretrans))
        
        self.fluxes.deadbranchn = (self.params.bdecay * self.state.branchn *
                                  (1.0 - self.params.bretrans))
        
        # N in stemwood litter - only mobile n is retranslocated
        self.fluxes.deadstemn = (self.params.wdecay * 
                                (self.state.stemnimm + self.state.stemnmob *
                                (1.0 - self.params.wretrans)))
        
        
        
        # Animal grazing...
        #
        # Daily...
        if self.control.grazing == 1: 
            (self.fluxes.ceaten, 
             self.fluxes.neaten) = self.daily_grazing_calc(fdecay)
        # Once annually
        elif self.control.grazing == 2 and self.params.disturbance_doy == doy: 
            (self.fluxes.ceaten, 
             self.fluxes.neaten) = self.annual_grazing_calc()
        else: # no grazing
            self.fluxes.ceaten = 0.0
            self.fluxes.neaten = 0.0

        return (fdecay, rdecay)

    def daily_grazing_calc(self, fdecay):
        """ daily grass grazing...

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
    
    def annual_grazing_calc(self):
        """ annual grass grazing...single one off event
        
        Returns:
        --------
        ceaten : float
            C consumed by grazers [tonnes C/ha/day]
        neaten : float
            N consumed by grazers [tonnes C/ha/day]
        """
        ceaten = self.state.shoot * self.params.fracteaten
        neaten = self.state.shootn * self.params.fracteaten
        return (ceaten, neaten)
    
    def decay_in_dry_soils(self, decay_rate, decay_rate_dry):
        """Decay rates (e.g. leaf litterfall) can increase in dry soil, adjust
        decay param. This is based on field measurements by F. J. Hingston 
        (unpublished) cited in Corbeels.

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
        
        Reference:
        ----------
        Corbeels et al. (2005) Ecological Modelling, 187, 449-474.
        
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
