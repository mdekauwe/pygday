import random
import sys
from math import log

class Disturbance(object):

    def __init__(self, control, params, state, fluxes, met_data):
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
        self.met_data = met_data
                
    def initialise(self, years):
        
        
        yrs_till_event = self.time_till_next_disturbance()
        year_of_disturbance = years[0]+yrs_till_event
        
        # figure out the years of the disturbance events 
        self.yrs = []
        while year_of_disturbance < years[-1]:
            
            index = self.index_id(years, year_of_disturbance)
            
            yrs_till_event = self.time_till_next_disturbance()
            self.yrs.append(year_of_disturbance)
            
            # See if there is another event?
            year_of_disturbance = years[index]+yrs_till_event
        
    
    def index_id(self, a_list, elem):
        return (index for index, item in enumerate(a_list) 
                if item == elem).next()
    
    def check_year(self, year):
        if year in self.yrs:
            return True
        else:
            return False
    
    def time_till_next_disturbance(self):
        """ calculate the number of years until a disturbance event occurs
        assuming a return interval of X years 
    
        - section 3.4.1 D. Knuth, The Art of Computer Programming.
        
        Parameters
        ----------
        return_interval : int/float
            interval disturbance return at in years
        """
        rate = 1.0 / self.params.return_interval
        
        return int(-log(1.0 - random.random()) / rate)

    def disturbance(self):
        """
        Fire...
        
        * 100 percent of aboveground biomass 
        * 100 percent of surface litter
        * 50 percent of N volatilized to the atmosphere
        * 50 percent of N returned to inorgn pool"
        
        vaguely following ...
        http://treephys.oxfordjournals.org/content/24/7/765.full.pdf
        """
        totaln = (self.state.branchn + self.state.shootn + self.state.stemn + 
                  self.state.structsurfn)
        self.state.inorgn += totaln / 2.0
        
        # re-establish everything with C/N ~ 25.
        if self.control.alloc_model == "GRASSES":
            self.state.branch = 0.0
            self.state.branchn = 0.0
            self.state.sapwood = 0.0
            self.state.stem = 0.0
            self.state.stemn = 0.0
            self.state.stemnimm = 0.0
            self.state.stemnmob = 0.0
        else:
            self.state.branch = 0.001
            self.state.branchn = 0.00004
            self.state.sapwood = 0.001
            self.state.stem = 0.001
            self.state.stemn = 0.00004
            self.state.stemnimm = 0.00004
            self.state.stemnmob = 0.0
        
        self.state.age = 0.0
        self.state.canht = 0.1
        self.state.metabsoil = 0.0
        self.state.metabsoiln = 0.0
        self.state.metabsurf = 0.0
        self.state.metabsurfn = 0.0
        self.state.prev_sma = 1.0
        self.state.root = 0.001
        self.state.rootn = 0.00004
        self.state.shoot = 0.001
        self.state.shootn = 0.00004
        self.state.structsurf = 0.001
        self.state.structsurfn = 0.00004  
        
        