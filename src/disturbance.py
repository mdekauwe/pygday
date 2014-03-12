import random

class Disturbance(object):

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
        
    def time_till_next_disturbance(self, return_interval):
        """ calculate the number of years until a disturbance event occurs
        assuming a return interval of X years 
    
        - section 3.4.1 D. Knuth, The Art of Computer Programming.
        
        Parameters
        ----------
        return_interval : int/float
            interval disturbance return at in years
        """
        rate = 1.0 / return_interval
        
        return int(-math.log(1.0 - random.random()) / rate)

    def fire(self):
        pass
        
    
    def hurricane(self):
        pass