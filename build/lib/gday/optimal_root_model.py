""" Really need to fix these units so that there are consistent with G'DAY!!!"""

from math import exp
from utilities import float_ne

__author__  = "Martin De Kauwe"
__version__ = "1.0 (26.10.2012)"
__email__   = "mdekauwe@gmail.com"


class RootingDepthModel(object):
    """ Ross's Optimal rooting depth model.
    
    Optimisation hypothesis = maximise total N uptake by determining optimal 
    distribution of root mass with soil depth and optimal root depth
    
    References:
    ----------
    * McMurtire, R. et al (2011) Increased nitrogen-uptake efficiency at 
      elevated  CO2 explained by an hypothesis of optimal root function. Ecology 
      and Evolution, 2, 1235--1250
    """
    def __init__(self, zval, r0, top_soil_depth):
        """
        Parameters:
        -----------
        zval : float
            Length scale for exponential decline of Umax(z)
        r0 : float
            root C at half-maximum N uptake (kg C/m3)
        top_soil_depth : float
            depth of soil assumed by G'DAY, note Ross comment about 20 cm 
            [email]
            
        Returns:
        --------
        rtot : float
            Total root C mass
        
        """
        self.zval = zval    
        self.r0 = r0      
        self.top_soil_depth = top_soil_depth
        
    def main(self, rtoti=None, nsupply=None, depth_guess=None):
        """
        Parameters:
        -----------
        rtoti : float
            Initial fine root C mass -> from G'DAY [kg m-2] -> paper says DM?!
        depth_guess : float
            Initial guess at the rooting depth, used as the first point in the 
            root depth minimisation scheme [m]. 
        
        Returns:
        --------
        nuptake : float
            N uptake from roots [gN m-2 yr-1]
        """
        # step 2: determine maximum rooting depth for model
        root_depth = self.estimate_max_root_depth(rtoti, depth_guess)
        
        # step 6: calculate plant N uptake -> eqn B8.
        nuptake = self.calc_plant_nuptake(root_depth, nsupply)
        
        # step 7: daily root litter calculation. G'DAY requires root litter
        # input to the top 30 cm of soil, so the G'DAY needs changing. So
        # mortality below 30 cm should be ignored. Ross assumes that root
        # longevity and root N:C are independent of depth, thus avoiding the
        # integrals
        rabove = self.calculate_root_mass_above_depth(rtoti, root_depth)
        
        return (root_depth, nuptake, rabove)
    
    def estimate_max_root_depth(self, rtoti, depth_guess):
        """ Determing the maximum rooting depth through solving Eqn. B6. for 
        rooting depth
        
        Parameters:
        -----------
        dmax_iteration : float
            An iteration of the rooting depth defined by the minimisation scheme
            [NOT USED]
        rtoti : float
            Initial fine root root C mass [from G'DAY] 
            [kg m-2] -> paper says DM?!
        r0 : float
            Root C at half-max N uptake. [NOT USED]
        zval : float
            Length scale for exponential decline of Umax(z). [NOT USED]
        depth_guess : float
            initial starting guess at the root depth [m]
            
        Returns:
        --------
        rooting_depth : float
            Minimised rooting depth [m]
        
        """
        root_depth = newton(self.rtot_wrapper, self.rtot_derivative, 
                            depth_guess, args=(rtoti, self.r0, self.zval))
        
        # check minmised value is sensible?
        min_rtot = round(self.rtot(root_depth, rtoti, self.r0, self.zval), 4)
        gday_rtot = round(rtoti, 4)
        if float_ne(min_rtot, gday_rtot):
            msg = "Error, rtot supplied = %f but min = %f" % (rtoti, smin_rtot)
            raise RuntimeError(msg)
            
        return root_depth
        
    def rtot_wrapper(self, *args):    
        """ Wrapper method that calls rtot, but subtracts rtoti from the result
        to give you the rooting depth
        
        Parameters:
        -----------
        args[1] : float
            Initial fine root root C mass [from G'DAY], rtoti
        self.rtot : function
            call rtot function to estimate a value of rtot given a rooting
            depth iteration
        
        Returns:
        --------
        val  : float
            A minimised rooting depth iteration
        """
        return self.rtot(*args) - args[1]  
    
    def rtot(self, *args):
        """ Estimate the total root biomass per unit ground area, i.e. the 
        integral of root mass per unit soil volume, R(z), over depth z from the
        soil surface to the maximim rooting depth, dmax.
        
        Parameters:
        -----------
        dmax : float
            Rooting depth [m]
        rtoti : float
            Initial fine root root C mass [from G'DAY] 
            [kg m-2] -> paper says DM?!
        r0 : float
            Root C at half-max N uptake.
        zval : float
            Length scale for exponential decline of Umax(z)
        
        Returns:
        --------
        rtot : float
            Total root C mass given a rooting depth
        
        """
        (dmax, rtoti, r0, zval) = args
        return (r0 * (2.0 * zval * exp(0.5 * dmax / zval) - 
                (dmax + 2.0 * zval)))      
    
    def rtot_derivative(self, *args):
        """ Derivative of maximum root depth equation, rtot
        
        Parameters:
        -----------
        dmax : float
            Rooting depth [m]
        rtoti : float
            Initial fine root root C mass [from G'DAY]
            [kg m-2] -> paper says DM?!
        r0 : float
            Root C at half-max N uptake.
        zval : float
            Length scale for exponential decline of Umax(z)
        
        Returns:
        --------
        val : float
            derivative of rtot
        
        """
        (dmax, rtoti, r0, zval) = args
        return r0 * (1.0 * exp(0.5 * dmax / zval) - 1.0)
    

    def calculate_root_mass_above_depth(self, rtoti, root_depth):
        """ Estimate cumulative root mass above depth, 30 cm for the G'DAY model
        
        Parameters
        ----------
        rtoti : float
            Initial fine root root C mass [from G'DAY]
            [kg m-2] -> paper says DM?!
        root_depth : float
            model rooting depth (m)

        Returns
        -------
        val : float
            cumulative root C mass above soil depth assumed by G'DAY model, 30cm
        """
        arg1 = rtoti + 2.0 * self.r0 * self.zval + root_depth * self.r0
        arg2 = 1.0 - exp(-self.top_soil_depth / (2.0 * self.zval))
        
        return arg1 * arg2 - self.r0 * self.top_soil_depth
        
    def calc_plant_nuptake(self, *args):
        """ Plant N uptake as a func of maximum rooting depth, eqn B8
        
        Parameters
        ----------
        root_depth : float
            max rooting depth [m]
        z : float
            incremental depth provided by integration func
        nsupply : float
            soil N supply rate to plants per day [N/m2]
        top_soil_depth : float
            Depth of soil assumed by G'DAY model [m]
            
        Returns
        -------
        nuptake : float
            plant N uptake
        """
        (root_depth, nsupply) = args
        arg1 = nsupply / (1.0 - exp(-self.top_soil_depth / self.zval))
        arg2 = (1.0 - exp(-root_depth / (2.0 * self.zval)))**2
        
        return arg1 * arg2   


def newton(f, fprime, x0, args=(), tol=1E-6, maxiter=250):
    """ Newton-Raphson: finds a zero of the func, given an inital guess
    
    Parameters
    ----------
    f : function
        The function whose zero is wanted. 
    x0 : float
        An initial guess.
    fprime : function
        The derivative of the function 
    args : tuple, optional
        Extra arguments to be used in the function call.
    tol : float, optional
        The allowable error of the zero value.
    maxiter : int, optional
        Maximum number of iterations.

    Returns
    -------
    val : float
        Estimated location where function is zero
        
    """
    for iter in xrange(maxiter):
        myargs = (x0,) + args
        dx = f(*myargs) / fprime(*myargs)
        x = x0 - dx
        if abs(x - x0) < tol: 
            return x
        x0 = x
    raise RuntimeError, "No minimum found after %d iterations" % maxiter
