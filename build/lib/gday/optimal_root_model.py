""" Really need to fix these units so that there are consistent with G'DAY!!!"""

from math import exp
from utilities import float_ne
import math

__author__  = "Martin De Kauwe"
__version__ = "1.0 (01.03.2013)"
__email__   = "mdekauwe@gmail.com"


class RootingDepthModel(object):
    """ Ross's Optimal rooting depth model.
    
    Optimisation hypothesis = maximise total N uptake (Utot) by determining 
    optimal distribution of root mass (R(z)) with soil depth and optimal
    root depth (D, m)
    
    References:
    ----------
    * McMurtire, R. et al (2011) Increased nitrogen-uptake efficiency at 
      elevated  CO2 explained by an hypothesis of optimal root function. Ecology 
      and Evolution, 2, 1235--1250
    """
    def __init__(self, d0x, r0, top_soil_depth):
        """
        Parameters:
        -----------
        d0x : float
            Length scale for exponential decline of Umax(z)
        r0 : float
            root C at half-maximum N uptake (kg C/m3)
        top_soil_depth : float
            depth of soil assumed by G'DAY, note Ross comment about 30 cm 
            [email]
            
        """
        self.d0 = d0x   
        self.r0 = r0      
        self.top_soil_depth = top_soil_depth
        
    def main(self, rtoti=None, nsupply=None, depth_guess=None):
        """
        Parameters:
        -----------
        rtoti : float
            Initial fine root C mass -> from G'DAY [kg m-2] -> paper says DM?!
        nsupply : float
            daily net N mineralisation in top soil layer from G'DAY
        depth_guess : float
            Initial guess at the rooting depth, used as the first point in the 
            root depth optimisation scheme [m]. 
        
        Returns:
        --------
        root_depth : float
            rooting depth [m]
        nuptake : float
            N uptake from roots [gN m-2 yr-1]
        rabove : float
            
        """
        # step 2: determine maximum rooting depth for model for a value rtoti
        # from G'DAY
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
            An iteration of the rooting depth defined by the optimisation scheme
            [NOT USED]
        rtoti : float
            Initial fine root root C mass [from G'DAY] 
            [kg m-2] -> paper says DM?! 
        depth_guess : float
            initial starting guess at the root depth [m]
            
        Returns:
        --------
        rooting_depth : float
            optimised rooting depth [m]
        
        """
        root_depth = newton(self.rtot_wrapper, self.rtot_derivative, 
                            depth_guess, args=(rtoti, self.r0, self.d0))
       
       
        
        #cons = ({'type': 'ineq', 'fun': self.constraint})
        #import scipy.optimize
        
        #args=(rtoti, self.r0, self.d0)
        #result = scipy.optimize.minimize(self.rtot_wrapper2, depth_guess, 
         #                                method="COBYLA",  constraints=cons,
        #                                 args=args)
        
        
       
        #print root_depth, result.x
        
       
        
       
       
        # check optimised value is sensible?
        min_rtot = round(self.rtot(root_depth, rtoti, self.r0, self.d0), 4)
        gday_rtot = round(rtoti, 4)
        if float_ne(min_rtot, gday_rtot):
            msg = "Error, rtot supplied = %f but min = %f" % (rtoti, smin_rtot)
            raise RuntimeError(msg)
            
        return root_depth
        
    def rtot_wrapper(self, *args):    
        """ Wrapper method that calls rtot, but subtracts rtoti from the result
        to give you the rooting depth (dmax)
        
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
            A optimised rooting depth iteration
        """
        return self.rtot(*args) - args[1]  
    
    def rtot_wrapper2(self, *args):
        sign = -1.0
        return (self.rtot(*args) - args[1]) * sign
    
    def rtot(self, *args):
        """ Estimate the total root biomass per unit ground area, i.e. the 
        integral of root mass per unit soil volume, R(z), over depth z from the
        soil surface to the maximim rooting depth, dmax. (Eqn 8, in McM 2012)
        
        Parameters:
        -----------
        dmax : float
            Rooting depth [m]
        rtoti : float
            Initial fine root root C mass [from G'DAY] 
            [kg m-2] -> paper says DM?!
        r0 : float
            Root C at half-max N uptake.
        d0 : float
            Length scale for exponential decline of Umax(z)
        
        Returns:
        --------
        rtot : float
            Total root C mass given a rooting depth
        
        """
        (dmax, rtoti, r0, d0) = args
        
        return r0 * (2.0 * d0 * (exp(dmax / (2.0 * d0)) - 1.0) - dmax) 
   
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
        d0 : float
            Length scale for exponential decline of Umax(z)
        
        Returns:
        --------
        val : float
            derivative of rtot
        
        """
        (dmax, rtoti, r0, d0) = args
        
        return r0 * (exp(0.5 * dmax / d0) - 1.0)
    
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
        arg1 = rtoti + 2.0 * self.r0 * self.d0 + root_depth * self.r0
        arg2 = 1.0 - exp(-self.top_soil_depth / (2.0 * self.d0))
        
        return arg1 * arg2 - self.r0 * self.top_soil_depth
        
    def calc_plant_nuptake(self, *args):
        """ Plant N uptake (Utot) as a func of maximum rooting depth
        
        This is the alternative eqn from McM word document
        
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
        arg1 = nsupply / (1.0 - exp(-self.top_soil_depth / self.d0))
        arg2 = (1.0 - exp(-root_depth / (2.0 * self.d0)))**2
        
        
        # Eqn B7
        #arg3 = nsup * d0
        #arg4 = (1.0 - exp(-root_depth/top_soil_depth))**2
        #print arg1 * arg2, arg3*arg4 
        
        return arg1 * arg2   
    
    def constraint(self, x):
        """ root_depth > 0.0 
        
        returns a positive number if within bound and 0.0 it is exactly on the 
        edge of the bound
        """
        return x
   

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
    
    
if __name__ == "__main__":
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    # Create the PdfPages object to which we will save the pages:
    pdf = PdfPages('/Users/mdekauwe/Desktop/root_model_test.pdf')
    
    #========Params=======
    d0x = 0.35
    r0 = .265#0.1325
    top_soil_depth = 0.3
    depth_guess=1.0
    #=====================
    rtot = np.linspace(0, 1.4, 100)
    nsupply = np.linspace(0.1, 10.0, 5)
    zval = np.linspace(0.05, .99, 5)
    
    RM = RootingDepthModel(d0x=d0x, r0=r0, top_soil_depth=top_soil_depth)
    for nsup in nsupply:
        dmax = np.zeros(0)
        for rt in rtot:
            (root_depth, nuptake, rabove) = RM.main(rt, nsup, 
                                                    depth_guess=depth_guess)
            dmax = np.append(dmax, root_depth)
        
        plt.plot(rtot, dmax, label="N supply = %.2f" %(nsup))
    plt.xlim(0, 1.0)
    plt.legend(numpoints=1, loc="best")
    plt.ylabel("Maximum rooting depth (D$_{max}$, m)")
    plt.xlabel("Total root mass (R$_{tot}$, kg DM m$^{-2}$)")
    plt.title("Fig2a: d$_0$ = %.2f, r$_0$ = %.4f, top_soil_depth = %.2f" % \
              (d0x, r0, top_soil_depth))
    plt.rcParams.update({'legend.fontsize': 8})
    pdf.savefig() 
    plt.clf()
    
    for nsup in nsupply:
        nup = np.zeros(0)
        for rt in rtot:
            (root_depth, nuptake, rabove) = RM.main(rt, nsup, 
                                                    depth_guess=depth_guess)
            nup = np.append(nup, nuptake)
        
        plt.plot(rtot, nup, label="N supply = %.2f" %(nsup))
    plt.xlim(0, 1.4)
    plt.legend(numpoints=1, loc="best")
    plt.ylabel("N Uptake (g N m$^{-2}$)")
    plt.xlabel("Total root mass (R$_{tot}$, kg DM m$^{-2}$)")
    plt.title("d$_0$ = %.2f, r$_0$ = %.4f, top_soil_depth = %.2f" % \
              (d0x, r0, top_soil_depth))
    plt.rcParams.update({'legend.fontsize': 8})
    pdf.savefig() 
    plt.clf()
    
    
    for zv in zval:
        RM = RootingDepthModel(d0x=zv, r0=r0, top_soil_depth=top_soil_depth)
        nup = np.zeros(0)
        for rt in rtot:
            (root_depth, nuptake, rabove) = RM.main(rt, nsup, 
                                                    depth_guess=depth_guess)
            nup = np.append(nup, nuptake)
        
        plt.plot(rtot, nup, label="soil depth (Z) = %.2f" %(zv))
    plt.xlim(0, 1.4)
    #plt.ylim(0, 1.0)
    plt.legend(numpoints=1, loc="best")
    plt.ylabel("N Uptake (g N m$^{-2}$)")
    plt.xlabel("Total root mass (R$_{tot}$, kg DM m$^{-2}$)")
    plt.title("d$_0$ = %.2f, r$_0$ = %.4f, top_soil_depth = %.2f" % \
              (d0x, r0, top_soil_depth))
    plt.rcParams.update({'legend.fontsize': 8})
    pdf.savefig() 
    plt.clf()
    #"""
    
    
    
    RM = RootingDepthModel(d0x=d0x, r0=r0, top_soil_depth=top_soil_depth)
    
    for nsup in nsupply:
        n_uptake_frac = np.zeros(0)
        for rt in rtot:
            (root_depth, nuptake, rabove) = RM.main(rt, nsup, 
                                                    depth_guess=depth_guess)
            
            # Total potential annual N uptake integrated over all soil depths
            Umax = nsup / (1.0 - exp(-top_soil_depth / d0x)) 
            
            # annual total N uptake per unit land area
            Utot = nuptake
           
            phi_N = Utot / Umax
            n_uptake_frac = np.append(n_uptake_frac, phi_N)
            
        plt.plot(rtot, n_uptake_frac, label="N supply = %.2f" %(nsup))
    
    plt.legend(numpoints=1, loc="best")
    plt.xlim(0, 1.4)
    plt.ylim(0, 1.0)
    plt.ylabel("Gross N-uptake fraction")
    plt.xlabel("Total root mass (R$_{tot}$, kg DM m$^{-2}$)")
    plt.title("Fig2b: d$_0$ = %.2f, r$_0$ = %.4f, top_soil_depth = %.2f" % \
              (d0x, r0, top_soil_depth))
    plt.rcParams.update({'legend.fontsize': 8})
    pdf.savefig() 
    plt.clf()
    
    # Remember to close the object
    pdf.close()
    
   