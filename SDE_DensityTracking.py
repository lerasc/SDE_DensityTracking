
import numpy            as np
import pandas           as pd

from warnings           import warn 
from scipy.integrate    import quad
from numbers            import Number
from inspect            import signature
from numpy 	            import sqrt, exp, pi
from scipy.interpolate  import interp1d, interp2d


class SDE_DensityTracking:
    """
    Consider a stochastic differential equation of the form
        dX_t 	= mu(X_t, t) dt + sigma(X_t, t) dW_t,    X_0 = x0 													(1)
    where mu and sigma are two functions of X and time t, and x0 is some initial value. In this class, we calculate
    numerically the probability p(x,t) of being at position x at time t. Instead of solving the Fokker-Planck equation
    associated with equation (1), we use a method called density tracking by quadrature (DTQ) [1], which first
    disretizes (1) in time as
        X_{i+1} = X_i + mu(X_t,t) h + sigma(X_i,t) sqrt(h) Z_{i+1},												    (2)
    with h the discretization time-step, and Z_{i+1} a Gaussian random variable with zero mean and unit variance.
    Then,  we apply the discretized Chapman-Kolmogorov equation to get
        p(x, t_{i+1}) = k sum_{j=-M}^M G(x, y_j) p(y_j, t_i) 														(3)
    where k is the spatial discretization width, G(x,y) is the density of being at position x at time i+1 when at
    position y at time t, and M is some finite size cut-off.
    See description.pdf and [1] for details. A similar implementation exists for R [2].
    Submit questions at sandrolera@gmail.com

    input:
    -----
    mu:			Drift function in equation (1), must be function of one variable (X) or two variables (X and t).
                First argument (location X) must always be present. Second argument (time t) is optional.

    sigma:		Same as mu, but the volatility function in equation (1).

    x0:         Special initial condition where p0 is a Dirac delta located at x=x0. Either p0 or x0 must be provided, 
                but not both.     

    p0:			Initial distribution, i.e. distribution at time t=0. The user must provide a triple that consist of an
                initial pdf (a functional object of one variable), as well as a lower/upper value, meaning that values
                below/above this value can be considered roughly equal to zero. For instance, if we have an initial 
                Gaussian distribution centred at 3 with standard deviation 1, we may assume that 5 standard deviations 
                away from 3, the distribution is rougly equal to zero. We then write: 
                    p0 = ( lambda x: np.exp( -(x-3)**2)/2 ) / np.sqrt( 2 * np.pi ), -2, 8 ) 

    x_:			Location of boundary, if any. Implementation is such that x0 > x_ is required, i.e. we assume the
                boundary is below the starting point. Currently, only one boundary is supported. See ref. [3] for a 
                generalization to two boundaries. 

    boundary:	Either 'reflective'	or 'absorbing', depending on whether the boundary at x_ is reflective or absorbing.
                Ignored if x_ is None. Note that the case boundary='reflective' is not yet implemented. But this is easy
                to do, cf. for instance equation (29) in [4]. 

    T:			Final time considered. Execution is stopped once t > T, i.e. after a total of T/h iterations. 

    h:			Temporal discretization step.

    k:			Spatial discretization step.

    tol:		Precision tolerance. At time t, we consider values between x_min(t) and x_max(t), which are determined
                as largest x value with  p(x_min,t) < tol and smallest x value with p(x_max,t) < tol. p(x) values with
                x < x_min or x > x_max are set to zero. 
                Note that this may cause problems if the initial or reflected density is bimodal with density < tol 
                between the two peaks. Because this is a rather rare case, we allow for this imperfection for 
                simplicity. Improving on this is straight forward. 
                Note that tol generally depends on the value range over which the x-values extend. Since the integral 
                of p between x_min and x_max is (roughly) equal to 1, we expect tol to be a value > 1 if 
                x_max-x_min << 1 and we expect tol << 1 if x_max - x_min >> 1. 


    references:
    ----------
    [1]	Bhat, H.S., and Madushani, R.W.M.A, "Density tracking by quadrature for stochastic differential equations."
        arXiv preprint arXiv:1610.09572 (2016).
    [2]	Rdtq: Density Tracking by Quadrature (https://cran.r-project.org/web/packages/Rdtq/index.html)
    [3] Veestraeten, D. "The conditional probability density function for a reflected Brownian motion."
        Computational Economics 24.2 (2004): 185-207.
    [4] Molini, A., et al. "First passage time statistics of Brownian motion with purely time dependent drift and
        diffusion." Physica A: Statistical Mechanics and its Applications 390.11 (2011): 1841-1852.        
    """

    ####################################################################################################################
    ### constructor
    ####################################################################################################################
    def __init__(self,
                 mu,
                 sigma,
                 p0			= None,
                 x0 		= None,
                 x_ 		= None,
                 boundary 	= None,
                 T 			= 10,
                 h 			= 0.05,
                 k 			= 0.05,
                 tol 		= 1e-4,
                 ):

        # check that provided input is correct
        ################################################################################################################
        assert( callable(mu) ),"mu must be function"
        no_mu_args  =  len(signature(mu).parameters) # number of arguments the mu function takes
        assert( no_mu_args in [1,2] ),"mu must be function of one or two variables"

        assert( callable(sigma) ),"sigma must be function"
        no_sig_args =  len(signature(sigma).parameters) # number of arguments the sigma function takes
        assert( no_sig_args in [1,2] ),"sigma must be function of one or two variables"

        if p0 is not None:
            assert(len(p0)==3),"p0 must be a tripple of values"
            p, lower, upper = p0 
            assert( callable(p) ),"the first value in the p0 triple must be a function"
            assert( len(signature(p).parameters)==1 ),"the first value in the p0 triple takes only one variable"
            if abs( quad(p,-np.inf,np.inf)[0] - 1 ) < 1e-3: warn("p0 not normalized to 1")
            assert(lower < upper),"lower boundary (2nd value in 3-tuple) must be larger than upper boundary (3rd value)"

        if x_ is not None:
            assert( isinstance(x_,Number) ),"x_ must be a number"

        if x0 is not None:
            assert( isinstance(x0,Number) ),"x0 must be a number"
            if x_ is not None: assert( x0 > x_ ),"x0 > x_ required"

        if p0 is None and x0 is None:
            raise ValueError("Either p0 or x0 must be specified.")

        if p0 is not None and x0 is not None:
            raise ValueError("Only p0 or x0 must be specified.")

        if boundary is not None:
            assert(boundary in ['absorbing','reflective']),"boundary must be 'reflective' or 'absorbing'."
            if x_ is None: boundary = None # for internal consistency

        assert(isinstance(T,Number)),"T must be a number"
        assert(isinstance(h,Number)),"h must be a number"    
        assert(isinstance(k,Number)),"k must be a number"    
        assert(isinstance(tol,Number)),"tol must be a number"

        # if mu/sigma was provided with only on variable, add additional t-dependence argument for internal consistency
        ################################################################################################################
        if no_mu_args==1:  mu_func = lambda x,t: mu(x)
        else:              mu_func = lambda x,t: mu(x,t)

        if no_sig_args==1: sig_func = lambda x,t: sigma(x)
        else:              sig_func = lambda x,t: sigma(x,t)


        # define class atributes
        ################################################################################################################
        self._mu 		    = mu_func
        self._sigma 	    = sig_func
        self._p0		    = p0
        self._x0 		    = x0
        self._x_ 		    = x_
        self._boundary	    = boundary
        self._T 		    = T
        self._h			    = h
        self._k 		    = k      
        self._tol 		    = tol
        self._ts  		    = np.arange(0,T+h,h) 	        # discrete times considered, incl. 0 and T.
        self._ns    	    = len(self._ts)			        # number of time slices
        self._tslices	    = np.arange(1,self._ns)	        # iteration of time slices
        self._calculated    = False					        # set two True once 'track_density' has been called


    ####################################################################################################################
    ### private methods
    ####################################################################################################################
    def _G(self,x,y,t):
        """
        Transition density of being at position x at time t+1 when at position y at time t.
        See e.g. [4] for a derivation of the density with absorbing and reflective barrier.
        """

        if self._boundary is None: 							# free transition density, equation (4) in [1]

            sqt 	= 2 * self._sigma(y,t)**2 * self._h
            p 		= exp(- (x-y-self._mu(y,t)*self._h)**2 / sqt )
            p 	   /= sqrt( pi * sqt )
            return p


        elif self._boundary=='absorbing':					# vanishing probability at x=x_, equation (26) in [3]

            if x <= self._x_: return 0.0                    # by definition, zero mass below x_

            y 	   -= self._x_								# rescale such that absorbing boundary is at y=0
            x 	   -= self._x_								# rescale such that absorbing boundary is at y=0
            sqt 	= 2 * self._sigma(y,t)**2 * self._h
            dr 		= self._mu(y,t) * self._h
            p1  	= exp( -( y+dr-x )**2 / sqt )
            p2	    = exp( - 2*y*self._mu(y) / self._sigma(y,t)**2 )
            p2     *= exp( -( y-dr+x )**2 / sqt )
            p       = (p1-p2) / sqrt( pi *sqt )
            return p

        else:												# vanishing derivative at x=x_, equation (29) in [3]

            raise NotImplementedError("Not yet implemented, but easy to do so.")


    def _CK(self, x, ys, t, ps, efficient=False):
        """
        Chapman-Kolmogorov step.
        Given a value x at which we want to calculate the density for the next time-step, integrate over all ys at the
        current time-step to obtain that value, using equation (7) in [1], with a natural truncation resulting from
        the finiteness of ys. The array ps stores the probabilities of the ys. If efficient is set to True, the
        convolution is not extended over all y values, but only the y-values that are within a range of 4 standard
        deviations from x, since G(x,y,t) is a Gaussian and hence strongly localized.
        """
        if not efficient:

            return self._k * sum( np.array([self._G(x,y,t) for y in ys]) * ps )       # integrate over all y values

        else: 
            std             = self._sigma(x,t) * np.sqrt(self._h)                     # one standard deviation
            intv            = pd.Interval( x-4*std, x+4*std )                         # outside this interval: approx. 0
            summands        = [ self._G(x,y,t)*p for y,p in zip(ys,ps) if y in intv ] # restricted integration range
            return            self._k * sum(summands)                              


    def _initialize(self):
        """
        Initialization of the t=0 probability array. Depending on whether this is an initial point or an intial
        distribution, we must initialize differently. 
        """
        if self._x0 is not None: 

            x_vals     = [ self._nx(self._x0) ]			                   # we have only one location
            p_vals     = [ 1/self._k ]						               # and one associated probability

        else:												       

            p0, x_min, x_max    = self._p0                                 # unpack the input tripple 
            x_min               = self._nx(x_min)                          # round to next multiple of k for consistency
            x_max               = self._nx(x_max)                          # round to next multiple of k for consistency
            x_vals              = np.arange(x_min, x_max+self._k, self._k) # from x_min to x_max in units of k
            x_vals              = [ self._nx(x)  for x in x_vals ]         # avoid rounding errors
            p_vals              = [ p0(x)        for x in x_vals ]         # disretized p-values at x values

            if len(x_vals)==0: raise ValueError("Couldn't find reasonable x_min and x_max for initial distribution.")

        self._xs   += [ x_vals ] 							        # x-values for first time slice
        self._ps   += [ p_vals ] 							        # probailities for first time slice
    

    def _reshape_results(self):
        """
        When running 'track_density', we collect lists of p(x,i) values, one such list per time step i. The different
        p(x,i) arrays or of different length, and generally increase in length (as the probability diffuses). In this
        method, we reshape this list of lists into a DataFrame and interpolate the data for a functional representation.
        """

        # create DataFrame of discrete values
        ################################################################################################################
        pxs   = []                                          # collects lists of p(x,i) arrays as Series
        for (x,p) in zip(self._xs, self._ps):               # for each time slice, transform the list of x an p values..
            pxs += [ pd.Series(p,index=x) ]                 # .. into a TimeSeries with x values as axis
        probs = pd.concat(pxs,axis=1).T                     # merge into one Frame, axis is time, columns are x-values
        probs = probs.replace(np.nan, 0.)                   # because NaN means prob < tol, i.e. approximately 0
        probs.index = self._ts                              # index is time axis
        probs.index.name = 'time'                           # label

        # create interpolated function p(x,t)
        ################################################################################################################
        self._probs_discrete    = probs
        try:  self._probs_functional  = interp2d(   x            =  probs.columns,
                                                    y            =  probs.index,
                                                    z            =  probs.values,
                                                    kind         = 'linear',
                                                    bounds_error =  False,
                                                    fill_value   =  0
                                                    )
        except: raise ValueError("Error in interp2d. Probably not enough datapoints. Decrease h and/or k.")


    def _interpolate_Dirac(self):
        """
        Special case: If initial condition is a Dirac Delta (x0 is not None), and the return value of 'get_density' is
        a functional expression at t=0, then we approximate the initial condition as a Gaussian centred around at x0
        with a standard devation equal to k/10.
        """

        sig  = 0.1* self._k
        def sharp_Gaussian(x): return np.exp(-0.5 * ((x-self._x0)/sig)**2) / (np.sqrt(2*np.pi) * 0.5 )

        return sharp_Gaussian


    def _i(self,x):
        """
        Given a location x, retrieve the associated index i, defined through x/k = i or the next rounded integer. 
        """
        return int(np.round( x/self._k ))


    def _nx(self, x):
        """
        Given location x, return the next closest value that is a multiple of k. Also used to avoid rounding errors.
        """
        return self._k * self._i(x)


    @staticmethod
    def _find_nearest(array, value):
        """
        Find and return value inside array that is closest to the value provided as input.
        """
        idx = (np.abs(array - value)).argmin()  # closest index
        return array[idx]


    ####################################################################################################################
    ### public methods
    ####################################################################################################################
    def track_density(self, print_status=False):
        """
        Core method that iterates over time and calculates the probability distribution for each time-slice.
        If print_status is True, the iteration progress is printed to the output
        """

        # initalize probability at t=0 and interate up to time i=T
        ################################################################################################################
        self._xs 		= []	# stores the considered x values between x_min(t) and x_max(t).
        self._ps 		= []	# stores lists of p(x)=p(x,t) for different t values.
        self._initialize()		# calculate i=0

        for i, t in enumerate(self._tslices): # for i=1,...,T do:

            # print the progress
            ############################################################################################################            
            percentage = 100 * (i+1) / len(self._tslices)
            if print_status: print("Simulation running, %.1f%% done."%percentage, end="\r")

            # extract x-values from previous time slice
            ############################################################################################################
            prev_x_vals     = self._xs[-1]                      # previous slice of x values
            prev_p_vals     = self._ps[-1]                      # previous slice of x values

            # if previous slice is empty, it means that in the previous step all probability densities were < tol
            ############################################################################################################
            if len(prev_x_vals)==0:                 # no probability mass (can happen if boundary='absorbing')
                self._xs   += [ [] ]                # dummy empty list that will propagate until t=T
                self._ps   += [ [] ]                # dummy empty list that will propagate until t=T
                continue                 

            # extract previously largest value. This is where we start the iteration for this time-slice
            ############################################################################################################
            prev_max        = np.array(prev_p_vals).argmax()    # index of previous maximum of bimodality
            x_start         = prev_x_vals[prev_max]             # where we start this iteration from            
                       # jump to next time-step

            # we start from x_start, and move downwards in steps of k until x_min(i) is reached
            ############################################################################################################
            x_vals_below 	= []
            p_vals_below    = []
            new_x 			= x_start
            new_p 			= self._CK(new_x, prev_x_vals, t, prev_p_vals)

            while new_p > self._tol:

                x_vals_below += [new_x]
                p_vals_below += [new_p]
                new_x         = self._nx( new_x - self._k )                     # move one step down and calculate...
                new_p 		  = self._CK(new_x, prev_x_vals, t, prev_p_vals)    # ... probability to reach there
                new_p         = max(0,new_p)                                    # to avoid numerical errors


            # we start from x_start+k, and move upwards in steps of k until x_max(i) is reached
            ############################################################################################################
            x_vals_above	= []
            p_vals_above    = []
            new_x 			= self._nx( x_start + self._k )                     # move one step up and calculate...
            new_p 			= self._CK(new_x, prev_x_vals, t, prev_p_vals)      # ... probability to reach there

            while new_p > self._tol:

                x_vals_above += [new_x]
                p_vals_above += [new_p]
                new_x         = self._nx( new_x + self._k )                     # move one step up and calculate...
                new_p 		  = self._CK(new_x, prev_x_vals, t, prev_p_vals)    # ... probability to reach there
                new_p         = max(0,new_p)                                    # to avoid numerical errors

            # merge the new values above and below x_min(i-1) into one new list and append it to the slices
            ############################################################################################################
            x_vals 		= [ x           for x in np.flipud(x_vals_below) ] + x_vals_above
            x_vals      = [ self._nx(x) for x in x_vals ]
            p_vals 		= [ p           for p in np.flipud(p_vals_below) ] + p_vals_above

            self._xs   += [ x_vals ] # x-values for first time slice
            self._ps   += [ p_vals ] # probailities for first time slice   

        self._reshape_results() # reshape the results such that they can be returned to the output
        self._calculated = True # because else 'get_density' cannot be called


    def get_density(self, x='all', t='all', functional=True ):
        """
        Return the calculated density p(x,t) in a desired format.

        input:
        -----
        x:			Position x at which we want the density. Either numeric value or 'all'. I 'all', the entire
                    range of x-values is returned for a specific time t, either as an array, if functional=False and
                    else as a function. If x is a value, the probability at that value is returned. 

        t:			Time t at which we want the density. Either numeric value 'all' or 'last'. If 'all', the entire 
                    range of t-values is returned, either as an array, if functional=False and else as a function. If t 
                    is a value, the probability at that time is returned. Setting t='all' is equivalent to t=T. 

        functional:	If True, return an interpolated functional object p is returned, which is a function of x, t
                    or both, depending on the arguments for x and t. If False, a value, Pandas Series or DataFrame
                    of discretized p values is returned, depending on the arguments for x and t.
        """

        assert(self._calculated),"Run 'track_density' before calling this method."
        if x is not 'all': assert(isinstance(x,Number)),"x must be number or 'all'."
        if t=='last': t = self._T
        if t is not 'all': assert(isinstance(t,Number)),"t must be number, 'all' or 'last'."
        if t is not 'all' and t < 0: raise ValueError("t < 0 not allowed")

        if x=='all' and t=='all':
            if functional:  return self._probs_functional
            else:           return self._probs_discrete

        elif x=='all':
            t_nearest = self._find_nearest(self._ts, t)
            slice     = self._probs_discrete.loc[t_nearest,:]
            x         = slice.index
            y         = slice.values
            special   = functional and (t_nearest==0) and ( self._x0 is not None ) # very special case
            if special:         return self._interpolate_Dirac()
            elif functional:    return interp1d(x,y, bounds_error=False, fill_value=0, kind='linear')
            else:               return slice

        elif t=='all':
            x_nearest = self._find_nearest(self._probs_discrete.columns, x)
            slice     = self._probs_discrete.loc[:,x_nearest]
            x         = slice.index
            y         = slice.values
            if functional: return interp1d(x,y, bounds_error=False, fill_value=0, kind='linear')
            else:          return slice

        else:
            val = float( self._probs_functional(x,t) )
            return val
