"""nonlinear.py -- nonlinear least squares fitting with SciPy

    Tests of using scipy.optimize.leastsq interface, and wrapper
    functions for exponential fit.

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    6/17/15 (mac): Initiated.
    Last modified 6/17/15.

"""

import numpy as np
import scipy.optimize

################################################################
# model functions
################################################################

def model_line(x_values,params):
    """ Model function for straight line fit.

    For fitting data as y = c0 + c1*x.

    Intended as sanity check case.

    Args:
        params (iterable): parameters [c0,c1]
        x_values (NumPy vector): independent variable values [x0,...,x_(M-1)]

    Returns:
        (NumPy vector): values [f0,...,f_(M-1)]
    """

    (c0,c1) = params
    x_vec = np.array(x_values,dtype=float)  # force to numpy float array to avoid type error surprises
    y_values = c0 + c1*float(x_vec)
    return y_values

def model_exp(x_values,params):
    """ Model function for exponential decay fit with baseline.

    For fitting data as y = c0 + c1*exp(-c2*x).

    Args:
        params (iterable): parameters [c0,c1,c2]
        x_values (NumPy vector): independent variable values [x0,...,x_(M-1)]

    Returns:
        (NumPy vector): values [f0,...,f_(M-1)]
    """

    (c0,c1,c2) = params

    # Debugging: Must force x_values to float.  If use x_values argument inside
    # exponential, this may come as a list of int, which, after np broadcast,
    # gives rise to a type error:
    #
    # TypeError: 'numpy.float64' object cannot be interpreted as an integer
    x_vec = np.array(x_values,dtype=float)  # force to numpy float array to avoid type error surprises
    y_values = c0 + c1*np.exp(-c2*x_vec)
    return y_values

################################################################
# residuals wrapper for least squares fit models
################################################################

def residuals(params,f_model,x_values,y_values):
    """Residual function for nonlinear fit.

    Calculates residuals d_i=y_i-f(x_i), as required for
    scipy.optimize.leastsq.  Typical signature for call will be

        fit = scipy.optimize.leastsq(residuals, c_guess, args=(f_model,x_values,y_values))

    where c_guess are the initial guesses for the parameters.

    The model function must "broadcast" over a vector of x arguments
    but need not be a full-fledged NumPy ufunc.

    Args:
        f_model (callable): model function f_model(x_values,params)
        params (NumPy vector): parameters [c0,c1]
        x_values (NumPy vector): independent variable values [x0,...,x_(M-1)]
        y_values (NumPy vector): dependent variable values [x0,...,x_(M-1)]

     Returns:
        (NumPy vector): residuals [d0,...,d_(M-1)]

    """

    d_values = y_values - f_model(x_values,params)
    return d_values

################################################################
# full wrapper for nonlinear fit
################################################################

def fit(f_model,x_values,y_values,c_guess):
    """Wrapper for scipy.optimize.leastsq least squares data fit.

    Args:
        f_model (callable): model function f_model(x_values,params)
        x_values (NumPy vector): independent variable values [x0,...,x_(M-1)]
        y_values (NumPy vector): dependent variable values [x0,...,x_(M-1)]
        c_guess (NumPy vector): initial guess for parameters [c0,c1,...]

     Returns:
        (tuple): return value from scipy.optimize.leastsq, typically 
             a tuple (params,success_code)

    """

    fit = scipy.optimize.leastsq(residuals, c_guess, args=(f_model,x_values,y_values))
    return fit

################################################################
# test control code
################################################################
    
def line_test():
    """ Test of simple straight-line fit.

    Result: (array([ -3.50005319e-26,   1.10000000e+00]), 2)
    """

    x_values = np.array([0.,1.,2.])
    y_values = np.array([0.,1.1,2.2])
    c_guess = np.array([0.,1.])
    fit = scipy.optimize.leastsq(residuals, c_guess, args=(model_line,x_values,y_values))
    print(fit)

def exp_test():
    """ Test of exponential fit.

    Reproduces test case from berotor2-fig.nb.

    Result: (array([-39.82661333,  37.74536425,   0.33765202]), 2)

    Compare Mathematica c0=-39.8266.
    """

    x_values = np.array([6,8,10])
    y_values = np.array([-34.849,-37.293,-38.537])
    c_guess = np.array([50,50,0.3])
    fit = scipy.optimize.leastsq(residuals, c_guess, args=(model_exp,x_values,y_values))
    print(fit)

################################################################
# main
################################################################

if (__name__ == "__main__"):
    
    line_test()
    exp_test()
