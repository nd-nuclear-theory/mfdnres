""" extrapolation.py

    Extrapolation analysis.

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    07/11/17 (mac): Split out from analysis.py.
    09/18/18 (mac): Clean up imports.
"""

import os
import math
import glob

import numpy as np

# intra-package references
from . import nonlinear

################################################################
# extrapolation
################################################################

def extrapolate_energies_exp(Nmax_values,E_values,c2_guess=0.3,verbose=False):
    """Obtain exponentially extrapolated energy.

    Args:
        Nmax_values (NumPy vector): vector of Nmax values for the fit
        E_values (NumPy vector): vector of energy values for the fit
        c2_guess (float, optional): assumed exponential decay constant for initial linear fit
        verbose (bool, optional): whether or not to display fit progress

    Returns:
        (float): extrapolated energy value

    Example:

        >>> Nmax_values = np.array([6,8,10])
        >>> E_values = np.array([-34.849,-37.293,-38.537])
        >>> mfdnres.analysis.extrapolate_energies_exp(Nmax_values,E_values,verbose=True)

        Exponential fit
        Nmax values: [ 6  8 10]
        E values: [-34.849 -37.293 -38.537]
        Linear fit assuming c2 = 0.3
        (c0,c1): (-40.15798483275664, 32.030180719893316)
        Nonlinear fit returns: (array([-39.82661333,  37.74536425,   0.33765202]), 2)
    """

    if (verbose):
        print("Exponential fit")
        print("Nmax values:",Nmax_values)
        print("E values:",E_values)

    # short circuit on missing input values
    for E in E_values:
        if (E is np.nan):
            return np.nan

    # do preliminary linear fit
    # based on typical assumed exponential decay constant
    A = np.array([
        [1,math.exp(-c2_guess*Nmax)]
        for Nmax in Nmax_values
    ])
    (c0_guess,c1_guess) = np.linalg.lstsq(A,E_values)[0]
    if (verbose):
        print("Linear fit assuming c2 = {}".format(c2_guess))
        print("(c0,c1):",(c0_guess,c1_guess))

    # do nonlinear fit
    c_guess = np.array([c0_guess,c1_guess,c2_guess])
    fit = nonlinear.fit(nonlinear.model_exp,Nmax_values,E_values,c_guess)
    if (verbose):
        print("Nonlinear fit returns:",fit)
    (c_values,_) = fit
    (c0,_,_) = c_values

    return c0

################################################################
# test control code
################################################################
    
def extrapolate_energies_exp_test():
    """ Test of exponentially extrapolated energy.
    """

    Nmax_values = np.array([6,8,10])
    E_values = np.array([-34.849,-37.293,-38.537])
    E_extrap = extrapolate_energies_exp(Nmax_values,E_values,verbose=True)
    print("Extrapolated energy: {}".format(E_extrap))

################################################################
# main
################################################################

if (__name__ == "__main__"):
    extrapolate_energies_exp_test()
