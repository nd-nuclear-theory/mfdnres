"""decomposition.py -- Analysis tools for Lanczos decompositions

Language: Python 3
Patrick J. Fasano
University of Notre Dame

+ 10/24/19 (pjf): Created, migrated from diagonalize_alphabeta.py.
"""

import numpy as np
import scipy.linalg as linalg
from . import histogram

def read_lanczos(filename="mfdn_alphabeta.dat", **kwargs):
    """Read and parse mfdn_alphabeta file into d and e arrays for eigh_tridiagonal

    Arguments:
        filename (str, default "mfdn_alphabeta.dat"): input filename
    Returns:
        (tuple of np.array): diagonal and off-diagonal matrix elements
    """

    alphabeta = np.loadtxt(filename, usecols=(1, 2))
    return (alphabeta[:, 0], alphabeta[:-1, 1])


def generate_decomposition(labels, casimir_eigenvalues, alpha, beta, **kwargs):
    """Generate decomposition from Lanczos alpha-beta matrix and expected eigenvalues.

    Arguments:
        labels (list): irrep labels
        casimir_eigenvalues: eigenvalues of Casimir (or linear combination thereof)
        alpha (np.array of float): alpha matrix elements
        beta (np.array of float): beta matrix elements

    """
    eigvals, eigvecs = linalg.eigh_tridiagonal(alpha, beta)

    bins = histogram.BinMapping.create_bisection_bins(casimir_eigenvalues)
    hist = histogram.BinMapping(keys=labels, bins=bins)
    for i, eigval in enumerate(eigvals):
        hist[eigval] += eigvecs[0, i]**2

    return hist