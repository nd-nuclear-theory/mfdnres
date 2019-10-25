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


def generate_decomposition(labels, expected_eigenvalues, alphabeta, **kwargs):
    """Generate decomposition from Lanczos alpha-beta matrix and expected eigenvalues.

    Arguments:
        labels (list): irrep labels
        expected_eigenvalues: expected eigenvalues of Casimir (or linear combination thereof)
        alphabeta (tuple): (alpha,beta)
            alpha (np.array of float): alpha matrix elements
            beta (np.array of float): beta matrix elements

    Returns:
        raw_decomposition (list of tuple): (eigenvalue,probability) pairs from Lanczos alphabeta diagonalization
        binned_decomposition (histogram.BinMapping): probabilities binned by expected labels

    """
    # generate Lanczos decomposition
    alpha, beta = alphabeta
    eigvals, eigvecs = linalg.eigh_tridiagonal(alpha, beta)
    raw_decomposition = [
        (eigval,eigvecs[0, i]**2)
        for i, eigval in enumerate(eigvals)
    ]
    
    # bin Lanczos decomposition
    bins = histogram.BinMapping.create_bisection_bins(expected_eigenvalues)
    binned_decomposition = histogram.BinMapping(keys=labels, bins=bins)
    for eigval, probability in raw_decomposition:
        binned_decomposition[eigval] += probability
        
    return (raw_decomposition,binned_decomposition)
