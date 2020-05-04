"""decomposition.py -- Analysis tools for Lanczos decompositions

Language: Python 3
Patrick J. Fasano
University of Notre Dame

+ 10/24/19 (pjf): Created, migrated from diagonalize_alphabeta.py.
+ 02/20/20 (mac): Provide eigenvalue input and overhaul to handle degenerate labels.
"""

import numpy as np
import scipy.linalg as linalg
from . import histogram

import mcscript.utils  # for value_range

def read_lanczos(filename="mfdn_alphabeta.dat", **kwargs):
    """Read and parse mfdn_alphabeta file into d and e arrays for eigh_tridiagonal

    Arguments:
        filename (str, default "mfdn_alphabeta.dat"): input filename
    Returns:
        (tuple of np.array): diagonal and off-diagonal matrix elements
    """

    alphabeta = np.loadtxt(filename, usecols=(1, 2))
    return (alphabeta[:, 0], alphabeta[:-1, 1])

def generate_decomposition_LEGACY(labels, expected_eigenvalues, alphabeta, **kwargs):
    """Generate decomposition from Lanczos alpha-beta matrix and expected eigenvalues.

    This is the initial implementation which returns "low level" decomposition,
    still by eigenvalues, as well as histogram.BinMapping object.

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

def generate_decomposition(alphabeta,eigenvalue_label_dict,verbose=False):
    """Generate decomposition from Lanczos alpha-beta matrix and expected eigenvalues.

    Arguments:
        alphabeta (tuple): (alpha,beta)
            alpha (np.array of float): alpha matrix elements
            beta (np.array of float): beta matrix elements
        eigenvalue_label_dict (dict): mapping from eigenvalue to labels

    Returns:
        decomposition (dict): probabilities binned by labels (as lists of degenerate labels)

    """

    # generate Lanczos decomposition
    alpha, beta = alphabeta
    eigvals, eigvecs = linalg.eigh_tridiagonal(alpha, beta)
    raw_decomposition = [
        (eigval,eigvecs[0, i]**2)
        for i, eigval in enumerate(eigvals)
    ]
    
    # bin Lanczos decomposition
    expected_eigenvalues = sorted(eigenvalue_label_dict.keys())
    label_groups = []
    for eigenvalue in expected_eigenvalues:
        label_groups.append(eigenvalue_label_dict[eigenvalue])
        
    bins = histogram.BinMapping.create_bisection_bins(expected_eigenvalues)
    binned_decomposition = histogram.BinMapping(keys=label_groups, bins=bins)
    for eigval, probability in raw_decomposition:
        binned_decomposition[eigval] += probability
        
    # convert to dict for well-behaved access using label (rather than eigenvalue) as key
    decomposition = binned_decomposition.as_dict()
        
    return decomposition

def eigenvalue_label_dict_am(am_max,verbose=False):
    """Generate eigenvalue dictionary for decomposition by angular momentum Casimir eigenvalue.

    Arguments:
        am_max (float): maximum expected angular momentum (integer or half-integer)

    Returns:
        (dict): mapping from eigenvalue to J (as float)

    """

    eigenvalue_label_dict = {
        float(J*(J+1)): float(J)
        for J in mcscript.utils.value_range(am_max%1,am_max,1)
        }
    return eigenvalue_label_dict

def eigenvalue_label_dict_Nex(Nmax,verbose=False):
    """Generate eigenvalue dictionary for decomposition by Nex.

    Arguments:
        Nmax (int): Nmax for calculation

    Returns:
        (dict): mapping from eigenvalue to Nex (as int)

    """

    eigenvalue_label_dict = {
        float(Nex): Nex
        for Nex in mcscript.utils.value_range(Nmax%2,Nmax,2)
        }
    return eigenvalue_label_dict

def read_eigenvalues(filename,verbose=False):
    """Read table mapping irrep labels to Casimir eigenvalues.

    Eigenvalue degeneracies are allowed.

    File format:
        label1 label2 ... eigenvalue

    Arguments:
        filename (str): Input filename

    Returns:
        (dict): Mapping from eigenvalue to list of degenerate labels
    """

    # read data
    table = np.loadtxt(filename)
    label_length = table.shape[1]-1  # all but last column constitutes label
    if (verbose):
        print(table)
        
    # collect labels by eigenvalue
    eigenvalue_label_dict = {}
    for row in table:
        label = tuple(row[:-1])
        eigenvalue = row[-1]
        eigenvalue_label_dict.setdefault(eigenvalue,[])
        eigenvalue_label_dict[eigenvalue].append(label)

    # convert list of labels to tuples
    #
    # This is necessary so that they can be used as an immutable binning key.
    for eigenvalue in eigenvalue_label_dict:
        eigenvalue_label_dict[eigenvalue] = tuple(eigenvalue_label_dict[eigenvalue])
    
    # diagnostic output
    if (verbose):
        for (eigenvalue,label_list) in eigenvalue_label_dict.items():
            print("{:12e} {:1d} {}".format(eigenvalue,len(label_list),label_list))

    return eigenvalue_label_dict

def mean_am_sqr(decomposition):
    """Calculate mean angular momentum squared from angular momentum distribution.

    May be used for check against MFDn observable mean am sqr.

    May use with mfdnres.tools.effective_am to extract effective am by:

        J*(J+1) = <J^2>

    Arguments:
        decomposition (dict): mapping from am to probability

    Returns
        (float): mean am squared

    """
    mean_am2=0.
    for am in sorted(list(decomposition.keys())):
        mean_am2+=am*(am+1)*decomposition[am]
    return mean_am2

def print_decomposition(decomposition,label_format="",probability_format="8.6f"):
    """Print diagnostic output of decomposition, sorted by labels.

    Note: Must use label_format="", not label_format="s", for tuple labels (or
    else get error "unsupported format string passed to tuple.__format__").

    """

    format_str = "{{:{}}} {{:{}}}".format(label_format,probability_format)
    for labels in sorted(decomposition.keys()):
        print(format_str.format(labels,decomposition[labels]))
    
