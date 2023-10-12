"""decomposition.py -- Analysis tools for Lanczos decompositions

Language: Python 3
Patrick J. Fasano
University of Notre Dame

+ 10/24/19 (pjf): Created, migrated from diagonalize_alphabeta.py.
+ 02/20/20 (mac): Provide eigenvalue input and overhaul to handle degenerate labels.
+ 09/02/20 (mac): Add labeled_decomposition and rebinned_decomposition.
+ 09/21/20 (mac): Provide decomposition based on truncated number of Lanczos iterations.
+ 01/25/21 (pjf): Pull raw decomposition generation into its own function.
+ 03/29/21 (zz): Add a function to canonicalize a nuclide and add support for reading canonicalized nuclei files in read_eigenvalues.
+ 05/05/21 (zz): Add rebinning functions for U3LS, Sp3R and Sp3RS.
+ 10/12/23 (mac): Provide decomposition rebinning based on namedtuple label types.

"""

import collections

import numpy as np
import scipy.linalg as linalg
from . import histogram

import mcscript.utils  # for value_range

################################################################
# raw decomposition generation
################################################################

def read_lanczos(filename="mfdn_alphabeta.dat"):
    """Read and parse mfdn_alphabeta file into d and e arrays for eigh_tridiagonal.

    Result of calculation with n Lanczos iterations is an alpha vector of length n and
    beta vector of length n-1.

    Arguments:
        filename (str, default "mfdn_alphabeta.dat"): input filename

    Returns:
        alpha (np.array): vectors of diagonal matrix elements
        beta (np.array): vectors of off-diagonal matrix elements

    """

    # extract raw vectors
    alphabeta_array = np.loadtxt(filename, usecols=(1, 2))
    alpha, beta = (alphabeta_array[:, 0], alphabeta_array[:-1, 1])

    return alpha, beta


def generate_raw_decomposition(alphabeta, lanczos_iterations=None):
    """Generate raw decomposition from Lanczos alpha-beta matrix.

    This does not perform any binning on eigenvalues.

    Arguments:
        alphabeta (tuple): (alpha,beta)
            alpha (np.array of float): alpha matrix elements
            beta (np.array of float): beta matrix elements
        lanczos_iterations (int, optional): number of effective Lanczos iterations to which to truncate

    Returns:
        raw_decomposition (list of tuple): (eigenvalue,probability) pairs from Lanczos alphabeta diagonalization
    """
    # extract matrix elements
    alpha, beta = alphabeta

    # trim vectors
    if (lanczos_iterations is not None):
        (alpha, beta) = (alpha[:lanczos_iterations], beta[:lanczos_iterations-1])

    # generate Lanczos decomposition
    eigvals, eigvecs = linalg.eigh_tridiagonal(alpha, beta)
    raw_decomposition = [
        (eigval,eigvecs[0, i]**2)
        for i, eigval in enumerate(eigvals)
    ]

    return raw_decomposition


################################################################
# decomposition binning
################################################################

def generate_decomposition(alphabeta, eigenvalue_label_dict, lanczos_iterations=None, verbose=False):
    """Generate decomposition from Lanczos alpha-beta matrix and expected eigenvalues.

    Arguments:
        alphabeta (tuple): (alpha,beta)
            alpha (np.array of float): alpha matrix elements
            beta (np.array of float): beta matrix elements
        eigenvalue_label_dict (dict): mapping of eigenvalue to labels (may be tuple of degenerate labels)
            eigenvalue (float)
            labels (int, tuple, etc.)
        lanczos_iterations (int, optional): number of effective Lanczos iterations to which to truncate

    Returns:
        decomposition (dict): probabilities binned by label (given as tuple of degenerate labels)

    """

    # generate Lanczos decomposition
    raw_decomposition = generate_raw_decomposition(alphabeta, lanczos_iterations=lanczos_iterations)
    if (verbose):
        print("Raw decomposition")
        print(np.array(raw_decomposition))

    # bin Lanczos decomposition
    expected_eigenvalues = sorted(eigenvalue_label_dict.keys())
    if (verbose):
        print("Expected eigenvalue -> label group")
        for eigenvalue in expected_eigenvalues:
            print("{:+8.3f} -> {}".format(eigenvalue, eigenvalue_label_dict[eigenvalue]))
    label_groups = []
    for eigenvalue in expected_eigenvalues:
        label_groups.append(eigenvalue_label_dict[eigenvalue])
    bins = histogram.BinMapping.create_bisection_bins(expected_eigenvalues)
    binned_decomposition = histogram.BinMapping(keys=label_groups, bins=bins)
    for eigenvalue, probability in raw_decomposition:
        binned_decomposition[eigenvalue] += probability
    if (verbose):
        print("Binned results (sorted by eigenvalue)")
        for eigenvalue in expected_eigenvalues:
            print("{:+8.3f} : {:8.6f}".format(eigenvalue, binned_decomposition[eigenvalue]))

    # convert to dict for well-behaved access using label (rather than eigenvalue) as key
    decomposition = binned_decomposition.as_dict()

    return decomposition


def eigenvalue_label_dict_am(am_max, verbose=False):
    """Generate eigenvalue dictionary for decomposition by angular momentum Casimir eigenvalue.

    Arguments:
        am_max (float): maximum expected angular momentum (integer or half-integer)

    Returns:
        (dict): mapping from eigenvalue to J (as float)

    """

    eigenvalue_label_dict = {
        float(J*(J+1)): float(J)
        for J in mcscript.utils.value_range(am_max%1, am_max, 1)
        }
    return eigenvalue_label_dict


def eigenvalue_label_dict_Nex(Nmax, verbose=False):
    """Generate eigenvalue dictionary for decomposition by Nex.

    Arguments:
        Nmax (int): Nmax for calculation

    Returns:
        (dict): mapping from eigenvalue to Nex (as int)

    """

    eigenvalue_label_dict = {
        float(Nex): Nex
        for Nex in mcscript.utils.value_range(Nmax%2, Nmax, 2)
        }
    return eigenvalue_label_dict


def read_eigenvalues(filename, swap=False, verbose=False):
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
        if swap:
            label = tuple(list(row[:3])+list(row[4:2:-1])+list(row[5:-1])) # swap row[3] and row[4]
        else:
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
        for (eigenvalue, label_list) in eigenvalue_label_dict.items():
            print("{:12e} {:1d} {}".format(eigenvalue, len(label_list), label_list))

    return eigenvalue_label_dict


################################################################
# decomposition postprocessing
################################################################

def mean_am_sqr(decomposition):
    """Calculate mean angular momentum squared from angular momentum distribution.

    May be used for check against MFDn observable mean am sqr.

    May be used with mfdnres.tools.effective_am to extract effective am by:

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

def print_decomposition(decomposition, label_format="", probability_format="8.6f"):
    """Print diagnostic output of decomposition, sorted by labels.

    Note: Must use label_format="", not label_format="s", for tuple labels (or
    else get error "unsupported format string passed to tuple.__format__").

    Arguments:
        decomposition (dict): mapping from label (typically int or tuple) to probability
        label_format (str,optional): format descriptor for label
        probability_format (str,optional): format descriptor for probability
    """

    format_str = "{{:{}}} {{:{}}}".format(probability_format, label_format)
    for labels in sorted(decomposition.keys()):
        print(format_str.format(decomposition[labels], labels))

        
################################################################
# labeling native decompositions
################################################################

def labeled_decomposition(label_list, decomposition):
    """Digest natively-calculated decomposition into dictionary.

    This function repackaged decompositions provided in the results from the
    many-body code (e.g., the decomposition by Nex from mfdn).

    Dictionary is of form:

        (label,)->probability

    Although the native decompositions have no label degeneracies, the label
    must still be wrapped in a tuple for compatibility with analyses of
    decompositions where multiple labels may be binned together (e.g., from
    Lanczos decompositions where different labels have degenerat eigenvalues).

    Arguments:
        label_list (list): list of individual bin labels (typically int or tuple)
        decomposition (np.array): one-dimensional array of probabilities

    """

    if (len(label_list) != len(decomposition)):
        raise ValueError("mismatched lengths for label_list and decomposition")

    decomposition_dict = {
        (label, ) : probability
        for (label, probability) in zip(label_list, decomposition)
    }

    return decomposition_dict


################################################################
# decomposition rebinning -- legacy subsetting functions
################################################################

# transformation functions to extract coarse-grained labeling from fine-grained labeling

# DEPRECATED in favor of using label_subsetting_function factory function

def label_transformation_u3spsns_to_nex(labels):
    """
    """
    (N,lam,mu,Sp,Sn,S) = labels
    return N

def label_transformation_u3spsns_to_s(labels):
    """
    """
    (N,lam,mu,Sp,Sn,S) = labels
    return S

def label_transformation_u3spsns_to_u3(labels):
    """
    """
    (N,lam,mu,Sp,Sn,S) = labels
    return (N,lam,mu)

def label_transformation_u3spsns_to_u3s(labels):
    """
    """
    (N,lam,mu,Sp,Sn,S) = labels
    return (N,lam,mu,S)

def label_transformation_u3lspsns_to_nex(labels):
    """
    """
    (N,lam,mu,Sp,Sn,S,L) = labels
    return N

def label_transformation_u3lspsns_to_s(labels):
    """
    """
    (N,lam,mu,Sp,Sn,S,L) = labels
    return S

def label_transformation_u3lspsns_to_l(labels):
    """
    """
    (N,lam,mu,Sp,Sn,S,L) = labels
    return L

def label_transformation_u3lspsns_to_u3ls(labels):
    """
    """
    (N,lam,mu,Sp,Sn,S,L) = labels
    return (N,lam,mu,S,L)

def label_transformation_baby_spncci_to_s(labels):
    """
    """
    (N_sigma,lambda_sigma,mu_sigma,N_omega,lambda_omega,mu_omega,Sp,Sn,S) = labels
    return S

def label_transformation_baby_spncci_to_nex(labels):
    """
    """
    (N_sigma,lambda_sigma,mu_sigma,N_omega,lambda_omega,mu_omega,Sp,Sn,S) = labels
    return N_omega

def label_transformation_baby_spncci_to_u3(labels):
    """
    """
    (N_sigma,lambda_sigma,mu_sigma,N_omega,lambda_omega,mu_omega,Sp,Sn,S) = labels
    return (N_omega,lambda_omega,mu_omega)

def label_transformation_baby_spncci_to_u3s(labels):
    """
    """
    (N_sigma,lambda_sigma,mu_sigma,N_omega,lambda_omega,mu_omega,Sp,Sn,S) = labels
    return (N_omega,lambda_omega,mu_omega,S)

def label_transformation_baby_spncci_to_sp3r(labels):
    """
    """
    (N_sigma,lambda_sigma,mu_sigma,N_omega,lambda_omega,mu_omega,Sp,Sn,S) = labels
    return (N_sigma,lambda_sigma,mu_sigma)

def label_transformation_baby_spncci_to_sp3rs(labels):
    """
    """
    (N_sigma,lambda_sigma,mu_sigma,N_omega,lambda_omega,mu_omega,Sp,Sn,S) = labels
    return (N_sigma,lambda_sigma,mu_sigma,S)

################################################################
# decomposition rebinning -- label subsetting
################################################################

# singlet labels
NexLabels = collections.namedtuple("NexLabels", ["N_omega"])
SLabels = collections.namedtuple("SLabels", ["S"])
LLabels = collections.namedtuple("LLabels", ["L"])

# U(3) but non-Sp(3,R) labels
U3Labels = collections.namedtuple("U3Labels", ["N_omega", "lambda_omega", "mu_omega"])
U3SLabels = collections.namedtuple("U3SLabels", ["N_omega", "lambda_omega", "mu_omega", "S"])
U3LSpSnSLabels = collections.namedtuple("U3LSpSnSLabels", ["N_omega", "lambda_omega", "mu_omega", "Sp", "Sn", "S", "L"])
U3SpSnSLabels = collections.namedtuple("U3SpSnSLabels", ["N_omega", "lambda_omega", "mu_omega", "Sp", "Sn", "S"])

# Sp(3,R) labels
Sp3RLabels = collections.namedtuple("Sp3RLabels", ["N_sigma", "lambda_sigma", "mu_sigma"])
Sp3RSLabels = collections.namedtuple("Sp3RSLabels", ["N_sigma", "lambda_sigma", "mu_sigma", "S"])
BabySpNCCILabels = collections.namedtuple("BabySpNCCILabels", ["N_sigma", "lambda_sigma", "mu_sigma", "N_omega", "lambda_omega", "mu_omega", "Sp", "Sn", "S"])

def namedtuple_subsetting_function(long_labels_type, short_labels_type):
    """Factory function to provide subsetting function between namedtuple types.

    Assumes fields in short_labels_type are subset of those in long_labels_type.
    Otherwise, use of resulting subsetting function will result in a TypeError
    exception.

    The subsetting function produced here is based on dictionary operations.  It
    could perhaps be made much more efficient if the mapping between long and
    short tuple fields were set up in terms of positional arguments at function
    production time.

    Example:

        >>> import collections
        >>> LongLabels = collections.namedtuple("LongLabels", ["a", "b", "c"])
        >>> ShortLabels = collections.namedtuple("ShortLabels", ["a", "c"])
        >>> long_labels = LongLabels(1,2,3)
        >>> f = mfdnres.decomposition.namedtuple_subsetting_function(LongLabels, ShortLabels)
        >>> short_labels = f(long_labels)
        >>> print("{} -> {}".format(long_labels, short_labels))
        
        LongLabels(a=1, b=2, c=3) -> ShortLabels(a=1, c=3)

    Arguments:

        long_labels_type (collections.namedtuple): long label tuple type

        short_labels_type (collections.namedtuple): short label tuple type

    Return:

        (callable): subsetting function (tuple or long_labels_type -> short_labels_type)

    """

    def the_subsetting_function(long_labels):
        # cast argument (which might just be plain tuple) to long_labels_type namedtuple
        long_labels = long_labels_type(*long_labels)

        # subset the key-value pairs from long_labels to those supported by short_labels_type
        filtered_dict = {
            k: v
            for k, v in long_labels._asdict().items()
            if k in short_labels_type._fields
        }

        # cast result to short_labels_type namedtuple
        short_labels = short_labels_type(**filtered_dict)
        
        return short_labels
    
    return the_subsetting_function

def rebinned_decomposition(decomposition, subsetting_specifier, verbose=False):
    """Rebin decomposition according to new labeling.

    E.g., may by used to rebin "U3S" to S, by transformation
    label_transformation_U3S_to_S.

    Arguments:

        decomposition (dict): mapping from tuple of degenerate labels label (typically int or tuple) to probability

        subsetting_specifier (tuple or callable): tuple (old_label_type,
        new_label_type) giving old (long) and new (short) label types; for
        legacy support, may be function mapping old label to new label

    Returns
        (dict): rebinned decomposition

    """

    if type(subsetting_specifier) is tuple:
        old_label_type, new_label_type = subsetting_specifier
        subsetting_function = namedtuple_subsetting_function(old_label_type, new_label_type)
    else:
        print("gh")
        subsetting_function = subsetting_specifier
        
    new_decomposition = {}
    for (label_list,probability) in decomposition.items():
        if (verbose):
            print(subsetting_function,label_list,probability)
        new_label_list = tuple(set(map(subsetting_function,label_list)))
        new_decomposition[new_label_list] = new_decomposition.get(new_label_list,0) + probability

    return new_decomposition
