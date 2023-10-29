"""decomposition.py -- Analysis tools for Lanczos decompositions

Language: Python 3
Patrick J. Fasano
University of Notre Dame

    - 10/24/19 (pjf): Created, migrated from diagonalize_alphabeta.py.
    - 02/20/20 (mac): Provide eigenvalue input and overhaul to handle degenerate labels.
    - 09/02/20 (mac): Add labeled_decomposition and rebinned_decomposition.
    - 09/21/20 (mac): Provide decomposition based on truncated number of Lanczos iterations.
    - 01/25/21 (pjf): Pull raw decomposition generation into its own function.
    - 03/29/21 (zz): Add a function to canonicalize a nuclide and add support for reading canonicalized nuclei files in read_eigenvalues.
    - 05/05/21 (zz): Add rebinning functions for U3LS, Sp3R and Sp3RS.
    - 10/12/23 (mac):
        + Provide decomposition rebinning based on namedtuple label types.
        + Add slurp_lanczos_filenames() to populate MFDnResultsData Lanczos decomposition filenames.
    - 10/25/23 (mac): Replace slurp_lanczos_filenames() with slurp_lanczos_files().
    - 10/26/23 (mac): Support use of eigenvalue dict as label list for labeled_decomposition().
    - 10/28/23 (mac):
        + Add eigenvalue filename search function decomposition_eigenvalue_filename().
        + Provide __str__ method to label classes.
        + Add plotting functions set_up_decomposition_axes() and add_decomposition_plot().
    - 10/29/23 (mac):
        + Add canonical sorting (of labels within degenerate set, and of bins by
        their label sets) to rebinned_decomposition().
        + Add filter_decomposition().

"""

import collections
import glob
import os

import numpy as np
import scipy.linalg as linalg
from . import data, histogram, input, mfdn_results_data, ticks

import mcscript.utils  # for value_range

################################################################
# lanczos data input
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
    alpha_beta_array = np.loadtxt(filename, usecols=(1, 2))
    alpha, beta = (alpha_beta_array[:, 0], alpha_beta_array[:-1, 1])

    return alpha, beta


def slurp_lanczos_files(
        directory_list,
        filename_format,
        glob_pattern="*.lanczos",
        verbose=False,
):
    """Read all lanczos files in given directories.

    The results will be a list of results data objects, one for
    each lanczos data set.

    Inspired by mfdnres.input.slurp_res_files().

    Arguments:

        directory_list (str or list of str): directory or list of directories
            containing files to import

        filename_format (str,optional): identifier string for the results
            filename parser to use

        glob_pattern (str,optional): glob pattern for results filenames to read
            within each directory

        verbose (bool,optional): enable debugging output

    Returns:
        (list of MFDnResultsData): list of mesh point data objects

    """

    # process argument: upgrade single directory to list
    if (type(directory_list) == str):
        directory_list = [directory_list]
    directory_list = sorted(list(set(directory_list)))  # remove duplicate input directories
    if (verbose):
        print("  slurp_lanczos_filenames: directory list {}".format(directory_list))

    # accumulate mesh points
    mesh_data = []
    for directory in directory_list:
        full_glob_pattern = os.path.join(directory, glob_pattern)
        if (verbose):
            print("  slurp_lanczos_filenames: searching for files {}...".format(full_glob_pattern))
        filename_list = glob.glob(full_glob_pattern)

        # accumulate parsed data from different res files
        for filename in filename_list:

            # set up container
            results = mfdn_results_data.MFDnResultsData()

            # parse lanczos filename for run parameters
            info_from_filename = input.parse_filename(filename, filename_format)
            results.params.update(info_from_filename)

            # save lanczos filename
            decomposition_type = results.params["decomposition_type"]
            qn = results.params["decomposition_state"]
            alpha, beta = read_lanczos(filename)
            results.mfdn_level_lanczos_decomposition_data = {
                decomposition_type: {
                    qn: (filename, alpha, beta)
                }
            }

            mesh_data.append(results)

    if (verbose):
        print("  slurp_lanczos_filenames: extracted mesh points {}".format(len(mesh_data)))

    return mesh_data


################################################################
# raw decomposition generation
################################################################

def generate_raw_decomposition(alpha_beta, lanczos_iterations=None):
    """Generate raw decomposition from Lanczos alpha-beta matrix.

    This does not perform any binning on eigenvalues.

    Arguments:
        alpha_beta (tuple): (alpha,beta)
            alpha (np.array of float): alpha matrix elements
            beta (np.array of float): beta matrix elements
        lanczos_iterations (int, optional): number of effective Lanczos iterations to which to truncate

    Returns:
        raw_decomposition (list of tuple): (eigenvalue,probability) pairs from Lanczos alphabeta diagonalization
    """
    # extract matrix elements
    alpha, beta = alpha_beta

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

def generate_decomposition(alpha_beta, eigenvalue_label_dict, lanczos_iterations=None, verbose=False):
    """Generate decomposition from Lanczos alpha-beta matrix and expected eigenvalues.

    Arguments:
        alpha_beta (tuple): (alpha,beta)
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
    raw_decomposition = generate_raw_decomposition(alpha_beta, lanczos_iterations=lanczos_iterations)
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

################################################################
# decomposition eigenvalues
################################################################

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


# decomposition eigenvalue/coefficient paths
decomposition_dir_list = [
]

def decomposition_eigenvalue_filename(
        nuclide, Nmax, decomposition_type, *,
        decomposition_base_path=None,
        verbose=False
):
    """Search for full path name of decomposition eigenvalue file.

    Searches under decompositoin_base_path, in subdirectories listed in
    mfdnres.decomposition.decomposition_dir_list.

    Arguments:

        nuclide (tuple): (Z,N)
        
        Nmax (int): Nmax

        decomposition_type (str): decomposition type ("U3SpSnS", etc.)

        decomposition_base_path (str, optional): base path under which subdirectories
        containing eigenvalue files are to be found; if None, defaults to standard path
        ${GROUP_HOME}/data/u3-subspaces/decomposition


    Returns:

        (str): qualified filename

    """

    if decomposition_base_path is None:
        decomposition_base_path = os.path.join(os.environ["GROUP_HOME"], "data", "u3-subspaces", "decomposition")
    eigenvalue_filename_format ="decomposition_Z{nuclide[0]:02d}_N{nuclide[1]:02d}_Nmax{Nmax:02d}_{decomposition_type}_eigenvalues.dat"
    decomposition_types_with_casimir_in_filename = ["SU3", "Sp3R"]
    if decomposition_type in decomposition_types_with_casimir_in_filename: 
        decomposition_type_for_filename = "C{}".format(decomposition_type)
    else:
        decomposition_type_for_filename = decomposition_type
    eigenvalue_filename = eigenvalue_filename_format.format(nuclide=nuclide, Nmax=Nmax, decomposition_type=decomposition_type_for_filename)
    qualified_filename = mcscript.utils.search_in_subdirectories(
        decomposition_base_path,
        decomposition_dir_list,
        eigenvalue_filename,
        error_message="decomposition eigenvalue file not found",
        verbose=verbose,
    )

    return qualified_filename

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
    Lanczos decompositions where different labels have degenerate eigenvalues).

    Example:

        >>> results_data = ...
        >>> Nmax = 4
        >>> state = (1.5,1,1)
        >>> native_decomposition = results_data.get_decomposition("Nex", state)
        >>> eigenvalue_label_dict = mfdnres.decomposition.eigenvalue_label_dict_Nex(Nmax)
        >>> decomposition = mfdnres.decomposition.labeled_decomposition(eigenvalue_label_dict, native_decomposition)

        {(0,): 0.7587, (2,): 0.1309, (4,): 0.1104}
  
    Arguments:

        label_list (list or dict): list of individual bin labels (typically int
        or tuple); if dict, as returned by eigenvalue_label_dict_am() or
        eigenvalue_label_dict_Nex(), the dict values are used

        decomposition (np.array): one-dimensional array of probabilities

    """

    # take list from eigenvalue_label_dict
    if isinstance(label_list, dict):
        label_list = label_list.values()

    # validate number of provided labels
    print(label_list, decomposition)
    if (len(label_list) != len(decomposition)):
        raise ValueError("mismatched lengths for label_list and decomposition")

    # pair labels with values
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
# decomposition label classes
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
Sp3RSpSnSLabels = collections.namedtuple("Sp3RSLabels", ["N_sigma", "lambda_sigma", "mu_sigma", "Sp", "Sn", "S"])
BabySpNCCILabels = collections.namedtuple("BabySpNCCILabels", ["N_sigma", "lambda_sigma", "mu_sigma", "N_omega", "lambda_omega", "mu_omega", "Sp", "Sn", "S"])

# lookup table for decomosition label classes
LABEL_CLASS_BY_DECOMPOSITION_TYPE = {
    "Nex": NexLabels,
    "S": SLabels,
    "L": LLabels,
    "U3": U3Labels,
    "U3S": U3SLabels,
    "U3LSpSnS": U3LSpSnSLabels,
    "U3SpSnS": U3SpSnSLabels,
    "Sp3R": Sp3RLabels,
    "Sp3RS": Sp3RSLabels,
    "Sp3RSpSnS": Sp3RSpSnSLabels,
    "BabySpNCCI": BabySpNCCILabels,
}

# string formatting

def format_u3_label(self):
    N, lam, mu = self
    label_text = "{:d}({:d},{:d})".format(int(N), int(lam), int(mu))
    return label_text

U3Labels.__str__ = format_u3_label
Sp3RLabels.__str__ = format_u3_label

def format_u3s_label(self):
    # TODO 10/28/23 (mac): upgrade spin formatting from float to solidus fraction
    N, lam, mu, S = self
    label_text = "{:d}({:d},{:d}){:s}".format(int(N), int(lam), int(mu), ticks.half_int_str(S))
    return label_text

U3SLabels.__str__ = format_u3s_label
Sp3RLabels.__str__ = format_u3s_label

def format_u3sss_label(self):
    # TODO 10/28/23 (mac): upgrade spin formatting from float to solidus fraction
    N, lam, mu, Sp, Sn, S = self
    label_text = "{:d}({:d},{:d}){:s},{:s},{:s}".format(int(N), int(lam), int(mu), ticks.half_int_str(Sp), ticks.half_int_str(Sn), ticks.half_int_str(S))
    return label_text

U3SpSnSLabels.__str__ = format_u3sss_label
Sp3RSpSnSLabels.__str__ = format_u3sss_label

################################################################
# decomposition rebinning -- label subsetting
################################################################

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

    TODO 10/28/23 (mac): Merge overlapping label lists.

    Arguments:

        decomposition (dict): mapping from tuple of degenerate labels label (typically int or tuple) to probability

        subsetting_specifier (tuple[str]): tuple (old_decomposition_type,
            new_decomposition_type) defining old (longer) and new (shorter)
            label types; for legacy support, may instead be a callable
            (function) mapping old label to new label

    Returns
        (dict): rebinned decomposition

    """

    # construct label subsetting function
    if type(subsetting_specifier) is tuple:
        old_decomposition_type, new_decomposition_type = subsetting_specifier
        old_label_type = LABEL_CLASS_BY_DECOMPOSITION_TYPE[old_decomposition_type]
        new_label_type = LABEL_CLASS_BY_DECOMPOSITION_TYPE[new_decomposition_type]
        subsetting_function = namedtuple_subsetting_function(old_label_type, new_label_type)
    else:
        subsetting_function = subsetting_specifier

    # rebin decomposition
    new_decomposition = {}
    for (label_list, probability) in decomposition.items():
        if (verbose):
            print(subsetting_function, label_list, probability)
        new_label_set = set(map(subsetting_function, label_list))  # downsample and remove redundancies (with set)
        new_label_list = tuple(sorted(tuple(new_label_set)))   # canonicalize order of labels
        new_decomposition[new_label_list] = new_decomposition.get(new_label_list, 0) + probability

    # sort new decomposition canonically by label sets
    new_decomposition = dict(sorted(new_decomposition.items()))

    return new_decomposition


def filter_decomposition(condition, decomposition, verbose=False):
    """Filter decomposition (by labels or value).

    Order of arguments is chosen to match that of built-in filter().

    Arguments:

        condition (callable): filter condition applied to decomposition
            dictionary items (should have signature ((labels_1, labels_2, ...),
            probability) -> bool)

        decomposition (dict): mapping from tuple of degenerate labels label
            (typically int or tuple) to probability

    Returns
        (dict): filtered decomposition

    """

    decomposition = dict(filter(condition, decomposition.items()))
    
    return decomposition


def Nex_filter(Nex_max):
    """ Generate filter condition to select by Nex of first decomposition label.

    For use with filter_decomposition().
    """
    
    return lambda decomposition_item : decomposition_item[0][0][0]<=Nex_max

################################################################
# decomposition plotting
################################################################

def set_up_decomposition_axes(
        ax, decomposition, *,
        decomposition_range_extension=(0.02,0.02),
        decomposition_axis_label_text=None,
        probability_range=(0.00,1.00),
        probability_range_extension=(0.00,0.00),
        probability_scale=None,
        probability_axis_label_text=None,
        decomposition_labelpad=None,
        probability_labelpad=None,
        probability_tick_specifier=None,
        labelrotation=0,
        verbose=False,
):
    """ Set up axes for probability decomposition plot.

    Arguments:

        ax (mpl.axes.Axes): axes object

        decomposition (dict): decomposition dictionary (with label sets as keys)
 
        decomposition_axis_label_text (str, optional): override for decomposition axis label text

        probability_range (tuple of float): y range, or None for matplotlib auto

        probability_scale (str): y scale ("linear" or "log")

        probability_axis_label_text (str, optional): override for probability axis label text

        decomposition_range_extension (tuple of float, optional): x range relative extension

        decomposition_range_extension (tuple of float, optional): y range relative extension

        decomposition_labelpad (scalar, optional): pass-though labelpad option for xlabel

        probability_labelpad (scalar, optional): pass-though labelpad option for ylabel

        probability_tick_specifier (tuple, optional): tick specification (min,max,step,num_subdivision) for probability ticks

        labelrotation (int): pass through option for ax.tick_params for x tick rotation

    """

    # construct labels
    def format_label_set(label_set):
        label_text = " / ".join(map(str, label_set))
        return label_text
    decomposition_labels = list(map(format_label_set, decomposition.keys()))
    num_labels = len(decomposition_labels)
    if verbose:
        print("Labels ({}): {}".format(num_labels, decomposition_labels))
    
    # set ticks
    #
    # Note that the tick specification must come *before* setting range limits,
    # since the range limits automatically readjust when you set the ticks.
    ax.set_xticks(range(num_labels), decomposition_labels)
    ax.tick_params(axis="x", bottom=False, labelrotation=labelrotation)
    
    if probability_tick_specifier is not None:
        y_ticks = ticks.linear_ticks(*probability_tick_specifier)
        ticks.set_ticks(ax,"y",y_ticks)

    # set limits
    decomposition_range = (-0.5, num_labels-0.5)  # histogram bar rendering range (for width=1.0)
    ax.set_xlim(*data.extend_interval_relative(decomposition_range, decomposition_range_extension))
    if probability_scale=="log":
        # Note: Override any range extension for log scale
        probability_range_extension=(0.,0.)
        ax.set_yscale("log")
    if (probability_range is not None) and np.isfinite(probability_range[0]).all():
        ax.set_ylim(*data.extend_interval_relative(probability_range, probability_range_extension))
        
    # set axis labels
    if probability_axis_label_text is not None:
        ax.set_xlabel(decomposition_axis_label_text, labelpad=decomposition_labelpad)
    if probability_axis_label_text is None:
        probability_axis_label_text = "P"
    ax.set_ylabel(
        r"${}$".format(probability_axis_label_text),
        labelpad=probability_labelpad,
    )


def add_decomposition_plot(
        ax, decomposition,
        **kwargs,
):
    """Add decomposition bar plot to axes.

    Arguments:

        ax (mpl.axes.Axes): axes object

        decomposition (dict): decomposition dictionary (with label sets as keys)

        kwargs (Line2D properties, optional): kwargs are used to specify line
        properties

    """

    kw_defaults = {
        "width": 1.0,
        "color": "white",
        "edgecolor": "black",
    }

    kw_full = {
        **kw_defaults,
        **kwargs,
    }
    
    ax.bar(
        range(len(decomposition.values())),
        decomposition.values(), 
        **kw_full,
    )

    
