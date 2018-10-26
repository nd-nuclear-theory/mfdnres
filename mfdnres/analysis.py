""" analysis.py -- analysis tools for processing MFDn results file

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    06/02/15 (mac): Initiated (as mfdn_analysis.py).
    06/05/15 (mac): Restructure as subpackage.
    07/07/17 (mac): Remove obsolete import_res_files.
    07/08/17 (mac): Start adding new tabulation functions which return
        arrays of data rather than writing to file.
    07/11/17 (mac): Pull out extrapolation function to remove scipy dependence.
    07/30/17 (mac): Pull out band analysis functions to band.py.
    06/01/18 (pjf): Add reference state option for energy table generation.
    09/05/18 (mac): Remove deprecated write_xxx_table functions.  Update dosctrings.
    09/08/18 (mac):
        - Add floor2 arithmetic function for odd-even truncation comparisons.
        - Add energy difference tabulation function.
    09/18/18 (mac): Add dictionary helper functions common_keys and operate_dictionaries.
    09/29/18 (mac): Add sorting and pruning of nan values to energy tabulation.

"""

import math
import functools

import numpy as np

################################################################
# global constants
################################################################

# electromagnetic g factors (glp,gln,gsp,gsn)
MU_EM = np.array([1,0,5.586,-3.826])

################################################################
# arithmetic helper functions
################################################################

def floor2(i):
    """Find greatest even integer less than or equal to given integer.

    Arguments:
        i (int): starting integer

    Returns:
        (int): greatest even integer
    """
    
    return i - (i%2)

################################################################
# dictionary helper functions
################################################################

def common_keys(dictionary_list):
    """ Identify keys common to a set of dictionaries.

    Arguments:
        dictionary_list (list of dict): dictionaries to analyze

    Returns:
        (list): sorted list of common keys
    """

    # find intersection of key sets
    common_key_set = None
    for current_dictionary in dictionary_list:
        current_key_set = set(current_dictionary.keys())
        if (common_key_set is None):
            # prime the intersection with initial key set
            common_key_set = current_key_set
        else:
            # find intersection with current key set
            common_key_set &= current_key_set

    # convert key set into a sorted list
    common_key_list = sorted(list(common_key_set))
    ## print("Common keys: {}".format(common_key_list))

    return common_key_list

def operate_dictionaries(dict1,dict2,op):
    """ Perform binary operation on common entries of two dictionaries.

    Note: Common operations are provided by the operator library.

    Arguments:
        dict1 (dict): dictionary providing operand 1
        dict2 (dict): dictionary providing operand 2
        op (callable): binary operator

    Returns:
        (dict): results dictionary
    """
    
    results = dict()
    for key in common_keys((dict1,dict2)):
        results[key] = op(dict1[key],dict2[key])
    return results

################################################################
# consolidate data into dictionary
################################################################

def extract_key(key_fields,mesh_point):
    """Generate key tuple from mesh point data's parameters.

    This function may be specialized to return a specific key, e.g.:

        key_function_Nsigmamax_Nmax_hw = functools.partial(
            extract_key,
            ("Nsigmamax","Nmax","hw")
        )

    This process would normally be done through make_key_function.

    Arguments:
        key_fields (tuple of str) : tuple of run descriptor keys to
            use to generate key to results
        mesh_data (list of ResultsData): data for mesh points

    Returns:
        (tuple): values of given parameters

    """
    key = tuple([mesh_point.params[key] for key in key_fields])
    return key

def make_key_function(key_descriptor):
    """Generate key function from dtype descriptor for parameter tuple.

    Example key descriptor:

       (("Nsigmamax",int),("Nmax",int),("hw",float))

    Arguments:
        key_descriptor (tuple of tuple): dtype descriptor for key

    """

    key_fields = [name for (name,_) in key_descriptor]
    key_function = functools.partial(extract_key,key_fields)
    return key_function

def make_results_dict(
        mesh_data,key_descriptor,
        verbose=False
):
    """Load mesh data into dictionary, using specified parameter tuple as
       key.

    Example key descriptor:

       (("Nsigmamax",int),("Nmax",int),("hw",float))

    Example:
        >>> KEY_DESCRIPTOR_NMAX_HW = (("Nmax",int),("hw",float))


    For now, in the even that the same mesh point arises multiple
    times on input (i.e., a given value for the key tuple is
    duplicated), the final occurrence overwrites any earlier
    occurrences.  In the future, a more sophisticated "merging"
    process might be appropriate.

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        verbose (bool,optional): verbose output

    Returns:
        (dict): mapping from key tuple to data object

    """

    key_function = make_key_function(key_descriptor)

    results_dict = dict()
    for mesh_point in mesh_data:

        # make key
        key = key_function(mesh_point)
        if (verbose):
            print("  make_results_dict: filename {} key {}".format(mesh_point.filename,key))

        # store data point
        if (key not in results_dict):
            # save mesh point
            results_dict[key] = mesh_point
        else:
            # TODO: do smart merge "update" on existing mesh point
            # overwrite mesh point
            results_dict[key] = mesh_point

    return results_dict

def selected_mesh_data(
        mesh_data,selection_dict,
        verbose=False
):
    """Select mesh points with specific values for given paremeters.

    Example key-value specification:

        {"nuclide":(3,4),"interaction":"JISP16"}

    which is then interpreted as

       (("nuclide",(3,4)),("interaction","JISP16"))

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        selection_dict (dict): key value pairs for parameter values to select
        verbose (bool,optional): verbose output

    Returns:
        (list of ResultsData): selected data for mesh points

    """

    # set up key extraction and selection
    key_value_pairs = selection_dict.items()  # convert to view as iterable of pairs for cleaner debugging
    key_function = make_key_function(key_value_pairs)
    selected_values = tuple([value for (_,value) in key_value_pairs])
    if (verbose):
        print("  selected_mesh_data: selection key-value pairs {}".format(tuple(key_value_pairs)))

    # select submesh
    new_mesh_data = [
        mesh_point
        for mesh_point in mesh_data
        if (key_function(mesh_point) == selected_values)
    ]
    if (verbose):
        print("  selected_mesh_data: selected mesh points {}".format(len(new_mesh_data)))
    return new_mesh_data

def sorted_mesh_data(
        mesh_data,key_descriptor,
        verbose=False
):
    """Sort and consolidate mesh data, using specified parameter tuple as
       key.

    See docstring for make_results_dict, which provides the underlying
    sorting and consolidation engine.

    Example:
        # import data
        slurp_directory = mfdnres.res.res_file_directory("mcaprio","spncci","mac0419")
        mesh_data = mfdnres.res.slurp_res_files(slurp_directory,"spncci",verbose=False)
        mesh_data = mfdnres.analysis.sorted_mesh_data(mesh_data,mfdnres.spncci.KEY_DESCRIPTOR_NNHW)

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        verbose (bool,optional): verbose output

    Returns:
        (list): sorted and consolidated list of data objects
    """

    results_dict = make_results_dict(
        mesh_data,key_descriptor,
        verbose=verbose
    )
    new_mesh_data = [
        results_dict[key]
        for key in sorted(results_dict.keys())
        ]
    return new_mesh_data

################################################################
# tabulation functions -- observable vs. parameters
################################################################

def make_energy_table(mesh_data,key_descriptor,qn,key_list=None,prune=False):
    """Generate energy tabulation.

    The key descriptor is used to provide the parameter columns of the
    tabulation (and may in general be different from the key which was
    used for sorting the data, e.g., it might add additional parameter
    columns).

    It is expected that the mesh data have already been consolidated to make the
    key values unique.  Only key tuples which are present in both sets of mesh
    data are retained.

    Results are then sorted by key tuple.

    Data format:
        param1 param2 ... E

    Arguments:
        mesh_data (list of ResultsData): data set
        key_descriptor (tuple of tuple): dtype descriptor for key
        qn (tuple): quantum number (J,g,n) of level of interest
        key_list (list of tuple, optional): key list for data points; defaults to
            keys for which qn found
        prune (bool, optional): whether or not to prune nan values from table

    Returns:
       (array): data table

    """

    # process results into dictionaries
    results_dict = make_results_dict(mesh_data,key_descriptor)

    # find common keys
    if (key_list is not None):
        common_key_list = key_list
    else:
        common_key_list = common_keys([results_dict])
    ## print("Common keys: {}".format(common_key_list))

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = []
    for key in common_key_list:
        mesh_point = results_dict[key]
        value = mesh_point.get_energy(qn)
        if (prune and np.isnan(value)):
            continue
        table_data += [
            key + (value,)
        ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = list(key_descriptor)+[("value",float)]
     )
    return table


def make_energy_difference_table(mesh_data_pair,key_descriptor,qn_pair,key_list=None,prune=False):
    """Generate energy tabulation.

    The key descriptor is used to provide the parameter columns of the
    tabulation (and may in general be different from the key which was
    used for sorting the data, e.g., it might add additional parameter
    columns).

    It is expected that the mesh data have already been consolidated to make the
    key values unique.  Only key tuples which are present in both sets of mesh
    data are retained.

    Results are then sorted by key tuple.

    Differences are taken as (E1-E2), i.e., the second data set is the reference
    data set.

    Data format:
        param1 param2 ... E

    Arguments:
        mesh_data_pair (list of ResultsData): pair of data sets
        key_descriptor (tuple of tuple): dtype descriptor for key
        qn_pair (tuple): pair of quantum numbers (J,g,n) of levels of interest
        key_list (list of tuple, optional): key list for data points; defaults to
            common keys from the two data sets
        prune (bool, optional): whether or not to prune nan values from table

    Returns:
       (array): data table

    """

    # unpack arguments
    (mesh_data1,mesh_data2) = mesh_data_pair
    (qn1,qn2) = qn_pair

    # process results into dictionaries
    results_dict1 = make_results_dict(mesh_data1,key_descriptor)
    results_dict2 = make_results_dict(mesh_data2,key_descriptor)

    # find common keys
    if (key_list is not None):
        common_key_list = key_list
    else:
        common_key_list = common_keys([results_dict1,results_dict2])
    ## print("Common keys: {}".format(common_key_list))

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = []
    for key in common_key_list:
        mesh_point1 = results_dict1[key]
        mesh_point2 = results_dict2[key]
        value1 = mesh_point1.get_energy(qn1)
        value2 = mesh_point1.get_energy(qn2)
        value = value1-value2
        if (prune and np.isnan(value)):
            continue
        table_data += [
            key + (value,)
        ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = list(key_descriptor)+[("value",float)]
     )
    return table

def make_radius_table(mesh_data,key_descriptor,radius_type,qn):
    """ Generate radius tabulation.

    See docstring for make_energy_table for sorting and tabulation conventions.

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        radus_type (str): radius type code (rp/rn/r)
        qn (tuple): quantum numbers (J,g,n) of level to retrieve

    Returns:
       (array): data table
    """

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = [
        key_function(mesh_point) + (mesh_point.get_radius(radius_type,qn),)
        for mesh_point in mesh_data
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = list(key_descriptor)+[("r",float)]
     )
    return table

def make_rme_table(mesh_data,key_descriptor,observable,qnf,qni):
    """ Generate reduced matrix element (RME) tabulation.

    Data format:
        <key> RME
    where actual label columns depend up on key_descriptor, e.g.:
        Nmax hw RME
        Nsigmamax Nmax hw RME

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        observable_name (str): key naming observable
        qnf, qni (tuple): quantum numbers (J,g,n) of final and initial levels

    Returns:
       (array): data table
    """

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = [
        key_function(mesh_point) + (mesh_point.get_rme(observable,(qnf,qni)),)
        for mesh_point in mesh_data
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = list(key_descriptor)+[("RME",float)]
     )
    return table

def make_rtp_table(mesh_data,key_descriptor,observable,qnf,qni):
    """ Generate reduced transition probability (RTP) tabulation.

    Data format:
        <key> RTP
    where actual label columns depend up on key_descriptor, e.g.:
        Nmax hw RTP
        Nsigmamax Nmax hw RTP

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        observable_name (str): key naming observable
        qnf, qni (tuple): quantum numbers (J,g,n) of final and initial levels

    Returns:
       (array): data table
    """


    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = [
        key_function(mesh_point) + (mesh_point.get_rtp(observable,(qnf,qni)),)
        for mesh_point in mesh_data
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = list(key_descriptor)+[("RTP",float)]
     )
    return table

################################################################
# tabulation functions -- listing by level
################################################################

def make_level_table(mesh_point,levels=None,energy_cutoff=None):
    """Generate listing of level energies from single run.

    Data format:
      J gex n E

    Legacy format (for comparison):
      seq J gex n T E

    Arguments:
        mesh_point (BaseResultsData): data for mesh point
        levels (list of tuples): levels to include or None for all (default: None)
        energy_cutoff (float,optional): energy cutoff to limit output size

    Returns:
       (array): data table

    >>> level_table = mfdnres.analysis.make_level_table(mesh_point)
    >>> mfdnres.tools.write_table(
    >>>     filename,"{:4.1f} {:1d} {:3d} {:7.3f}",
    >>>     level_table
    >>> )


    """

    # determine list of levels to use
    if (levels is None):
        levels = mesh_point.levels

    # tabulate values
    table_data = [
        qn + (mesh_point.get_energy(qn),)
        for qn in levels
        if not ((energy_cutoff is not None) and (mesh_point.get_energy(qn)>energy_cutoff))  # clumsy!
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = [("J",float),("gex",int),("n",int),("E",float)]
     )
    return table

################################################################
# calculate effective angular momentum
################################################################

def effective_am(J_sqr):
    """  Convert mean square angular momentum to effective angular momentum.

    Args:
        J_sqr (float): value representing <J.J>

    Returns:
        (float): effective J
    """

    J = (math.sqrt(4*J_sqr+1)-1)/2
    return J

if (__name__ == "__main__"):
    pass
