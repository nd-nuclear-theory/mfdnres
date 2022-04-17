"""analysis.py -- analysis tools for processing MFDn results file

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
    10/27/18 (mac): Add ability to transform keys in results dictionary.
    02/25/19 (pjf):
        - Add sorting and pruning of nan values to radius tabulation.
        - Add moment tabulation with sorting and pruning of nan values.
        - Add __MatchAny to act as wildcard for matching.
    02/27/19 (mac): Extend make_energy_difference_table to handle opposite parities via key
        transformation on reference state.
    03/31/19 (mac): Add format string for level table.
    04/02/19 (mac): Update am handling (add make_am_table and move effective_am to tools).
    06/02/19 (mac): Update am handling (add make_am_table and move effective_am to tools).
    06/24/19 (mac):
        - Add merged_mesh and ancillary dictionary utilities for manipulating
            params (refactored from tabulate_band_11be_11be-xfer.py).
        - Revise and extend make_level_table_* family of functions.
    07/11/19 (mac): Add mesh_key_listing for mesh diagnostic listing.
    06/18/20 (pjf): Make sorted_mesh_data only sort, not consolidate.
    03/16/21 (mac): Standardize structured array observable column name as "value".
    04/02/21 (mac): Add make_obs_table for generic single-observable tabulation.
    10/12/21 (pjf):
        - Add reverse option to sorted_mesh_data.
        - Add nonspurious_to_spurious_qn.
        - Use dict comprehension in subdict to preserve order.
        - Add preprocessor option to merged_mesh; keep common params in merged mesh point.
    04/17/22 (mac): Make make_obs_table provide np.nan value when extractor raises exception.
"""

import copy
import functools
import itertools
import more_itertools
import math

import numpy as np

################################################################
# global constants
################################################################

# electromagnetic g factors (glp,gln,gsp,gsn)
MU_EM = np.array([1,0,5.586,-3.826])

# singleton representing match any
class __MatchAny(tuple):
    def __str__(self):
        return "ANY"
    def __repr__(self):
        return "ANY"
    def __eq__(self, other):
        return True
    def __ne__(self, other):
        return False
    def __hash__(self):
        return id(self)
ANY = __MatchAny()

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
# key-based operations on mesh
################################################################

def extract_key(key_fields,results_data):
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
    key = tuple([results_data.params[key] for key in key_fields])
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
        key_transformation=None,
        verbose=False
):
    """Load mesh data into dictionary, using specified parameter tuple as
       key.

    Example key descriptor:

       (("Nsigmamax",int),("Nmax",int),("hw",float))

    Example:
        >>> KEY_DESCRIPTOR_NMAX_HW = (("Nmax",int),("hw",float))

    For now, in the event that the same mesh point arises multiple
    times on input (i.e., a given value for the key tuple is
    duplicated), the final occurrence overwrites any earlier
    occurrences.  In the future, a more sophisticated "merging"
    process might be appropriate.

    An optional key transformation is useful for, e.g., shifting the Nmax value
    stored in the dictionary when results are to be used as reference results
    for the space of opposite parity.

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        key_transformation (callable,optional): transformation function to apply to key tuple
        verbose (bool,optional): verbose output

    Returns:
        (dict): mapping from key tuple to data object

    """

    key_function = make_key_function(key_descriptor)

    results_dict = dict()
    for results_data in mesh_data:

        # make key
        key = key_function(results_data)
        if (key_transformation is not None):
            key = key_transformation(key)
        if (verbose):
            print("  make_results_dict: filename {} key {}".format(results_data.filename,key))

        # store data point
        if (key not in results_dict):
            # save mesh point
            results_dict[key] = results_data
        else:
            # TODO: do smart merge "update" on existing mesh point
            # overwrite mesh point
            results_dict[key] = results_data

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
        results_data
        for results_data in mesh_data
        if (key_function(results_data) == selected_values)
    ]
    if (verbose):
        print("  selected_mesh_data: selected mesh points {}".format(len(new_mesh_data)))
    return new_mesh_data

def sorted_mesh_data(mesh_data, key_descriptor, reverse=False, verbose=False):
    """Sort mesh data, using specified parameter tuple as key.

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
        (list): sorted list of data objects
    """

    key_function = make_key_function(key_descriptor)
    new_mesh_data = sorted(mesh_data, key=key_function, reverse=reverse)

    return new_mesh_data

################################################################
# quantum numbers within nonspurious space
################################################################

def nonspurious_to_spurious_qn(mesh_point, qn, threshold=0.9):
    """Converts qn in non-spurious space to qn in spurious space

    Arguments:
        mesh_point (MFDnResultsData): mesh point
        qn (tuple): desired quantum numbers in nonspurious space
        threshold (float, optional): threshold value of <Ncm> for
            spuriosity

    Returns:
        (tuple): qn in the full space, None if not found
    """
    (J,g,target_n) = qn
    n = 0
    for level in mesh_point.levels:
        if (J,g) == level[0:2] and mesh_point.mfdn_tb_expectations["Ncm"][level] < threshold:
            n += 1
        if n == target_n:
            return level
    return None


################################################################
# consolidation of results data objects by shared parameters
################################################################

def subdict(d,keys):
    """Subset a dictionary to only include given keys.

    Arguments:
        d (dict): dictionary to subset
        keys (list): list of key values

    Returns:
        subdict (dict): dictionary with reduced set of keys

    """
    subitems = {
        key: value
        for key,value in d.items()
        if key in keys
    }
    return subitems

def dict_items(d):
    """ Convert dictionary to sorted tuple of items.

    Arguments:
        d (dict): dictionary

    Returns:
        dict_items (tuple of tuple): tuple of items
    """
    return tuple(sorted(d.items()))

## def make_subdict_function(keys):
##     """ Generate a function to subset a dictionary to the given keys.
##
##     Arguments:
##         keys (list): list of key values
##
##     Returns:
##         subdict_function (callable): function to extract subsetted dictionary
##     """
##     return functools.partial(subdict,keys=keys)

def make_params_subdict_items_function(keys):
    """ Generate a function to subset a dictionary to the given keys.

    Arguments:
        keys (list): list of keys

    Returns:
        params_subdict_function (callable): function to extract subsetted dictionary from ResultsData params
    """

    return lambda results : dict_items(subdict(results.params,keys))

def mesh_key_listing(mesh,keys,verbose=False):
    """Provide sorted list of key tuples appearing in mesh.

    In verbose operation, can discared result and simply use to print
    diagnostic information on mesh.

    Arguments:

        mesh (list of ResultsData): mesh to be merged (entrywise)

        keys (list): list of keys

    Return:

        mesh (list of tuples): listing of key-value tuples

    """

    keyfunc = make_params_subdict_items_function(keys)
    mesh_keys = map(keyfunc,mesh)
    mesh_keys = sorted(list(mesh_keys))

    if (verbose):
        for key_tuple in mesh_keys:
            print(key_tuple)

    return mesh_keys

def merged_mesh(mesh,keys,preprocessor=None,postprocessor=None,verbose=False):
    """Obtain results mesh in which results data objects sharing same parameter
    values are merged.

    Params: The "params" for the resulting object are the minimal common set of
    parameters used in the sorting key which defines the merge (see "Merge key"
    below).

    Merge key: The key function should return a tuple of key-value pairs, used
    to distinguish equivalent (for merger) vs. nonequivalent mesh points.  The
    params dictionary for each new ResultsData is constructed from the key-value
    items in the tuble returned by the key function.  Therefore, the key-value
    items should include not only the parameters required in the merger, to
    distinguish equivalent vs. nonequivalent mesh points, but also any ancillary
    key-value pairs which might be required in later processing of the mesh
    point, even if it is known a priori that all mesh elements share the same
    value for this parameter making it irrelevant to the merging process (e.g.,
    the "nuclide" parameter should be included if that parameter is later used
    in filename construction, or the "parity" parameter may be used in later
    mesh selection, even when the "Nmax" parameter would be sufficient to
    distinquish nonequivalent mesh points).

    Preservation of original objects: Each ResultsData object in the returned
    mesh is obtained by first defining a new empty object.  Then all ResultsData
    objects to be merged are copied into this empty object by the update method.
    This ensures that no original ResultsData object is modified in the
    process.

    Example:

    >>> # merge results data (of different M) by shared mesh parameters
    >>> mesh_data = mfdnres.analysis.merged_mesh(
    >>>     mesh_data,
    >>>     ("nuclide","interaction","coulomb","hw","Nmax","parity"),
    >>>     verbose=True
    >>> )

    Arguments:

        mesh (list of ResultsData): mesh to be merged (entrywise)

        keys (list): list of keys

        preprocessor (callable): callable to apply to all mesh points in initial
           un-merged mesh (e.g., to define extra parameters)

        postprocessor (callable): callable to apply to all mesh points in final
           merged mesh (e.g., to define extra parameters)

    Return:

        merged_mesh (list of ResultsData): merged mesh

    """

    try:
        results_type = type(mesh[0])
    except IndexError:  # empty mesh has no result type
        return []

    # preprocess mesh
    source_mesh = copy.deepcopy(mesh)
    if (preprocessor is not None):
        for results_data in source_mesh:
            preprocessor(results_data)

    # presort mesh (required for groupby)
    keyfunc = make_params_subdict_items_function(keys)
    sorted_mesh = sorted(source_mesh,key=keyfunc)

    # group and merge mesh points
    target_mesh = []
    for group_key, group in itertools.groupby(sorted_mesh,keyfunc):

        if (verbose):
            print()
            print("merging key: {}".format(group_key))

        # initialize new mesh point
        results = results_type()
        results.params=dict(group_key)

        # save params for later intersection
        params_list = []

        # merge in data
        for other_results in group:
            if (verbose):
                print("  contributing params: {}".format(other_results.params))
                print("  contributing levels: {}".format(other_results.levels))
            results.update(other_results)
            params_list += [other_results.params]

        # use merged params
        shared_param_keys = set.intersection(*(set(params.keys()) for params in params_list))
        for key in shared_param_keys:
            if more_itertools.all_equal(params[key] for params in params_list):
                results.params[key] = params_list[0][key]

        # save new mesh point
        if (verbose):
            print("  merged params: {}".format(results.params))
            print("  merged levels: {}".format(results.levels))

        target_mesh.append(results)

    # postprocess mesh
    if (postprocessor is not None):
        for results_data in target_mesh:
            postprocessor(results_data)

    return target_mesh

################################################################
# tabulation functions -- observable vs. parameters
################################################################

# "make table" functions now DEPRECATED in favor of pandas-based analysis in
# mfdnres.data.

def make_obs_table(mesh_data,key_descriptor,obs_extractor,key_list=None,prune=False):
    """Generate tabulation of generic observable on parameter mesh.

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
        obs_extractor (callable): function to extract observable value from results data object
        key_list (list of tuple, optional): key list for data points; defaults to
            keys for which qn found
        prune (bool, optional): whether or not to prune nan values from table

    Returns:
       (np.ndarray): data table as structured array

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
        results_data = results_dict[key]
        try:
            value = obs_extractor(results_data)
        except:
            value = np.nan
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

def make_energy_table(mesh_data,key_descriptor,qn,key_list=None,prune=False):
    """Generate energy tabulation on parameter mesh.

    NOTE: May now be trivially reimplemented in terms of make_obs_table.

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
       (np.ndarray): data table as structured array

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
        results_data = results_dict[key]
        value = results_data.get_energy(qn)
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

def make_energy_difference_table(
        mesh_data_pair,key_descriptor,qn_pair,
        reference_key_transformation=None,key_list=None,prune=False,
        verbose=False
):
    """Generate energy tabulation on parameter mesh.

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

    A reference_key_transformation is necessary if the excited state and
    reference state are expected to live on different mesh points, e.g., for
    states of opposite parity.

    Data format:
        param1 param2 ... E

    Arguments:
        mesh_data_pair (list of ResultsData): pair of data sets
        key_descriptor (tuple of tuple): dtype descriptor for key
        qn_pair (tuple): pair of quantum numbers (J,g,n) of levels of interest
        reference_key_transformation (callable,optional): transformation function to apply
            to key tuple for reference data set
        key_list (list of tuple, optional): key list for data points; defaults to
            common keys from the two data sets
        prune (bool, optional): whether or not to prune nan values from table

    Returns:
       (np.ndarray): data table as structured array

    """

    # unpack arguments
    (mesh_data1,mesh_data2) = mesh_data_pair
    (qn1,qn2) = qn_pair

    # process results into dictionaries
    results_dict1 = make_results_dict(mesh_data1,key_descriptor)
    results_dict2 = make_results_dict(mesh_data2,key_descriptor,key_transformation=reference_key_transformation)

    # find common keys
    if (key_list is not None):
        common_key_list = key_list
    else:
        common_key_list = common_keys([results_dict1,results_dict2])
    if (verbose):
        print("results_dict2 keys: {}".format(sorted(results_dict2.keys())))
        print("Common keys: {}".format(common_key_list))

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = []
    for key in common_key_list:
        results_data1 = results_dict1[key]
        results_data2 = results_dict2[key]
        value1 = results_data1.get_energy(qn1)
        value2 = results_data2.get_energy(qn2)
        value = value1-value2
        if (prune and np.isnan(value)):
            continue
        table_data += [
            key + (value,)
        ]
    if (verbose):
        print("Table data: {}".format(table_data))

    # convert to structured array
    table = np.array(
        table_data,
        dtype = list(key_descriptor)+[("value",float)]
     )
    return table

def make_radius_table(mesh_data,key_descriptor,radius_type,qn,key_list=None,prune=False):
    """ Generate radius tabulation on parameter mesh.

    NOTE: May now be trivially reimplemented in terms of make_obs_table.

    See docstring for make_energy_table for sorting and tabulation conventions.

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        radus_type (str): radius type code ("rp","rn","r")
        qn (tuple): quantum numbers (J,g,n) of level to retrieve

    Returns:
       (np.ndarray): data table as structured array
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
        results_data = results_dict[key]
        value = results_data.get_radius(radius_type,qn)
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

def make_am_table(mesh_data,key_descriptor,qn):
    """ Generate effective angular momentum tabulation on parameter mesh.

    TODO: Update to sort keys with common_key_list API.

    Data format:
        <key> L Sp Sn S
    where actual label columns depend up on key_descriptor, e.g.:
        Nmax hw ...
        Nsigmamax Nmax hw ...

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        qn (tuple): quantum numbers (J,g,n) of level to retrieve

    Returns:
       (np.ndarray): data table as structured array
    """

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = [
        key_function(results_data) + (
            results_data.get_am("L",qn),
            results_data.get_am("Sp",qn),
            results_data.get_am("Sn",qn),
            results_data.get_am("S",qn),
        )
        for results_data in mesh_data
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = list(key_descriptor)+[("L",float),("Sp",float),("Sn",float),("S",float)]
     )
    return table

def make_rme_table(mesh_data,key_descriptor,observable,qnf,qni):
    """ Generate reduced matrix element (RME) tabulation on parameter mesh.

    NOTE: May now be trivially reimplemented in terms of make_obs_table.

    TODO: Update to sort keys with common_key_list API.

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
       (np.ndarray): data table as structured array
    """

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = [
        key_function(results_data) + (results_data.get_rme(observable,(qnf,qni)),)
        for results_data in mesh_data
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = list(key_descriptor)+[("value",float)]
     )
    return table

def make_moment_table(mesh_data,key_descriptor,observable,qn,key_list=None,prune=False):
    """ Generate moment tabulation on parameter mesh.

    NOTE: May now be trivially reimplemented in terms of make_obs_table.

    Data format:
        <key> moment
    where actual label columns depend up on key_descriptor, e.g.:
        Nmax hw moment
        Nsigmamax Nmax hw moment

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        observable_name (str): key naming observable
        qn (tuple): quantum numbers (J,g,n) of level

    Returns:
       (np.ndarray): data table as structured array
    """
    # process results into dictionaries
    results_dict = make_results_dict(mesh_data,key_descriptor)

    # find common keys
    if (key_list is not None):
        common_key_list = key_list
    else:
        common_key_list = common_keys([results_dict])

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = []
    for key in common_key_list:
        results_data = results_dict[key]
        value = results_data.get_moment(observable,qn)
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

def make_rtp_table(mesh_data,key_descriptor,observable,qnf,qni,key_list=None,prune=False,verbose=False):
    """ Generate reduced transition probability (RTP) tabulation on parameter mesh.

    NOTE: May now be trivially reimplemented in terms of make_obs_table.

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
       (np.ndarray): data table as structured array
    """

    # process results into dictionaries
    results_dict = make_results_dict(mesh_data,key_descriptor)

    # find common keys
    if (key_list is not None):
        common_key_list = key_list
    else:
        common_key_list = common_keys([results_dict])

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = []
    for key in common_key_list:
        results_data = results_dict[key]
        value = results_data.get_rtp(observable,(qnf,qni))
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

################################################################
# tabulation functions -- listing by level
################################################################

FORMAT_STRING_LEVEL_TABLE = "{:4.1f} {:1d} {:3d} {:7.3f}"
def make_level_table_energy(results_data,levels=None):
    """Make tabulation of energies by level, from single calculation mesh point.

    Data format (when written to file):
      J g n E

    Legacy format (for pre-2018 files, for comparison):
      seq J gex n T E

    Example (writing out level table):

    >>> level_table = mfdnres.analysis.make_level_table(results_data)
    >>> mfdnres.tools.write_table(
    >>>     filename,
    >>>     mfdnres.analysis.FORMAT_STRING_LEVEL_TABLE,
    >>>     level_table
    >>> )

    Limitation: All levels must be from the same results data object.  This
    presupposes that calculations from different M have already been merged.

    Arguments:
        results_data (BaseResultsData): data for mesh point
        levels (list of tuples): levels to include or None for all (default: None)

    Returns:
       (array): table of (J,g,n,E)

    """

    # determine list of levels to use
    if (levels is None):
        levels = results_data.levels

    # tabulate values
    qn_dtype = [("J",float),("gex",int),("n",int)]
    values_dtype = [("E",float)]
    def values(qn):
        return (results_data.get_energy(qn),)
    table_data = [
        qn + values(qn)
        for qn in levels
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = qn_dtype+values_dtype
    )
    return table

FORMAT_STRING_AM_TABLE = "{:4.1f} {:1d} {:3d} {:7.3f} {:7.3f} {:7.3f} {:7.3f}"
def make_level_table_am(results_data,levels=None):
    """Make tabulation of angular momenta by level, from single calculation mesh point.

    Arguments:
        results_data (BaseResultsData): data for mesh point
        levels (list of tuples): levels to include or None for all (default: None)

    Returns:
       (array): table of (J,g,n,L,Sp,Sn,S)

    """

    # determine list of levels to use
    if (levels is None):
        levels = results_data.levels

    # tabulate values
    qn_dtype = [("J",float),("gex",int),("n",int)]
    values_dtype = [("L",float),("Sp",float),("Sn",float),("S",float)]
    def values(qn):
        return (
            results_data.get_am("L",qn),
            results_data.get_am("Sp",qn),
            results_data.get_am("Sn",qn),
            results_data.get_am("S",qn),
        )
    table_data = [
        qn + values(qn)
        for qn in levels
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = qn_dtype+values_dtype
     )
    return table

def format_string_decomposition_table(decomposition_length):
    """Generate variable-length format string for decomposition table."""
    return "{:4.1f} {:1d} {:3d}"+decomposition_length*" {:7.5f}"
def make_level_table_decomposition(results_data,decomposition_length,levels=None):
    """Make tabulation of Nex decomposition by level, from single calculation mesh point.

    Decomposition length is typically

        decomposition_length = (Nmax - (Nmax%2))//2 + 1

    Arguments:
        results_data (BaseResultsData): data for mesh point
        levels (list of tuples): levels to include or None for all (default: None)
        decomposition_length (int): length to which to pad or truncate all
            decompositions, required since used in dtype construction

    Returns:
       (array): table of (J,g,n,P0,P1,...)

    """

    # determine list of levels to use
    if (levels is None):
        levels = results_data.levels

    # tabulate values
    qn_dtype = [("J",float),("gex",int),("n",int)]
    values_dtype = decomposition_length*[("",float),]
    def values(qn):
        decomposition=results_data.get_decomposition("Nex",qn)
        if (decomposition is None):
            value_list = [np.nan]*decomposition_length
        else:
            value_list = decomposition[:decomposition_length]
        return tuple(value_list)
    table_data = [
        qn + values(qn)
        for qn in levels
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = qn_dtype+values_dtype
     )
    return table


if (__name__ == "__main__"):
    pass
