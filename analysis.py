""" analysis.py -- analysis tools for processing MFDn results file

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    6/2/15 (mac): Initiated (as mfdn_analysis.py).
    6/5/15 (mac): Restructure as subpackage.
    7/7/17 (mac): Remove obsolete import_res_files.
    7/8/17 (mac): Start adding new tabulation functions which return
      arrays of data rather than writing to file.
    7/11/17 (mac): Pull out extrapolation function to remove scipy dependence.
    7/30/17 (mac): Pull out band analysis functions to band.py.
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
            print("  filename {} key {}".format(mesh_point.filename,key))

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

    key_value_pairs = selection_dict.items()  # convert to view as iterable of pairs for cleaner debugging
    if (verbose):
        print("selection key-value pairs: {}".format(tuple(key_value_pairs)))
    key_function = make_key_function(key_value_pairs)
    selected_values = tuple([value for (_,value) in key_value_pairs])

    new_mesh_data = [
        mesh_point
        for mesh_point in mesh_data
        if (key_function(mesh_point) == selected_values)
    ]
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

def make_energy_table(mesh_data,key_descriptor,qn):

    """Generate energy tabulation.

    The key descriptor is used to provide the parameter columns of the
    tabulation (and may in general be different from the key which was
    used for sorting the data, e.g., it might add additional parameter
    columns).

    It is expected that the mesh data have already been sorted and
    consolidated to make the key values unique, or else the table will
    contain duplicate entries.

    Data format:
        param1 param2 ... E

    Arguments:
        mesh_data (list of ResultsData): data for mesh points
        key_descriptor (tuple of tuple): dtype descriptor for key
        qn (tuple): quantum numbers (J,g,n) of level to retrieve

    Returns:
       (array): data table

    """

    # tabulate values
    key_function = make_key_function(key_descriptor)
    table_data = [
        key_function(mesh_point) + (mesh_point.get_energy(qn),)
        for mesh_point in mesh_data
    ]

    # convert to structured array
    table = np.array(
        table_data,
        dtype = list(key_descriptor)+[("E",float)]
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
    """ Generate RME tabulation.

    Data format:
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
    """ Generate RTP tabulation.

    Data format:
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
        levels = mesh_point.get_levels()

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


################################################################
################################################################
# LEGACY output table writing functions
################################################################
################################################################
#
# The following LEGACY functions actual *write* tables.  Instead, new
# functions should just *return* tables.
#
# They also include MFDn-specific assumptions, such as the presence of
# the isospin label, and use the deprecated get_property accessor.

################################################################
# LEGACY -- output tabulation: basic levels
################################################################

def write_level_table(results,filename,levels=None):
    """Writes basic level data.

        "seq"=0 J gex n T Eabs

    The field seq is included for historical reasons but written as a
    dummy zero.

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        levels (list of tuples): levels to include or None for all (default: None)

    """

    # determine list of levels to use
    if (levels is None):
        levels = results.get_levels()

    # assemble lines
    lines = []
    for qn in levels:
        line = (
            "{:1d} {:6.3f} {:1d} {:2d} {:6.3f} {:8.3f}"
            "\n"
        ).format(
            0,
            results.get_property(qn,"J"),
            results.get_property(qn,"g"),
            results.get_property(qn,"n"),
            results.get_property(qn,"T"),
            results.get_energy(qn)
        )
        lines.append(line)

    # write to file
    with open(filename,"wt") as fout:
        fout.writelines(lines)

################################################################
# LEGACY -- output tabulation: radii
################################################################

def write_radii_table(results,filename,levels=None):
    """ LEGACY -- Writes radius data.

        "seq"=0 J gex n T rp rn r

    The field seq is included for historical reasons but written as a
    dummy zero.

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        levels (list of tuples): levels to include or None for all (default: None)

    """

    # determine list of levels to use
    if (levels is None):
        levels = results.get_levels()

    # assemble lines
    lines = []
    for qn in levels:
        line = (
            "{:1d} {:6.3f} {:1d} {:2d} {:6.3f} {:8.3f} {:8.3f} {:8.3f}"
            "\n"
        ).format(
            0,
            results.get_property(qn,"J"),
            results.get_property(qn,"g"),
            results.get_property(qn,"n"),
            results.get_property(qn,"T"),
            results.get_tbo(qn,"rp"),
            results.get_tbo(qn,"rn"),
            results.get_tbo(qn,"r")
        )
        lines.append(line)

    # write to file
    with open(filename,"wt") as fout:
        fout.writelines(lines)

################################################################
# LEGACY -- output tabulation: angular momentum
################################################################

def write_level_am_table(results,filename,levels=None,default=np.nan):
    """LEGACY  -- Writes basic level data.

        seq J gex n T Eabs

    The field seq is included for historical reasons but written as a
    dummy zero.

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        levels (list of tuples): levels to include or None for all (default: None)
        default (float): value to use for missing numerical entries (default: np.nan)

    """

    # determine list of levels to use
    if (levels is None):
        levels = results.get_levels()

    # assemble lines
    lines = []
    for qn in levels:
        # basic properties
        line = ("{:1d} {:6.3f} {:1d} {:2d} {:6.3f} {:8.3f}").format(
            0,
            results.get_property(qn,"J"),
            results.get_property(qn,"g"),
            results.get_property(qn,"n"),
            results.get_property(qn,"T"),
            results.get_energy(qn)
        )

        # angular momenta
        tbo = results.tbo[qn]
        value_sqr_list = [tbo.get("L2",default),tbo.get("Sp2",default),tbo.get("Sn2",default),tbo.get("S2",default)]
        value_list = list(map(effective_am,value_sqr_list))
        value_format = " {:6.3f}"
        line += (4*value_format).format(*value_list)

        # finalize line
        line += "\n"
        lines.append(line)

    # write to file
    with open(filename,"wt") as fout:
        fout.writelines(lines)

if (__name__ == "__main__"):
    pass
