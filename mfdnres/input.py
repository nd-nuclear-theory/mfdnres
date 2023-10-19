"""input.py

    Import control code for results files.

    Caveat: "input" is also Python built-in function name.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    05/31/15 (mac): Initiated (as mfdn_res.py).
    06/05/15 (mac): Allow user-supplied res file parser.
    06/05/15 (mac): Restructure as subpackage.
    06/29/17 (jbutler): Added in inheritance for SpNCCI, updated documentation
    07/07/17 (mac):
        - Generalize slurp_res_files to take list of directory names (after old
        analysis.import_res_files.
        - Add directory name generation utility res_file_directory.
    07/09/17 (mac):
        - Restore read_file to be simple dispatch function.
        - Extract SpNCCIMeshPointData.
    10/06/17 (mac): Extract MFDnRunData subclass to mfdn.py.
    10/10/17 (mac): Extract results data base class to results_data.py.
    09/20/18 (pjf): Store filename in ResultsData.
    02/22/19 (pjf): Store only basename of filename in ResultsData.
    06/24/19 (mac): Update res_file_directory construction for new default
        location "results/res".
    06/26/19 (mac): Rename from res.py to input.py, and incorporate filename
        parsing control from descriptor.py.
    05/05/20 (mac): Suppress duplicate input directories in slurp_res_files.
    06/17/20 (pjf): Add code registration and detection from filename.
    09/02/20 (pjf): Add autodetection of filename format.
    09/07/20 (pjf): Fix filename parsing.
    10/12/21 (pjf): Print filename info if parser throws error.
    06/27/23 (mac): Add basic mesh data caching facility, based on code
        from pjf lenpic-analysis-2022.
    07/08/23 (mac): Add run_stem option for res_file_directory().
"""

import glob
import os
import pickle

import numpy as np

################################################################
# filename utility
################################################################

def res_file_directory(
        username, code, run_number, *,
        run_stem="run", results_dir="results", res_file_subdir=os.path.join("results","res"),
):
    """Construct full path to res file directory, given user, code, and run.

        This function assumes directory naming conventions appropriate
        to mcscript archive files.

        Arguments:

            username (str): user name (e.g., "mcaprio")

            code (str): code name (e.g., "spncci")

            run_number (str): run name "tail" (e.g., "mac0424")

            run_stem (str, optional): run name "stem" (normally "run")

            results_dir (str,optional): name of top-level results directory within GROUP_HOME

            res_file_subdir (str,optional): name of subdirectory within results directory (can be None for flat structure)

        Environment:
            GROUP_HOME: directory name for group top-level results directory
              (e.g., "/afs/crc.nd.edu/group/nuclthy" for shared group results directory,
               or, for local work in your home directory, you may set equal to HOME)

        >>> res_file_directory("mcaprio","mfdn","mac0563")

            /afs/crc.nd.edu/group/nuclthy/results/mcaprio/mfdn/runmac0563/results/res

        >>> res_file_directory("amccoy","spncci","aem0097",res_file_subdir="results")

            /afs/crc.nd.edu/group/nuclthy/results/amccoy/spncci/runaem0097/results

    """

    group_home = os.environ.get("GROUP_HOME")
    if (type(group_home) is not str):
        raise(ValueError("Need to set environment variable GROUP_HOME"))

    res_directory = os.path.join(group_home, results_dir, username, code, run_stem+run_number)
    if (res_file_subdir is not None):
        res_directory = os.path.join(res_directory,res_file_subdir)

    return res_directory

################################################################
# filename parser registry
################################################################

# global registration variables
filename_format_parser = {}

def register_filename_format(format_name,parser):
    """Register information for parsing filename.

    Args:
        format_name (str): name for filename format
        parser (callable): function for parsing filename

    """
    if format_name == "ALL":
        raise ValueError("filename format code ALL is reserved")

    filename_format_parser[format_name] = parser

################################################################
# filename parser control
################################################################

def parse_filename(filename, filename_format="ALL"):
    """Parse results filename.

    Only the basename is considered, extracted via os.path.basename,
    i.e., any preceding path is ignored.

    A filename parsing function is assumed to provide the following
    mandatory fields:

        "run" (str): run name (may be null)
        "descriptor" (str): the part of the filename which
             describes the run parameters
        "Z", "N" (int): proton and neutron numbers

    E.g., under "format5ho", the filename

        "run0364-mfdn-Z4-N3-JISP16-1-hw20.000-aL100-Nmax10-MM1-lan1000.res"

    yields

        "run" : "0364"
        "descriptor" : "Z4-N3-JISP16-1-hw20.000-aL100-Nmax10-MM1-lan1000"
        "Z" : 4
        "N" : 3
        "interaction" : "JISP16"
        ...

    The wrapper parse_filename will add the field "nuclide" as a
    tuple of int, e.g.,

        "nuclide" : (4,3)

    Args:
        filename (str): filename to parse
        filename_format (str, optional): filename format to match

    Returns: (dict) : dictionary with keys for parameters ("run",
        "descriptor", "Z", "N", ...) parsed from filename, plus
        "nuclide" as a tuple of int

    """

    # parse filename
    basename = os.path.basename(filename)

    # disable parsing if filename_format is None
    if filename_format is None:
        return {"filename": filename}

    # try all filename formats for special value ALL
    if filename_format == "ALL":
        for parser in filename_format_parser.values():
            try:
                info = parser(basename)
            except ValueError:
                info = {}
                continue
            else:
                break
    elif filename_format in filename_format_parser:
        parser = filename_format_parser[filename_format]
        info = parser(basename)
    else:
        raise KeyError("unknown filename_format={}".format(filename_format))


    # define nuclide tuple
    info["filename"] = filename
    if ("Z" in info) and ("N" in info):
        info["nuclide"] = (info["Z"],info["N"])

    return info

#################################################
# data parser registry
#################################################

# global registration variables
data_format_parser = {}

def register_data_format(format_name,parser):
    """Register information for parsing res file.

    Args:
        format_name (str): name for res file format
        parser (callable): function for parsing file stream

    """

    data_format_parser[format_name] = parser

################################################################
# code name registry
################################################################

# global registration variables
code_name_map = {None: None}

def register_code_name(code_name,format_name):
    """Register information for deducing res format.

    Args:
        code_name (str): name for code
        format_name (str): name for filename format

    """
    if format_name not in data_format_parser:
        raise ValueError("unknown format_name: {:s}".format(format_name))
    code_name_map[code_name] = format_name

##################################################
# data file import control
##################################################

def read_file(filename, *, res_format=None, filename_format=None, params=None, verbose=False):
    """Extract results from single results file.

    Dispatches filename to appropriate filename parser.  Dispatches
    file contents to appropriate results file parser.  Parameter
    values obtained from the file name are merged into the parameter
    dictionaries stored with each mesh point.

    The results will be a list of results data objects, one for each
    mesh point within the results file.  (Most commonly, the results
    file contains contains the results for just a single mesh point,
    so this will be a list containing just one object, but, e.g.,
    spncci can calculate multiple hw mesh points in a single run.)
    The results data objects will be children of the interface class
    BaseResultsData.

    Arguments:
        filename (str): filename for results file

        res_format (str, optional): identifier string for the results file
            parser to use

        filename_format (str, optional): identifier string for the results
            filename parser to use

        params (dict, optional): supplementary parameters to append to params attribute

        verbose (bool, optional): enable debugging output

    Returns:
        (list of ResultsData): list of mesh point data objects

    """

    # parse results filename for any supplementary run parameters
    info_from_filename = parse_filename(filename, filename_format)

    if res_format is None:
        if info_from_filename.get("code_name") is not None:
            res_format = code_name_map[info_from_filename["code_name"]]
        else:
            raise ValueError("unable to deduce res_format")

    # parse results file contents for run parameters and data
    if (verbose):
        print("  read_file: filename {}".format(filename))
    with open(filename,'rt') as fin:
        try:
            results_list = data_format_parser[res_format](fin, verbose=verbose)
        except Exception as e:
            print("filename {} filename_format {} res_format {}".format(filename, filename_format, res_format))
            raise e
    if (verbose):
        print("  read_file: mesh points {:d}".format(len(results_list)))

    # augment parameters with those obtained from filename
    #
    # Note: The parameter values obtained from the filename will
    # *override* any parameter values obtained by parsing the results
    # file.  So beware that parameter values formatted for the
    # filename might have lower precision than those stored in the
    # results file.

    for results in results_list:
        results.params.update(info_from_filename)
        results.filename = os.path.basename(filename)

    # augment parameters with those explicitly given
    if params is not None:
        for results in results_list:
            results.params.update(params)
        
    return results_list

def slurp_res_files(
        directory_list,
        *,
        res_format=None,
        filename_format=None,
        params=None,
        glob_pattern="*.res",
        verbose=False,
):
    """Read all results files in given directories.

    The results will be a list of results data objects, one for
    each mesh point within the results file.

    Arguments:
        directory_list (str or list of str): directory or list of directories
            containing files to import
        res_format (str, optional): identifier string for the results file parser to use
        filename_format (str,optional): identifier string for the results
            filename parser to use
        glob_pattern (str,optional): glob pattern for results filenames to read
            within each directory
        verbose (bool,optional): enable debugging output

    Returns:
        (list of ResultsData): list of mesh point data objects

    """

    # process argument: upgrade single directory to list
    if (type(directory_list) == str):
        directory_list = [directory_list]
    directory_list = sorted(list(set(directory_list)))  # remove duplicate input directories
    if (verbose):
        print("  slurp_res_files: directory list {}".format(directory_list))

    # accumulate mesh points
    mesh_data = []
    for directory in directory_list:
        full_glob_pattern = os.path.join(directory, glob_pattern)
        if (verbose):
            print("  slurp_res_files: searching for files {}...".format(full_glob_pattern))
        filename_list = glob.glob(full_glob_pattern)

        # accumulate parsed data from different res files
        for filename in filename_list:
            new_mesh_data = read_file(
                filename,
                res_format=res_format,
                filename_format=filename_format,
                params=params,
                ##verbose=(verbose=="verbose_by_file")
                verbose=False,  # disabled file-by-file verbosity
            )
            mesh_data += new_mesh_data

    if (verbose):
        print("  slurp_res_files: extracted mesh points {}".format(len(mesh_data)))

    return mesh_data

################################################################
# data pickling utility
################################################################

def read_data_with_caching(read_function, pickle_filename="mesh_data.pickle", **read_function_kw):
    """Read pickled mesh data with fallback to fresh read.

    To purge cached data, manually delete pickle file.

    Arguments:

        read_function (callable): Function to read mesh data from data files

        pickle_filename (str, optional): Pickle filename

        **read_function_kw (dict, optional): Keyword arguments for read function,
        e.g., verbose=True.

    Returns:

        mesh_data (list of ResultsData): Mesh data

    """

    # attempt to read cached data
    try:
        print("Attempting to read pickled mesh data from {}...  ".format(pickle_filename), end="", flush=True)
        with open(pickle_filename, 'rb') as fp:
            mesh_data = pickle.load(fp)
        print("Done.", end="", flush=True)
        return mesh_data
    except:
        print("Failed.", flush=True)
    
    # fall back on fresh read
    print("Reading mesh_data afresh...", flush=True)
    mesh_data = read_function(**read_function_kw)

    # cache this data
    try:
        with open(pickle_filename, 'wb') as fp:
            print("Attempting to write pickled mesh data to {}...  ".format(pickle_filename), end="", flush=True)
            pickle.dump(mesh_data, fp, pickle.HIGHEST_PROTOCOL)
            print("Done.", flush=True)
    except:
        # Not sure why caching might fail, other than, say, an invalid path for
        # the pickle file?
        print("Failed.", flush=True)

    return mesh_data

#################################################
# test code                                     #
#################################################

if (__name__ == "__main__"):
    pass
