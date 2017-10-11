"""res.py

    Import control code for results files.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    5/31/15 (mac): Initiated (as mfdn_res.py).
    6/5/15 (mac): Allow user-supplied res file parser.
    6/5/15 (mac): Restructure as subpackage.
    6/29/17 (jbutler): Added in inheritance for SpNCCI, updated documentation
    7/7/17 (mac):
        - Generalize slurp_res_files to take list of directory names (after old
        analysis.import_res_files.
        - Add directory name generation utility res_file_directory.
    7/9/17 (mac):
        - Restore read_file to be simple dispatch function.
        - Extract SpNCCIMeshPointData.
    10/6/17 (mac): Extract MFDnRunData subclass to mfdn.py.
    10/10/17 (mac): Extract results data base class to results_data.py.
"""

import glob
import os

import numpy as np

#intra-packages references
import mfdnres.descriptor

################################################################
# filename utility
################################################################

def res_file_directory(username,code,run_number,results_dir="results",run_results_are_in_subdir=True):
    """Construct full path to res file directory, given user, code, and run.

        This function assumes directory naming conventions appropriate
        to mcscript archive files.

        Arguments:
            username (str): user name (e.g., "mcaprio")
            code (str): code name (e.g., "spncci")
            run_number (str): run name "tail" (e.g., "mac0424")
            results_dir (str,optional): name of top-level results directory within GROUP_HOME
            run_results_are_in_subdir (bool,optional): if results are in subdirectory "results"
              of run directory (as in an mcscript multi-task run)

        Environment:
            GROUP_HOME: directory name for group top-level results directory
              (e.g., "/afs/crc.nd.edu/group/nuclthy" for shared group results directory,
               or, for local work in your home directory, you may set equal to HOME)

        >>> res_file_directory("mcaprio","spncci","mac0417")

            /afs/crc.nd.edu/group/nuclthy/results/mcaprio/spncci/runmac0423/results

    """

    group_home = os.environ.get("GROUP_HOME")
    if (type(group_home) is not str):
        raise(ValueError("Need to set environment variable GROUP_HOME"))

    res_directory = os.path.join(group_home,results_dir,username,code,"run"+run_number)
    if (run_results_are_in_subdir):
        res_directory = os.path.join(res_directory,"results")

    return res_directory

#################################################
# parser registry
#################################################

# global registration variables
res_format_parser = {}

def register_res_format(format_name,parser):
    """Register information for parsing res file.

    Args:
        format_name (str): name for res file format
        parser (callable): function for parsing file stream

    """

    res_format_parser[format_name] = parser

##################################################
# res file import
##################################################

def read_file(filename,res_format,filename_format=None,verbose=False):
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
        res_format (str): identifier string for the results file parser to use
        filename_format (str,optional): identifier string for the results
            filename parser to use
        verbose (bool,optional): enable debugging output

    Returns:
        (list of ResultsData): list of mesh point data objects

    """

    # temporary special-case trap for legacy parsers
    #
    # These parsers populated a single MFDnRunData object.  Now
    # all parsers should return a list of results data objects.
        
    ## if (res_format == 'mfdn_v14b05'):
    ##     with open(filename, 'rt') as fin:
    ##         data = MFDnRunData()
    ##         res_format_parser[res_format](data, fin, verbose=verbose)
    ##         full_data = [data]
    ##     return full_data

    # parse results filename for any supplementary run parameters
    if (filename_format is None):
        info_from_filename = {}
    else:
        info_from_filename = mfdnres.descriptor.parse_res_filename(filename,filename_format)

    # parse results file contents for run parameters and data
    with open(filename,'rt') as fin:
        results_list = res_format_parser[res_format](fin,verbose=verbose)

    # augment parameters with those obtained from filename
    #
    # Note: The parameter values obtained from the filename will
    # *override* any parameter values obtained by parsing the results
    # file.  So beware that parameter values formatted for the
    # filename might have lower precision than those stored in the
    # results file.

    for results in results_list:
        results.params.update(info_from_filename)

    return results_list

def slurp_res_files(
        res_directory_list,res_format,
        filename_format=None,
        glob_pattern="*.res",verbose=False
):
    """Read all results file in given directories.

    The results will be a list of results data objects, one for
    each mesh point within the results file.
    
    Arguments:
        res_directory_list (str or list of str): directory or list of directories
            containing files to import
        res_format (str): identifier string for the results file parser to use
        filename_format (str,optional): identifier string for the results
            filename parser to use
        glob_pattern (str,optional): glob pattern for results filenames to read
            within each directory
        verbose (bool,optional): enable debugging output

    Returns:
        (list of ResultsData): list of mesh point data objects

    """

    # process argument: upgrade single directory to list
    if (type(res_directory_list) == str):
        res_directory_list = [res_directory_list]

    # accumulate mesh points
    mesh_data = []
    for res_directory in res_directory_list:
        full_glob_pattern = os.path.join(res_directory,glob_pattern)
        if (verbose):
            print("Searching for {}...".format(full_glob_pattern))
        res_filename_list = glob.glob(full_glob_pattern)

        # accumulate parsed data from different res files
        for res_filename in res_filename_list:
            if (verbose):
                print("Reading {}...".format(res_filename))
            new_mesh_data = read_file(
                res_filename,
                res_format=res_format,filename_format=filename_format,
                verbose=verbose
            )
            if (verbose):
                print("  {:d} mesh points".format(len(new_mesh_data)))
            mesh_data += new_mesh_data

    return mesh_data

#################################################
# test code                                     #
#################################################

if (__name__ == "__main__"):
    pass
