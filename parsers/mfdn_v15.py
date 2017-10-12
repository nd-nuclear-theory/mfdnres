""" res_parser_mfdn_v15

    Provides parser for mfdn v15 results files.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    10/6/17 (mac): Created, from res_parser_spncci.py
"""

import itertools

import numpy as np

import mfdnres.tools
import mfdnres.mfdn_results_data


################################################################
# section handlers
################################################################

def parse_params(self,tokenized_lines):
    """
    Parse any section containing key-value pairs to add to params dictionary.

    TODO: extract "observables" from numbered TBMEfile entries
    """

    # extract key-value pairs
    conversions = {
        # MFDn -- only selected numerical fields
        "ndiags" : mfdnres.tools.singleton_of(int),
        "MPIranks" : mfdnres.tools.singleton_of(int),
        "OMPthreads" : mfdnres.tools.singleton_of(int),
        # Basis
        "Nprotons" : mfdnres.tools.singleton_of(int),
        "Nneutrons" : mfdnres.tools.singleton_of(int),
        "TwoMj" : mfdnres.tools.singleton_of(int),
        "parity" : mfdnres.tools.singleton_of(int),
        "Nmin" : mfdnres.tools.singleton_of(int),
        "Nmax" : mfdnres.tools.singleton_of(int),
        "DeltaN" : mfdnres.tools.singleton_of(int),
        "WTmax" : mfdnres.tools.singleton_of(float),
        # Many-body matrix
        "dimension" : mfdnres.tools.singleton_of(int),
        "numnonzero" : mfdnres.tools.singleton_of(int),
        # Interaction
        "Hrank" : mfdnres.tools.singleton_of(int),
        "hbomeg" : mfdnres.tools.singleton_of(float),
        "fmass" : mfdnres.tools.singleton_of(float),
        # Observables
        "numTBops" : mfdnres.tools.singleton_of(int),
        # Calculation
        "hw" : mfdnres.tools.singleton_of(float)
    }
    key_value_dict = mfdnres.tools.extract_key_value_pairs(
        tokenized_lines,conversions
    )

    # augment with float Mj (from TwoMj)
    if ("TwoMj" in key_value_dict):
        key_value_dict["Mj"] = key_value_dict["TwoMj"]/2

    # augment with alias hw (from hbomeg)
    if ("hbomeg" in key_value_dict):
        key_value_dict["hw"] = key_value_dict["hbomeg"]

    # augment with nuclide and mass number
    if (("Nprotons" in key_value_dict) and ("Nneutrons" in key_value_dict)):
        key_value_dict["nuclide"] = (key_value_dict["Nprotons"],key_value_dict["Nneutrons"])
        key_value_dict["A"] = key_value_dict["Nprotons"]+key_value_dict["Nneutrons"]

    # update to params dictionary
    self.params.update(key_value_dict)

def parse_energies(self,tokenized_lines):
    """ Parse energies.
    """

    # import energy tabulation
    table = np.array(
        tokenized_lines,
        dtype=[
            ("seq",int),("J",float),("n",int),("T",float),
            ("E",float),  # "Eabs" in section header comment
            ("Eex",float),  # "Eexc" in section header comment
            ("Error",float),("J-full",float)
        ]
        )

    # sort energies into energies dictionary
    #
    # how do we want to deal with parity label, or lack thereof?
    for entry in table:
        (_,J,n,_,E,_,_,_)=entry
        ## self.energies[(J,gex,n)]=E
        ## self.num_eigenvalues[(J,gex)]=self.num_eigenvalues.setdefault((J,gex),0)+1
        self.energies[(J,n)]=E
        self.num_eigenvalues[(J,)]=self.num_eigenvalues.setdefault((J,),0)+1

section_handlers = {
    # [CODE]
    "MFDn" : parse_params,
    # PARAMETERS
    "Basis" : parse_params,
    "Many-body matrix" : parse_params,
    "Interaction" : parse_params,
    "Observables" : parse_params,
    # RESULTS
    "Energies" : parse_energies,
    ## "Decompositions: Nex" : parse_decompositions_Nex,
    ## "Decompositions: BabySpNCCI" : parse_decompositions_baby_spncci,
    ## "Observable RMEs" : parse_observable_rmes
}

################################################################
# parsing control code
################################################################

def parse_mesh_point(self,sections,section_handlers):
    """ Parse single mesh point into results object.

    Arguments:
        self (ResultsData): results object to populate
        sections (list): data for sections to parse, as (section_name,tokenized_lines) tuples
        section_handlers (dict): dictionary of section handlers
    """

    for (section_name,tokenized_lines) in sections:
        if section_name in section_handlers:
            section_handlers[section_name](self,tokenized_lines)

def parser(in_file,verbose):
    """ Parse full spncci results file, into list of one or more results objects.

    Arguments:
        in_file (stream): input file stream (already opened by caller)
        verbose (bool,optional): enable verbose output
    """

    # perform high-level parsing into sections
    res_file_lines = [row for row in in_file]
    tokenized_lines = mfdnres.tools.split_and_prune_lines(res_file_lines)
    sections = mfdnres.tools.extracted_sections(tokenized_lines)

    # set up container
    results = mfdnres.mfdn_results_data.MFDnResultsData()

    # parse sections
    parse_mesh_point(results,sections,section_handlers)

    # package results
    mesh_data = [results]
    return mesh_data

# register the parser
mfdnres.res.register_res_format('mfdn_v15',parser)


if (__name__=="__main__"):
    pass
