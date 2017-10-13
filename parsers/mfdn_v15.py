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
# global state :(
################################################################

# k_parameter_g (int): parity grade of run (0 for positive parity, 1
#     for negative parity); required for generating (J,g,n) quantum
#     number labels

k_parameter_g = None

################################################################
# parsing utility
################################################################

def split_mfdn_results_line(tokenized_line):
    """ Recover qn and numerical data from standard MFDn output line.

    The isospin field T is discarded.

    Standard line format:
        # Seq    J    n      T        <data>
            1   1.5   1   0.500      0.7949      0.1272      0.7793E-01

    This line would yield (assuming g=1 for illustration):
        (
            (1.5, 1, 1),
            np.array([ 0.7949   0.1272   0.07793])
        )

    Arguments:
        tokenized_line (list of str): tokens in line

    Globals:
        k_parameter_g (input)

    Returns:
        qn (tuple): (J,g,n)
        data (np.array): floating point data as np vector
    """
    qn = (float(tokenized_line[1]),k_parameter_g,int(tokenized_line[2]))
    data_iterable = list(map(float,tokenized_line[4:]))
    data = np.array(data_iterable,dtype=float)
    return (qn,data)


################################################################
# section handlers
################################################################

def parse_params(self,tokenized_lines):
    """
    Parse any section containing key-value pairs to add to params dictionary.

    Includes special handling of some keys:

      - Add derived keys to dictionary:

        "Mj", "hw", "nuclide", "A"

      - Store global state based on key:

        k_parameter_g

    Globals:
        k_parameter_g (output)

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

    # do special handling of keys

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

    # store global parity grade
    if ("parity" in key_value_dict):
        global k_parameter_g
        k_parameter_g = (1-key_value_dict["parity"])//2

    # update to params dictionary
    self.params.update(key_value_dict)

def parse_energies(self,tokenized_lines):
    """ Parse energies.

    Globals:
        k_parameter_g (input)
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

    # set up container for native-calculated isospins
    self.native_static_properties["T"] = {}

    # sort energies into energies dictionary
    for entry in table:

        # retrieve entry
        (_,J,n,T,E,_,_,_)=entry
        g = k_parameter_g
        qn = (J,g,n)
        Jg_pair = (J,g)

        # store data from entry
        self.energies[qn] = E
        self.num_eigenvalues[Jg_pair] = self.num_eigenvalues.setdefault(Jg_pair,0)+1
        self.native_static_properties["T"][qn] = T

def parse_decompositions_Nex(self,tokenized_lines):
    """Parse Nex decomposition.

    Unstable feature: Note only alternate Nex values are output by
    MFDn (assuming Nstep=2).  Currently, these only these values are
    stored, verbatim, with no zero padding for the missing values.
    Beware this behavior is different from the spncci parser, for
    which all Nex values (step 1) are interleaved, so that the vector
    index is Nex.  Consider interleaving 0 values here to make this
    true for MFDn parser as well.

    """
    self.decompositions["Nex"] = {}
    for tokenized_line in tokenized_lines:
        (qn,data) = split_mfdn_results_line(tokenized_line)
        self.decompositions["Nex"][qn]=data

def parse_generic_static_properties(self,tokenized_lines,container,property_names):
    """Parse generic moments given list of property names for the data columns.

    Arguments:
        ...
        container (dict): dictionary to which properties should be added
        property_names (list of str): names for these properties
    """

    for property_name in property_names:
        container[property_name] = {}

    for tokenized_line in tokenized_lines:
        (qn,data) = split_mfdn_results_line(tokenized_line)
        for property_index in range(len(property_names)):
            property_name = property_names[property_index]
            container[property_name][qn]=data[property_index]

def parse_M1_moments(self,tokenized_lines):
    """Parse M1 moments.
    """
    property_names = ["M1mu","M1lp","M1ln","M1sp","M1sn"]
    parse_generic_static_properties(self,tokenized_lines,self.native_static_properties,property_names)

def parse_E2_moments(self,tokenized_lines):
    """Parse E2 moments.
    """
    property_names = ["E2p","E2n"]
    parse_generic_static_properties(self,tokenized_lines,self.native_static_properties,property_names)

def parse_angular_momenta(self,tokenized_lines):
    """Parse angular momenta.
    """
    property_names = ["L","S","Sp","Sn","J"]
    parse_generic_static_properties(self,tokenized_lines,self.native_static_properties,property_names)

def parse_radii(self,tokenized_lines):
    """Parse radii.
    """
    property_names = ["rp","rn","r"]
    parse_generic_static_properties(self,tokenized_lines,self.two_body_static_observables,property_names)

def parse_other_tbo(self,tokenized_lines):
    """Parse other two-body observables.  (WIP)
    """
    for tokenized_line in tokenized_lines:
        (qn,data) = split_mfdn_results_line(tokenized_line)
        # TODO
        pass

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
    "Oscillator quanta" : parse_decompositions_Nex,
    "M1 moments" : parse_M1_moments,
    "E2 moments" : parse_E2_moments,
    "Angular momenta" : parse_angular_momenta,
    "Relative radii" : parse_radii,
    "Other 2-body observables" : parse_other_tbo

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
    """ Parse full results file.

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
