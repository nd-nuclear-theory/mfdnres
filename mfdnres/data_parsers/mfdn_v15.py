""" res_parser_mfdn_v15

    Provides parser for mfdn v15 results files.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    10/06/17 (mac): Created, from res_parser_spncci.py
    04/27/18 (mac): Rename parameter Mj to M.
    09/06/18 (pjf): Added initial built-in transition extraction.
    12/12/18 (mac): Update handling of "Angular momenta" section for v15b01.
    02/22/19 (pjf): Handle early MFDn v15b00 files.
    04/02/19 (mac): Update naming of angular momenta observables.
    04/05/19 (pjf): Parse obscalc-ob one-body static and transition output.
    10/10/19 (pjf): Store Lanczos residuals.
    01/04/20 (mac): Remove num_eigenvalues as static data (now available as property).
    06/17/20 (pjf): Add code registration.
    07/08/20 (pjf): Fix quantum number types.
    09/17/20 (zz): Add one body observable parser.
"""

import itertools

import numpy as np

from .. import (
    mfdn_results_data,
    input,
    tools,
    )


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

    # recover qn
    qn = (float(tokenized_line[1]),k_parameter_g,int(tokenized_line[2]))

    # trap overflow values
    #
    # NAIVE: data_list = list(map(float,tokenized_line[4:]))
    #
    # Note: FORTRAN may output "NaN" or "*****".  The former is handled
    # gracefully by float as float("NaN") => nan, but float("*****") crashes, so
    # we must trap this case.
    data_list = [
        float(token) if (token[0]!="*") else np.nan
        for token in tokenized_line[4:]
    ]

    data = np.array(data_list,dtype=float)
    return (qn,data)


################################################################
# section handlers
################################################################

def parse_params(self,tokenized_lines):
    """
    Parse any section containing key-value pairs to add to params dictionary.

    Includes special handling of some keys:

      - Add derived keys to dictionary:

        "M", "hw", "nuclide", "A"

      - Store global state based on key:

        k_parameter_g

    Globals:
        k_parameter_g (output)

    TODO: extract "observables" from numbered TBMEfile entries
    """

    # extract key-value pairs
    conversions = {
        # MFDn
        # not yet parsed: Platform, Username
        "Version" : tools.singleton_of(int),
        "Revision" : tools.singleton_of(str),
        "ndiags" : tools.singleton_of(int),
        "MPIranks" : tools.singleton_of(int),
        "OMPthreads" : tools.singleton_of(int),
        # Basis
        "Nprotons" : tools.singleton_of(int),
        "Nneutrons" : tools.singleton_of(int),
        "TwoMj" : tools.singleton_of(int),
        "TwoM" : tools.singleton_of(int),
        "parity" : tools.singleton_of(int),
        "Nmin" : tools.singleton_of(int),
        "Nmax" : tools.singleton_of(int),
        "DeltaN" : tools.singleton_of(int),
        "WTmax" : tools.singleton_of(float),
        # Many-body matrix
        "dimension" : tools.singleton_of(int),
        "numnonzero" : tools.singleton_of(int),
        # Interaction
        "Hrank" : tools.singleton_of(int),
        "hbomeg" : tools.singleton_of(float),
        "fmass" : tools.singleton_of(float),
        # Observables
        "numTBops" : tools.singleton_of(int),
        "TBMEfile" : tools.singleton_of(str),
        "names": tools.list_of(str),
        # Calculation
        "hw" : tools.singleton_of(float)
    }
    key_value_dict = tools.extract_key_value_pairs(
        tokenized_lines,conversions
    )

    # do special handling of keys

    # augment with float M (from TwoM)
    if ("TwoMj" in key_value_dict):
        key_value_dict["M"] = key_value_dict["TwoMj"]/2
    if ("TwoM" in key_value_dict):
        key_value_dict["M"] = key_value_dict["TwoM"]/2

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

    if ("TBMEfile" in key_value_dict):
        key_value_dict["tbo_names"] = [filename.replace('.bin', '').replace('tbme-', '') for filename in key_value_dict["TBMEfile"]]

    if "names" in key_value_dict:
        key_value_dict["one_body_observable_names"] = key_value_dict["names"]

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

    # set up container for Lanczos residuals
    self.residuals = {}

    # sort energies into energies dictionary
    for entry in table:

        # retrieve entry
        (_,J,n,T,E,error,_,_)=entry
        g = k_parameter_g
        qn = (float(J),int(g),int(n))  # 20/07/08 (pjf): Cast explicitly to Python types.
        Jg_pair = (J,g)

        # store data from entry
        self.energies[qn] = E
        self.residuals[qn] = error
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
    """Parse generic static properties given list of property names for the data columns.

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
    property_names = ["M1","Dlp","Dln","Dsp","Dsn"]
    parse_generic_static_properties(self,tokenized_lines,self.native_static_properties,property_names)

def parse_E2_moments(self,tokenized_lines):
    """Parse E2 moments.
    """
    property_names = ["E2p","E2n"]
    parse_generic_static_properties(self,tokenized_lines,self.native_static_properties,property_names)

def parse_angular_momenta(self,tokenized_lines):
    """Parse squared angular momenta.
    """

    # parse raw (squared observable) values
    if (self.params.get("Revision","beta00")=="beta00"):  # early runs with beta00 left "Revision" field blank
        property_names = ["L_sqr","S_sqr","Sp_sqr","Sn_sqr","J_sqr"]
    else:
        # nonsensical Lp^2 and Ln^2 observables (unclear definition and sometimes negative) added with v15beta01
        property_names = ["L_sqr","S_sqr","Lp_sqr","Sp_sqr","Ln_sqr","Sn_sqr","J_sqr"]
    parse_generic_static_properties(self,tokenized_lines,self.native_static_properties,property_names)

def parse_radii(self,tokenized_lines):
    """Parse radii.
    """
    property_names = ["rp","rn","r"]
    parse_generic_static_properties(self,tokenized_lines,self.two_body_static_observables,property_names)

def parse_other_tbo(self,tokenized_lines):
    """Parse other two-body observables.
    """
    property_names = self.params["tbo_names"][1:]
    parse_generic_static_properties(self,tokenized_lines,self.two_body_static_observables,property_names)

def parse_transitions(self, tokenized_lines):
    """Parse native transition output.
    """
    transition_classes = {
        "GT": ["GTi(Z+1,N-1)", "GTi(Z-1,N+1)", "GTf(Z+1,N-1)", "GTf(Z-1,N+1)"],
        "M1": ["M1", "Dlp", "Dln", "Dsp", "Dsn"],
        "E2": ["E2p", "E2n"],
    }

    for tokenized_line in tokenized_lines:
        qnf = (float(tokenized_line[1]), k_parameter_g, int(tokenized_line[2]))
        transition_class = tokenized_line[3]
        qni = (float(tokenized_line[5]), k_parameter_g, int(tokenized_line[6]))
        data_iterable = list(map(float, tokenized_line[7:]))
        data = np.array(data_iterable, dtype=float)

        # only store canonical pair
        (Jgn_pair_canonical, flipped, canonicalization_factor) = tools.canonicalize_Jgn_pair(
            (qnf, qni), tools.RMEConvention.kEdmonds
        )

        # break out each component as a separate transition type
        transition_types = transition_classes[transition_class]
        for (transition_type, value) in zip(transition_types, data):
            transition_dict = self.native_transition_properties.setdefault(transition_type, dict())
            transition_dict[Jgn_pair_canonical] = canonicalization_factor*value

def parse_one_body_static_properties(self, tokenized_lines):
    """Parse obscalc-ob output for static properties.
    """
    names = self.params["one_body_observable_names"]
    for tokenized_line in tokenized_lines:
        qn = (float(tokenized_line[0]), int(tokenized_line[1]), int(tokenized_line[2]))
        data_iterable = map(float, tokenized_line[6:])

        for (property_type, value) in zip(names, data_iterable):
            property_dict = self.one_body_static_properties.setdefault(property_type, dict())
            property_dict[qn] = value


def parse_one_body_transitions(self, tokenized_lines):
    """Parse obscalc-ob output for transitions.
    """
    names = self.params["one_body_observable_names"]
    for tokenized_line in tokenized_lines:
        qnf = (float(tokenized_line[0]), int(tokenized_line[1]), int(tokenized_line[2]))
        qni = (float(tokenized_line[3]), int(tokenized_line[4]), int(tokenized_line[5]))
        data_iterable = map(float, tokenized_line[6:])

        # only store canonical pair
        (Jgn_pair_canonical, flipped, canonicalization_factor) = tools.canonicalize_Jgn_pair(
            (qnf, qni), tools.RMEConvention.kEdmonds
        )

        for (transition_type, value) in zip(names, data_iterable):
            transition_dict = self.one_body_transition_properties.setdefault(transition_type, dict())
            transition_dict[Jgn_pair_canonical] = canonicalization_factor*value


def parse_one_body_observable(self, tokenized_lines):
    """Parse obscalc-ob output for an observable
    """
    (J0, g0, Tz0, name) = tokenized_lines[0]
    transition_dict = self.one_body_transition_properties.setdefault(name, dict())
    for tokenized_line in tokenized_lines[1:]:
        qnf = (float(tokenized_line[0]), int(tokenized_line[1]), int(tokenized_line[2]))
        qni = (float(tokenized_line[3]), int(tokenized_line[4]), int(tokenized_line[5]))
        rme = float(tokenized_line[6])

        # only store canonical pair
        (Jgn_pair_canonical, flipped, canonicalization_factor) = tools.canonicalize_Jgn_pair(
            (qnf, qni), tools.RMEConvention.kEdmonds
        )

        transition_dict[Jgn_pair_canonical] = canonicalization_factor*rme


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
    "Other 2-body observables" : parse_other_tbo,
    "Transitions": parse_transitions,
    "Transition one-body observables": parse_one_body_transitions,
    ## "Transition one-body observables": parse_one_body_static_properties,  # WRONG?
    "One-body observable": parse_one_body_observable,

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
    tokenized_lines = tools.split_and_prune_lines(res_file_lines)
    sections = tools.extracted_sections(tokenized_lines)

    # set up container
    results = mfdn_results_data.MFDnResultsData()

    # parse sections
    parse_mesh_point(results,sections,section_handlers)

    # package results
    mesh_data = [results]
    return mesh_data

# register the parser
input.register_data_format('mfdn_v15',parser)
input.register_code_name('mfdn15', 'mfdn_v15')
input.register_code_name('obscalc', 'mfdn_v15')
input.register_code_name('obscalc-ob', 'mfdn_v15')


if (__name__=="__main__"):
    pass
