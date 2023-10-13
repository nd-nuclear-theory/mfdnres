""" mfdn_v15

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
    09/17/20 (zz): Add one-body observable parser.
    09/17/20 (mac): Update data attribute names.  Add two-body obervable parser.
    09/25/21 (mac): Parse extended relative radii observables.
    10/12/21 (pjf): Warn if file is empty.
    08/01/22 (pjf): Use RMEData for RME storage, storing operator quantum numbers.
    01/11/23 (mac): Parse Diagonalization parameters and truncate parsing to states within requested neivals.
    01/15/23 (mac): Support res file output from MFDn menj variant.
    09/30/23 (mac):
        - Provide more informative error output on parsing failure, indicating which input file section yields parsing error.
        - Update handling of angular momenta for mfdn GPU version.
    10/05/23 (mac): Update handling of ambiguous radius observables.
"""

from __future__ import annotations

import itertools
import traceback

import numpy as np

from .. import (
    results_data,
    mfdn_results_data,
    input,
    tools,
    )

from ..mfdn_results_data import (
    MFDnResultsData,  # for typing
)


################################################################
# global state :(
################################################################

# k_parameter_g (int): parity grade of run (0 for positive parity, 1
#     for negative parity); required for generating (J,g,n) quantum
#     number labels

k_parameter_g:int = None


################################################################
# parsing utilities
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

            
def count_mfdn_results_line_properties(self:MFDnResultsData,tokenized_lines):
    """Peek at generic mfdn results line, to count data properties.

    This is necessary when columns are added to MFDn's generic static properties
    output without notice (except perhaps in a the human-only hash line
    comment).

    Arguments:
        ...

    Returns:
        (int); number of columns, or None if section is empty
    """

    # trap empty section
    if len(tokenized_lines)==0:
        return None
    
    # count properties
    tokenized_line = tokenized_lines[0]
    (qn,data) = split_mfdn_results_line(tokenized_line)
    return len(data)


################################################################
# section handlers
################################################################

def parse_params(self:MFDnResultsData,tokenized_lines):
    """
    Parse any section containing key-value pairs to add to params dictionary.

    Includes special handling of some keys:

      - Add derived keys to dictionary:

        "M", "hw", "nuclide", "A"

      - Store global state based on key:

        k_parameter_g

    Globals:
        k_parameter_g (output)

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
        # Diagonalization
        "neivals" : tools.singleton_of(int),
        "maxits" : tools.singleton_of(int),
        "startit" : tools.singleton_of(int),
        "selectpiv" : tools.singleton_of(int),
        "tol" : tools.singleton_of(float),
        # Many-body matrix
        "dimension" : tools.singleton_of(int),
        "numnonzero" : tools.singleton_of(int),
        # Interaction
        "Hrank" : tools.singleton_of(int),
        "hbomeg" : tools.singleton_of(float),
        "fmass" : tools.singleton_of(float),
        # Interaction -- menj variant
        "EMax" : tools.singleton_of(int),
        "MEID" : tools.singleton_of(str),
        "E3Max" : tools.singleton_of(int),
        "ME3ID" : tools.singleton_of(str),
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

    ## self.params["neivals"] = None
    
def parse_energies(self:MFDnResultsData,tokenized_lines):
    """ Parse energies.

    Globals:
        k_parameter_g (input)
    """

    # truncate to requested states
    neivals = self.params.get("neivals")
    tokenized_lines = tokenized_lines[:neivals]
    
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
    self.mfdn_level_properties["T"] = {}

    # sort energies into energies dictionary
    for entry in table:

        # retrieve entry
        (_,J,n,T,E,error,_,_)=entry
        g = k_parameter_g
        qn = (float(J),int(g),int(n))  # 20/07/08 (pjf): Cast explicitly to Python types.
        Jg_pair = (J,g)

        # store data from entry
        self.energies[qn] = E
        self.mfdn_level_residuals[qn] = error
        self.mfdn_level_properties["T"][qn] = T

        
def parse_decompositions_Nex(self:MFDnResultsData,tokenized_lines):
    """Parse Nex decomposition.

    Unstable feature: Note only alternate Nex values are output by
    MFDn (assuming Nstep=2).  Currently, these only these values are
    stored, verbatim, with no zero padding for the missing values.
    Beware this behavior is different from the spncci parser, for
    which all Nex values (step 1) are interleaved, so that the vector
    index is Nex.  Consider interleaving 0 values here to make this
    true for MFDn parser as well.

    """

    # truncate to requested states
    neivals = self.params.get("neivals")
    tokenized_lines = tokenized_lines[:neivals]

    # store decompositions
    self.mfdn_level_decompositions["Nex"] = {}
    for tokenized_line in tokenized_lines:
        (qn,data) = split_mfdn_results_line(tokenized_line)
        self.mfdn_level_decompositions["Nex"][qn]=data

        
def parse_generic_static_properties(self:MFDnResultsData,tokenized_lines,container,property_names):
    """Parse generic static properties given list of property names for the data columns.

    Any "extra" values in the parsed line, beyond the number of property names
    given, are silently ignored.

    Arguments:

        ...

        container (dict): dictionary to which properties should be added

        property_names (list of str): names for these properties (or None for a field to be ignored)

    """

    # truncate to requested states
    neivals = self.params.get("neivals")
    tokenized_lines = tokenized_lines[:neivals]

    # set up containers for properties
    for property_name in property_names:
        container[property_name] = {}

    # store properties
    for tokenized_line in tokenized_lines:
        (qn,data) = split_mfdn_results_line(tokenized_line)
        for property_index in range(len(property_names)):
            property_name = property_names[property_index]
            if property_name is not None:
                container[property_name][qn]=data[property_index]

    
def parse_M1_moments(self:MFDnResultsData,tokenized_lines):
    """Parse M1 moments.
    """
    property_names = ["M1","Dlp","Dln","Dsp","Dsn"]
    parse_generic_static_properties(self,tokenized_lines,self.mfdn_ob_moments,property_names)

    
def parse_E2_moments(self:MFDnResultsData,tokenized_lines):
    """Parse E2 moments.
    """
    property_names = ["E2p","E2n"]
    parse_generic_static_properties(self,tokenized_lines,self.mfdn_ob_moments,property_names)

    
def parse_angular_momenta(self:MFDnResultsData,tokenized_lines):
    """Parse squared angular momenta.
    """

    # parse raw (squared observable) values
    revision = self.params.get("Revision","beta00")
    num_properties = count_mfdn_results_line_properties(self,tokenized_lines)
    if num_properties==5:
        # mfdn v15 beta00
        property_names = ["L_sqr", "S_sqr", "Sp_sqr", "Sn_sqr", "J_sqr"]
    elif num_properties==6:
        # 09/30/23 (mac): GPU version also has T_sqr (e.g., "v15b01-205-gdb2400a")
        property_names = ["L_sqr", "S_sqr", "Sp_sqr", "Sn_sqr", "J_sqr", "T_sqr"]
    elif num_properties==7:
        # v15beta01 added nonsensical Lp^2 and Ln^2 observables (unclear definition and sometimes negative)
        #
        #    ["L_sqr", "S_sqr", "Lp_sqr", "Sp_sqr", "Ln_sqr", "Sn_sqr", "J_sqr"]
        #
        # These nonsensical results should therefore be *ignored* on input, to
        # prevent later misuse.
        property_names = ["L_sqr", "S_sqr", None, "Sp_sqr", None, "Sn_sqr", "J_sqr"]
    else:
        raise ValueError("unexpected number of columns in angular momentum section ({} observables)".format(num_properties))
    parse_generic_static_properties(self,tokenized_lines,self.mfdn_level_properties,property_names)

    
def parse_radii(self:MFDnResultsData,tokenized_lines):
    """Parse radii.
    """
    # Originally, there were just three radius columns:
    #
    #     ["rp", "rn", "r"]
    #
    # An extended set of observables was introduced with some version of
    # MFDn v15beta01 (<= v15b01-39-g18c712b):
    #
    #     ["rp", "rn", "r", "rpp", "rnn", "rpn"]
    #
    # See cshalo [PRC 90, 034305 (2014)] (A5) for definitions of rpp, rnn,
    # and rpn.
    #
    # Then an extra, untitled 11th column (populated with numerical zeros)
    # appears to have been added sometime <=v15b01-37-gb46e061.
    #
    # 10/02/23 (mac): However, pm indicates different versions of MFDn have
    # output ["rpp", "rnn", "rpn"] or ["rpp", "rpn", "rnn"] for the last three
    # columns, so interpretation of the last two columns is ambiguous.  These
    # last two columns should therefore be *ignored* on input, to prevent later misuse.
    num_properties = count_mfdn_results_line_properties(self,tokenized_lines)
    if num_properties==3:
        property_names = ["rp", "rn", "r"]
    elif num_properties in [6, 7]:
        property_names = ["rp", "rn", "r", "rpp"]
    else:
        raise ValueError("unrecognized number of columns in radii section ({} observables)".format(num_properties))
    parse_generic_static_properties(self,tokenized_lines,self.mfdn_tb_expectations,property_names)

    
def parse_other_tbo(self:MFDnResultsData,tokenized_lines):
    """Parse other two-body observables.

    Requires "tbo_names" to have been parsed from TBMEfile entries in MFDn
    output (as for MFDn h2).  Otherwise quietly skips parsing "other TBOs" (as
    for MFDn menj variant).

    """

    if "tbo_names" in self.params:
        property_names = self.params["tbo_names"][1:]
    else:
        property_names = []
    parse_generic_static_properties(self,tokenized_lines,self.mfdn_tb_expectations,property_names)

    
def parse_mfdn_ob_rmes(self:MFDnResultsData, tokenized_lines):
    """Parse native transition output.

    Note: Native transitions were removed long before neivals parameter was
    provided in MFDn res file output, so these need not be truncated to neivals
    states.
    """
    transition_classes = {
        "GT": ["GTi(Z+1,N-1)", "GTi(Z-1,N+1)", "GTf(Z+1,N-1)", "GTf(Z-1,N+1)"],
        "M1": ["M1", "Dlp", "Dln", "Dsp", "Dsn"],
        "E2": ["E2p", "E2n"],
    }
    transition_type_qn:dict[str,tools.OperatorQNType] = {
        "GTi(Z+1,N-1)": (1,0,-1), "GTi(Z-1,N+1)": (1,0,+1),
        "GTf(Z+1,N-1)": (1,0,+1), "GTf(Z-1,N+1)": (1,0,-1),
        "M1": (1,0,0), "Dlp": (1,0,0), "Dln": (1,0,0), "Dsp": (1,0,0), "Dsn": (1,0,0),
        "E2p": (2,0,0), "E2n": (2,0,0),
    }

    for tokenized_line in tokenized_lines:
        qnf = (float(tokenized_line[1]), k_parameter_g, int(tokenized_line[2]))
        transition_class = tokenized_line[3]
        qni = (float(tokenized_line[5]), k_parameter_g, int(tokenized_line[6]))
        data_iterable = list(map(float, tokenized_line[7:]))
        data = np.array(data_iterable, dtype=float)

        # break out each component as a separate transition type
        transition_types = transition_classes[transition_class]
        for (transition_type, value) in zip(transition_types, data):
            operator_qn = transition_type_qn[transition_type]
            transition_dict = self.mfdn_ob_rmes.setdefault(
                transition_type,
                results_data.RMEData(
                    qn=operator_qn,
                    rme_convention=tools.RMEConvention.kEdmonds,
                )
            )
            transition_dict[(qnf,qni)] = value


def parse_postprocessor_ob_rmes_legacy(self:MFDnResultsData, tokenized_lines):
    """Parse obscalc-ob output for transitions (legacy tabular format).

    LEGACY: Old obscalc-ob format with multiple RMEs per line.
    """
    transition_type_qn:dict[str,tools.OperatorQNType] = {
        "Dlp": (1,0,0), "Dln": (1,0,0), "Dsp": (1,0,0), "Dsn": (1,0,0),
        "E2p": (2,0,0), "E2n": (2,0,0),
    }

    names = self.params["one_body_observable_names"]
    for tokenized_line in tokenized_lines:
        qnf = (float(tokenized_line[0]), int(tokenized_line[1]), int(tokenized_line[2]))
        qni = (float(tokenized_line[3]), int(tokenized_line[4]), int(tokenized_line[5]))
        data_iterable = map(float, tokenized_line[6:])

        for (transition_type, value) in zip(names, data_iterable):
            try:
                operator_qn = transition_type_qn[transition_type]
            except KeyError:
                raise ValueError("unsupported observable {} in legacy postprocessor file".format(transition_type))
            transition_dict = self.postprocessor_ob_rmes.setdefault(
                transition_type,
                results_data.RMEData(
                    qn=operator_qn,
                    rme_convention=tools.RMEConvention.kEdmonds,
                )
            )
            transition_dict[(qnf,qni)] = value


def parse_postprocessor_generic_rmes(self:MFDnResultsData,tokenized_lines,container:dict[str,results_data.RMEData]):
    """Parse generic (ob or tb) postprocessor rmes.

    Arguments:
        ...
        container (dict): dictionary to which rmes should be added
    """
    tokenized_line = tokenized_lines[0]
    (J0, g0, Tz0, name) = (int(tokenized_line[0]), int(tokenized_line[1]), int(tokenized_line[2]), tokenized_line[3])
    rme_dict = container.setdefault(
        name,
        results_data.RMEData(qn=(J0,g0,Tz0), rme_convention=tools.RMEConvention.kEdmonds)
    )
    for tokenized_line in tokenized_lines[1:]:
        qnf = (float(tokenized_line[0]), int(tokenized_line[1]), int(tokenized_line[2]))
        qni = (float(tokenized_line[3]), int(tokenized_line[4]), int(tokenized_line[5]))
        rme = float(tokenized_line[6])

        rme_dict[(qnf,qni)] = rme


def parse_postprocessor_ob_rmes(self:MFDnResultsData, tokenized_lines):
    """Parse obscalc-ob output for a one-body observable.
    """
    parse_postprocessor_generic_rmes(self, tokenized_lines, self.postprocessor_ob_rmes)

    
def parse_postprocessor_tb_rmes(self:MFDnResultsData, tokenized_lines):
    """Parse postprocessor output (as digested by the scripting) for a two-body observable.
    """
    parse_postprocessor_generic_rmes(self, tokenized_lines, self.postprocessor_tb_rmes)

section_handlers = {
    # [CODE]
    "MFDn" : parse_params,
    # PARAMETERS
    "Basis" : parse_params,
    "Diagonalization" : parse_params,
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
    "Transitions": parse_mfdn_ob_rmes,
    "Transition one-body observables": parse_postprocessor_ob_rmes_legacy,
    ## "Transition one-body observables": parse_one_body_static_properties,  # WRONG?
    "One-body observable": parse_postprocessor_ob_rmes,
    "Two-body observable": parse_postprocessor_tb_rmes,
}


################################################################
# parsing control code
################################################################

def parse_mesh_point(self:MFDnResultsData, sections, section_handlers):
    """ Parse single mesh point into results object.

    Arguments:
        self (ResultsData): results object to populate
        sections (list): data for sections to parse, as (section_name,tokenized_lines) tuples
        section_handlers (dict): dictionary of section handlers
    """

    for (section_name,tokenized_lines) in sections:
        if section_name in section_handlers:
            try:
                section_handlers[section_name](self,tokenized_lines)
            except Exception as err:
                print("ERROR: Unexpected content in results file section '{}'".format(section_name))
                raise err
            
def parser(in_file, verbose):
    """ Parse full results file.

    Arguments:
        in_file (stream): input file stream (already opened by caller)
        verbose (bool,optional): enable verbose output
    """

    # perform high-level parsing into sections
    res_file_lines = [row for row in in_file]
    tokenized_lines = tools.split_and_prune_lines(res_file_lines)
    sections = tools.extracted_sections(tokenized_lines)

    # handle empty files
    if len(res_file_lines)==0:
        print("WARNING: file {} is empty!".format(in_file.name))
        return []

    # set up container
    results = mfdn_results_data.MFDnResultsData()

    # initialize data attributes for "global" use within parser
    ## results._parameter_g:int = None  # TODO (mac): use to replace global k_parameter_g
    
    # parse sections
    parse_mesh_point(results,sections,section_handlers)

    # package results
    mesh_data = [results]
    return mesh_data


# register the parser
input.register_data_format('mfdn_v15', parser)
input.register_code_name('mfdn15', 'mfdn_v15')
input.register_code_name('obscalc', 'mfdn_v15')
input.register_code_name('obscalc-ob', 'mfdn_v15')
input.register_code_name('transitions-ob', 'mfdn_v15')
input.register_code_name('transitions-tb', 'mfdn_v15')


if (__name__=="__main__"):
    pass
