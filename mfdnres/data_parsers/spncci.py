""" spncci

    Provides parser for spncci results files.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    07/09/17 (mac): Created.
    09/17/17 (mac): Fix bool conversion on input.
    01/04/20 (mac): Remove num_eigenvalues as static data (now available as property).
    05/21/20 (mac): Fix canonicalization factor for storage in parse_observable_rmes.
"""

import itertools

import numpy as np

from .. import (
    spncci_results_data,
    input,
    tools,
    )

################################################################
# section handlers
################################################################

def parse_params(self,tokenized_lines):
    """
    Parse any section containing key-value pairs to add to params dictionary.
    """

    # extract key-value pairs
    conversions = {
        # Space
        "nuclide" : tools.tuple_of(int),  # use tuple so parameter is hashable when used as analysis key
        "A" : tools.singleton_of(int),
        "Nsigma" : tools.singleton_of(float),
        "Nsigmamax" : tools.singleton_of(int),
        "N1v" : tools.singleton_of(int),
        "Nmax" : tools.singleton_of(int),
        # Interaction
        "interaction" : tools.singleton_of(str),
        "use_coulomb" : tools.singleton_of(tools.bool_from_str),
        # Relative observables
        "observable_names" : tools.list_of(str),
        # Calculation
        "hw" : tools.singleton_of(float)
    }
    key_value_dict = tools.extract_key_value_pairs(
        tokenized_lines,conversions
    )

    # legacy support: force interaction to "JISP16" for early runs where
    # interaction field was provided as reserved field but not set to "JISP16"
    if ("interaction" in key_value_dict):
        if (key_value_dict["interaction"] == "RESERVED"):
            key_value_dict["interaction"]="JISP16"

    # provide "coulomb" as preferred field name to match mfdn results analysis
    if ("use_coulomb" in key_value_dict):
        key_value_dict["coulomb"] = key_value_dict["use_coulomb"]
        
    # update to params dictionary
    self.params.update(key_value_dict)


def parse_observables(self,tokenized_lines):
    """Legacy support: Parse ambiguously named "Observables" sections.

    These are to be renamed in future spncci runs!

    """

    #trap unfortunate overload of Observables section name
    if (tokenized_lines[0][0]!="filenames"):
        parse_observable_rmes(self,tokenized_lines)
        return

    conversions = {
        "filenames" : tools.list_of(str)
    }
    key_value_dict = tools.extract_key_value_pairs(
        tokenized_lines,conversions
    )
    self.params["observable_names"] = key_value_dict["filenames"]


def parse_spj_listing(self,tokenized_lines):
    """ Parse ...

    Future: May be adding gex quantum number.

    Warning: Hard-coded gex=0.
    """

    self.spj_listing = np.array(
        tokenized_lines,
        dtype=[("subspace_index",int),("J",float),("dim",int)]
        )
    gex = 0
    self.Jgex_values = list((J,gex) for J in self.spj_listing["J"])

def parse_baby_spncci_listing(self,tokenized_lines):
    """ Parse ...
    """
    table = np.array(
        tokenized_lines,
        dtype=[
            ("subspace_index",int),("irrep_family_index",int),
            ("Nsigmaex",int),("sigma.N",float),("sigma.lambda",int),("sigma.mu",int),
            ("Sp",float),("Sn",float),("S",float),
            ("Nex",int),("omega.N",float),("omega.lambda",int),("omega.mu",int),
            ("gamma_max",int),("upsilon_max",int),("dim",int)
        ]
    )
    self.baby_spncci_listing = table

def parse_decompositions_Nex(self,tokenized_lines):
    """ Parse matrices of RMEs.
    """

    self.decompositions["Nex"] = {}
    tokenized_lines_iterator = iter(tokenized_lines)  # so that we can read through sequentially
    for (J,gex) in self.Jgex_values:

        # skip empty eigenspace
        if ((J,gex) not in self.num_eigenvalues):
            continue

        # read table
        lines = itertools.islice(tokenized_lines_iterator,self.params["Nmax"]+1)
        numbers = [[float(x) for x in row] for row in lines]
        self.decompositions["Nex"][(J,gex)]=np.array(numbers,dtype=float)
        ## print("Numbers:",numbers)
        ## print((J,gex),self.decompositions["Nex"][(J,gex)])

def parse_decompositions_baby_spncci(self,tokenized_lines):
    """ Parse ...
    """

    self.decompositions["BabySpNCCI"] = {}
    baby_spncci_dim = len(self.baby_spncci_listing)
    tokenized_lines_iterator = iter(tokenized_lines)  # so that we can read through sequentially
    for (J,gex) in self.Jgex_values:

        # skip empty eigenspace
        if ((J,gex) not in self.num_eigenvalues):
            continue

        # read table
        lines = itertools.islice(tokenized_lines_iterator,baby_spncci_dim)
        numbers = [[float(x) for x in row] for row in lines]
        self.decompositions["BabySpNCCI"][(J,gex)]=np.array(numbers,dtype=float)
        ##print(self.decompositions["BabySpNCCI"][(J,gex)][:,0])

def parse_energies(self,tokenized_lines):
    """ Parse energies.
    """

    # import energy tabulation
    table = np.array(
        tokenized_lines,
        dtype=[("J",float),("gex",int),("n0",int),("E",float)]
        )

    # sort energies into energies dictionary
    for entry in table:
        (J,gex,n0,E)=entry
        n = n0 + 1  # convert to 1-based spectroscopic numbering of states
        self.energies[(J,gex,n)]=E

def parse_observable_rmes(self,tokenized_lines):
    """Parse ...

    Matrix is canonicalized (assuming RMEs are in group theory
    convention and operator has spherical-harmonic-like conjugation
    properties).

    """

    tokenized_lines_iterator = iter(tokenized_lines)  # so that we can read through sequentially

    observable_matrix_header = next(tokenized_lines_iterator,None)
    while (observable_matrix_header):

        # parse header
        conversions = (int,int,float,int,float,int,int,int)
        (observable_index,sector_index,J_bra,gex_bra,J_ket,gex_ket,rows,cols)=[
            conversion(x)
            for (x,conversion) in zip(observable_matrix_header,conversions)
        ]

        # determine canonicalization
        Jg_pair = ((J_bra,gex_bra),(J_ket,gex_ket))
        (Jg_pair_canonical,flipped,canonicalization_factor) = tools.canonicalize_Jg_pair(
            Jg_pair,tools.RMEConvention.kRose
        )

        # prepare matrix key
        observable_name = self.params["observable_names"][observable_index]

        # read matrix
        lines = itertools.islice(tokenized_lines_iterator,rows)
        numbers = [[float(x) for x in row] for row in lines]
        matrix = np.array(numbers,dtype=float)

        # canonicalize matrix for storage
        if (flipped):
            matrix = (1/canonicalization_factor)*matrix.transpose()

        # store matrix
        observable_dict = self.observables.setdefault(observable_name,dict())
        observable_dict[Jg_pair_canonical] = matrix
        ## print("Key:",key)
        ## print(self.observables[key])

        # attempt to read next header
        observable_matrix_header = next(tokenized_lines_iterator,None)


section_handlers = {
    # PARAMETERS
    "Space" : parse_params,
    "Interaction" : parse_params,
    "Relative observables" : parse_params,
    "Observables" : parse_observables,  # legacy support for early runs
    # BASIS
    "SpJ (listing)" : parse_spj_listing,
    "BabySpNCCI (listing)" : parse_baby_spncci_listing,
    # RESULTS
    "Calculation" : parse_params,
    "Energies" : parse_energies,
    "Decompositions: Nex" : parse_decompositions_Nex,
    "Decompositions: BabySpNCCI" : parse_decompositions_baby_spncci,
    "Observable RMEs" : parse_observable_rmes

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
    tokenized_lines = tools.split_and_prune_lines(res_file_lines)
    sections = tools.extracted_sections(tokenized_lines)

    # split out common sections and subsequent groups of results sections
    def is_results_sentinel_section(section):
        """ Identify mesh point separator "pseudo-section" header.

        (Helper function for res_parser_spncci.)
        """
        (section_name,_) = section
        return (section_name == "RESULTS")

    grouped_sections = tools.split_when(is_results_sentinel_section,sections)
    common_sections = list(next(grouped_sections))
    grouped_results_sections = [list(section_group) for section_group in grouped_sections]

    if (verbose):
        print("Section counts")
        print("  Common sections:",len(common_sections))
        for results_section_group in grouped_results_sections:
            print("  Results sections (by group):",len(results_section_group))

    # generate results objects by mesh point
    mesh_data = []
    if (grouped_results_sections):
        # there are results sections: actual mesh, not counting run
        for results_section_group in grouped_results_sections:
            full_section_group = common_sections + results_section_group
            results = spncci_results_data.SpNCCIResultsData()
            parse_mesh_point(results,full_section_group,section_handlers)
            mesh_data.append(results)
    else:
        # no results sections: counting run
        results = spncci_results_data.SpNCCIResultsData()
        parse_mesh_point(results,common_sections,section_handlers)
        mesh_data.append(results)

    return mesh_data

# register the parser
input.register_data_format('spncci',parser)

if (__name__=="__main__"):
    pass
