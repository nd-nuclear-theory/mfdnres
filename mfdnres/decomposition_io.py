"""decomposition_io.py -- Input and output for decomposition data file

Language: Python 3
Mark A. Caprio
University of Notre Dame

    - 08/24/26 (mac): Created, refactoring from parse_decomp.py and eigenvalues2decomp.py.
"""

import numpy as np

import mfdnres.tools

################################################################
# int/float casting helper function
################################################################

def int_or_float(x):
    """ Cast value to float, then to int if it is equal to its integer part.
    """

    x = float(x)
    if x==int(x):
        x = int(x)
    return x
    

################################################################
# key-value I/O helper functions
################################################################

def extract_variables(tokenized_lines, conversion_type):
    """ Parse tokenized lines as key-value pairs with values of given type.

    Arguments:

        tokenized_lines (list[list[str]]): Tokenized input lines.

        conversion_type (callable): Type casting function for values.
    """
    
    data = dict()
    for tokenized_line in tokenized_lines:

        # validate line format
        valid_line = (len(tokenized_line)==3) and (tokenized_line[1]=="=")
        if (not valid_line):
            raise ValueError("expected key-value line but found {}".format(tokenized_line))

        # extract line parts
        key = tokenized_line[0]
        value_string = tokenized_line[2]
        value = conversion_type(value_string)

        data[key] = value
        
    return data


def generate_key_value_list(key_value_dict):
    """Generate key-value entries from dictionary representation.

    Inspired by ncci.utils.write_namelist().

    Arguments:
        key_value_dict (dict): Dictionary of int or float values by keyword.
    """
    formatters = {
        int:   (lambda n: "{:d}".format(n)),
        float: (lambda x: "{:e}".format(x)),
    }

    def format_val(x):
        if type(x) not in formatters.keys():
            raise exception.ScriptError(
                "{} of type {} cannot be written".format(x, type(x))
                )
        return formatters[type(x)](x)

    lines = []

    # loop over contents
    for key, val in key_value_dict.items():
        # sanity check
        if type(key) is not str:
            raise exception.ScriptError("invalid key: {}".format(key))

        lines.append("{key:s} = {value}".format(key=key, value=format_val(val)))

    return lines


################################################################
# decomposition file parsing
################################################################

def parse_coefs(decomp_data, tokenized_lines):
    """ Parse coefficient key-value pairs.
    """

    coef_by_operator = extract_variables(tokenized_lines, float)
    decomp_data["coefficients"] = coef_by_operator


def parse_statistics(decomp_data, tokenized_lines):
    """ Parse statistics key-value pairs.
    """

    statistics = extract_variables(tokenized_lines, int)
    decomp_data["statistics"] = statistics


def parse_labels(decomp_data, tokenized_lines):
    """ Parse label identifiers.
    """

    label_identifiers = []
    for tokenized_line in tokenized_lines:

        # validate line format
        valid_line = len(tokenized_line)==1
        if (not valid_line):
            raise ValueError("expected label identifier line but found {}".format(tokenized_line))

        # store identifier
        label_identifiers.append(tokenized_line[0])
        
    decomp_data["labels"] = label_identifiers


def parse_eigenvalues(decomp_data, tokenized_lines):
    """ Parse eigenvalue tabulation.
    """

    eigenvalue_by_labels = dict()
    for tokenized_line in tokenized_lines:
        labels = tuple(map(int_or_float, tokenized_line[:-1]))
        value = float(tokenized_line[-1])
        eigenvalue_by_labels[labels] = value
    
    decomp_data["eigenvalues"] = eigenvalue_by_labels


def parse_lanczos(decomp_data, tokenized_lines):
    """ Parse lanczos alpha-beta tabulation.
    """

    rejoined_lines = [" ".join(tokenized_line) for tokenized_line in tokenized_lines]
    alpha_beta_table = np.loadtxt(rejoined_lines, ndmin=2)
    
    decomp_data["lanczos"] = alpha_beta_table
    
    
def parse_decomp_file(decomp_filename):
    """Write decomposition file lines to output stream.

    Arguments:

        decomp_filename (str): Filename.

     Returns:

        decomp_data (dict): Dictionary containing data by section
        ("coefficients", "statistics", "labels", "eigenvalues", "lanczos").

            "coefficients" (dict[str, float]): Coefficient by operator identifier string.

            "statistics" (dict[str, float]): Statistics on eivenvalues and labels.

            "labels" (list[str]): Identifiers for label "quantum numbers".

            "eigenvalues" (dict[tuple, float]): Eigenvalue by symmetry label tuple to eigenvalue.

            "lanczos" (ndarray): Array (two-column) of Lanczos alpha and beta coefficients.

    """

    decomp_data = {}

    # perform high-level parsing into sections
    decomp_file = open(decomp_filename, "r")
    lines = [line for line in decomp_file]
    tokenized_lines = mfdnres.tools.split_and_prune_lines(lines)
    sections = mfdnres.tools.extracted_sections(tokenized_lines)
    decomp_file.close()

    # warn of empty file
    if len(lines)==0:
        print("WARNING: file {} is empty!".format(decomp_file.name))

    # parse sections
    section_handlers = {
        "Coefficients": parse_coefs,
        "Statistics": parse_statistics,
        "Labels": parse_labels,
        "Eigenvalues": parse_eigenvalues,
        "Lanczos": parse_lanczos,
    }
    for section_name, tokenized_lines in sections:
        if section_name in section_handlers:
            try:
                section_handlers[section_name](decomp_data, tokenized_lines)
            except Exception as err:
                print("ERROR: Unexpected content in decomp file section '{}'".format(section_name))
                raise err
    return decomp_data

    
################################################################
# decomposition file generation
################################################################

def generate_decomp_file(decomp_data, header_comment_lines=[]):
    """Write decomposition file lines to output stream.

    Generates "statistics" data from "eigenvalues" data, if "statistics" not
    provided.

    Arguments:

        decomp_data (dict): Dictionary containing data by section
        ("coefficients", ["statistics"], "labels", "eigenvalues", ["lanczos"]).
        See docstring for parse_decomp_file.

        header_comment_lines (list[str]): Text for lines to include in header
        comment (sans leading hashes).

     Returns:

        lines (list[str]): Output lines for decomp file.

    """

    lines = []

    # generate header comment
    for line in header_comment_lines:
        lines.append("# {}".format(line))
    lines.append("")

    # generate coefficient section
    lines.append("[Coefficients]")
    for operator, coef in decomp_data["coefficients"].items():
        line = "{} = {:+e}".format(operator, coef)
        lines.append(line)
    lines.append("")

    # generate statistics section
    lines.append("[Statistics]")
    unique_labels = set()
    unique_eigenvalues = set()
    statistics = decomp_data.get("statistics")
    if statistics is None:
        for labels, eigenvalue in decomp_data["eigenvalues"].items():
            unique_labels.add(labels)
            unique_eigenvalues.add(eigenvalue)
        statistics = {
            "num_entries": len(unique_labels),
            "num_eigenvalues": len(unique_eigenvalues),
        }
    lines += generate_key_value_list(statistics)
    lines.append("")

    # generate labels section
    lines.append("[Labels]")
    lines += decomp_data["labels"]
    lines.append("")

    # generate eigenvalues section
    lines.append("[Eigenvalues]")
    ## lines.append(" ".join(label_identifiers))
    label_format_str_by_type = {
        int: "{:6d} ",
        float: "{:8.1f} ",
    }
    for labels, eigenvalue in decomp_data["eigenvalues"].items():
        line = ""
        for label in labels:
            line += label_format_str_by_type[type(label)].format(label)
        line += "   {:+e}".format(eigenvalue)
        lines.append(line)
    lines.append("")

    # generate lanczos section (optional)
    if "lanczos" in decomp_data:
        lines.append("[Lanczos]")
        for row in decomp_data["lanczos"]:
            lines.append("{:+e} {:+e}".format(*row))
        lines.append("")
    
    return lines

