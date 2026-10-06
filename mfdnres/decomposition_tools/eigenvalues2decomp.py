"""eigenvalues2decomp.py

    Generate decomposition data files from legacy eigenvalues and coefs files.

    Processes all eigenvalues files in current working directory.

    Mark A. Caprio
    University of Notre Dame

    + 08/21/26 (mac): Created.
    + 08/24/26 (mac): Refactor output routines to mfdnres.decomposition_io.
    + 10/05/26 (mac): Take directory names as arguments.
"""

import argparse
import glob
import re
import os.path

import numpy as np

import mfdnres.decomposition
import mfdnres.decomposition_io
   
    
################################################################
# parse arguments
################################################################

def parse_args():
    """Parse arguments.

    Returns:
        (argparse.Namespace) parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Generate decomposition data files from legacy eigenvalues and coefs files.",
        usage="%(prog)s eigenvalues_dir decomp_dir\n",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

    # positional arguments
    parser.add_argument("eigenvalues_dir", help="Legacy eigenvalue and coef file directory")
    parser.add_argument("decomp_dir", help="Decomposition data file directory")

    # options
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress output text")
    
    return parser.parse_args()


################################################################
# process coefficients
################################################################

def obtain_coefs(decomposition_type, coefs_filename):
    """ Read coefficients from coefs file, or else generate unit coefficient.

    Arguments:

        decomposition_type (str): Identifier for standard decomposition type.

        coefs_filename (str): Filename for coefs file.

    Returns:

        coef_by_operator (dict[str, float]): Dict mapping operator identifier to coefficient.
    """

    operators_by_decomposition_type = {
        "LS": ["S2", "L2"],
        "SU3": ["CSU3"],
        "Sp3R": ["CSp3R"],
        "U3S": ["Nex", "CSU3", "S2"],
        "U3LS": ["Nex", "CSU3", "S2", "L2"],
        "U3SpSnS": ["Nex", "CSU3", "Sp2", "Sn2", "S2"],
        "U3LSpSnS": ["Nex", "CSU3", "Sp2", "Sn2", "S2","L2"],
        "Sp3RS": ["CSp3R", "S2"],
        "Sp3RSpSnS": ["CSp3R", "Sp2", "Sn2", "S2"],
    }
    operators = operators_by_decomposition_type[decomposition_type]

    if len(operators)==1:
        # handle case of single operator (no coefs file)
        coef_by_operator = {operators[0]: 1.}
    else:
        coefs = np.loadtxt(coefs_filename)
        coef_by_operator = {
            operator: coef
            for operator, coef in zip(operators, coefs)
        }
        
    return coef_by_operator
    
    
################################################################
# process eigenvalues
################################################################

def obtain_eigenvalues(decomposition_type, eigenvalue_filename):
    """ Read eigenvalues from eigenvalues file, and coarse-grain labels.

    Arguments:

        decomposition_type (str): Identifier for standard decomposition type.

        eigenvalue_filename (str): Filename for eigenvalue file.

    Returns:

        label_identifiers (list[str]): Identifiers for label "quantum numbers".

        eigenvalue_by_labels (dict[tuple, float]): Dict mapping symmetry label tuple to eigenvalue.
    """

    # define labeling
    labels_type = mfdnres.decomposition.LABEL_CLASS_BY_DECOMPOSITION_TYPE[decomposition_type]
    label_identifiers = labels_type._fields
    labels_subsetting_function = mfdnres.decomposition.labels_subsetting_function(labels_type)
    
    # read data
    table = np.loadtxt(eigenvalue_filename, ndmin=2)  # need to force 2-dim array in case there is just one row
    label_length = table.shape[1]-1  # all but last column constitutes label
    if label_length == 6:
        source_labels_type = mfdnres.decomposition.SOURCE_LABEL_CLASS_BY_DECOMPOSITION_TYPE_SANS_L[decomposition_type]
    elif label_length == 7:
        source_labels_type = mfdnres.decomposition.SOURCE_LABEL_CLASS_BY_DECOMPOSITION_TYPE_WITH_L[decomposition_type]

    #  process data
    eigenvalue_by_labels = {}
    ##print(table)
    for row in table:

        # extract source labels
        ##print(row)
        source_labels = row[:-1]
        source_labels = source_labels_type(*map(mfdnres.decomposition_io.int_or_float, source_labels))

        # deduce target labels
        labels = labels_subsetting_function(source_labels)

        # save eigenvalue
        #
        # A given coarse grained label may occur more than once.  Only add
        # labels to dict if not already present (preserves order of first
        # occurrence).  And flag if eigenvalue is inconsistent with earlier
        # occurrence.
        eigenvalue = row[-1]
        if eigenvalue != eigenvalue_by_labels.setdefault(labels, eigenvalue):
            print("WARN: Inconsistent eigenvalue for labels {}".labels)
    
    return label_identifiers, eigenvalue_by_labels


################################################################
# testing
################################################################

def generate_decomp_file_test():
    """ Generate example decomp file data.
    """
    header_info = ["generated by eigenvalues2decomp"]
    coef_by_operator = {
        "Nex": 6.529787991265456526e+01,
        "CSU3": 2.389423454890835075e+00,
        "S2": 1.000296887033979253e-02,
        }
    label_identifiers = ["N_omega", "lambda_omega", "mu_omega", "S"]
    eigenvalue_by_labels = {
        (0, 0, 0, 0.0): +0.000000e+00,
        (0, 0, 0, 1.0): +2.000594e-02,
        }

    lines = mfdnres.decomposition_io.generate_decomp_file(header_info, coef_by_operator, label_identifiers, eigenvalue_by_labels)
    output_str = "\n".join(lines) + "\n"
    print(output_str)

################################################################
# main
################################################################

def main():

    args = parse_args()
    eigenvalues_dir = args.eigenvalues_dir
    decomp_dir = args.decomp_dir

    # create target directory for decomp files
    os.makedirs(decomp_dir, exist_ok=True)

    # process eigenvalues (and coefs) files
    eigenvalues_filenames = glob.glob(os.path.join(eigenvalues_dir, "*_eigenvalues.dat"))
    for eigenvalues_filename in eigenvalues_filenames:

        # parse eigenvalue filename
        eigenvalues_basename = os.path.basename(eigenvalues_filename)

        regex = re.compile(
            r"(?P<base>(decomposition_Z(?P<Z>\d+)_N(?P<N>\d+)_Nmax(?P<Nmax>\d+)_(?P<decomposition_type>.+)))_eigenvalues.dat"
        )
        match = regex.match(eigenvalues_basename)
        if (match == None):
            raise ValueError("bad form for eigenvalue filename: {}".format(eigenvalues_basename))
        info = match.groupdict()
        base = info["base"]
        if not args.quiet:
            print(base)

        # determine decomposition type
        decomposition_type_overrides = {
            # map special identifiers appearing in eigenvalue filenames to
            # corresponding decomposition type
            "CSU3": "SU3",
            "CSp3R": "Sp3R",
        }
        decomposition_type = info["decomposition_type"]
        if decomposition_type in decomposition_type_overrides:
            decomposition_type = decomposition_type_overrides[decomposition_type]

        # generate derived filenames
        coefs_filename = os.path.join(eigenvalues_dir, "{}_coefs.dat".format(base))
        decomp_filename_template = "Z{Z:s}-N{N:s}-Nmax{Nmax:s}-{decomposition_type:s}.decomp"
        decomp_filename = os.path.join(
            decomp_dir,
            decomp_filename_template.format(
                Z=info["Z"],
                N=info["N"],
                Nmax=info["Nmax"],
                decomposition_type=decomposition_type,  # override decomposition_type from input filename
            ),
        )

        # generate header
        header_comment_lines = [
            "generated by eigenvalues2decomp",
            "source: {}".format(base),
            "decomposition type: {}".format(decomposition_type),
        ]

        # obtain coefs
        coef_by_operator = obtain_coefs(decomposition_type, coefs_filename)

        # obtain eigenvalues
        label_identifiers, eigenvalue_by_labels = obtain_eigenvalues(decomposition_type, eigenvalues_filename)

        # write decomp file
        decomposition_data = {
            "coefficients": coef_by_operator,
            "labels": label_identifiers,
            "eigenvalues": eigenvalue_by_labels,
            }
        lines = mfdnres.decomposition_io.generate_decomp_file(decomposition_data, header_comment_lines)
        output_str = "\n".join(lines) + "\n"
        data_file = open(decomp_filename, "w")
        data_file.write(output_str)
        data_file.close()


if (__name__ == "__main__"):
    
    main()
        
