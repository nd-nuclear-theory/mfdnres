""" spncci.py -- declares descriptor parser

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    05/26/20 (mac): Initiated (based on mfdn_format_7_ho.py).

"""

import re

# intra-package references
from .. import input

def parser(filename):
    """Parse results filename in format 7, restricted to the ho basis
    special case, but allowing for natural orbitals built on this basis.

    Args:
        filename (string) : filename (as basename)

    Returns:
        (dict) : info parsed from filename

    """

    regex = re.compile(
        # prolog
        r"run(?P<run>\w+)"
        ##r"\-(?P<code_name>((mfdn)|(obscalc-ob))[^\-]*)"
        r"\-(?P<descriptor>"
        # descriptor contents
        r"Z(?P<Z>\d+)\-N(?P<N>\d+)"
        r"\-(?P<interaction>.+)\-(?P<coulomb>\d)"
        r"\-(?P<truncation_descriptor>.+)"
        ## r"\-Nmax(?P<Nmax>\d+)"
        # epilog
        r").res"
    )

    conversions = {
        "Z" : int,
        "N" : int,
        "interaction" : str,
        "coulomb" : int,
        }

    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for spncci results filename: " + filename)
    info = match.groupdict()

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key]) if (info[key] is not None) else None

    return info

input.register_filename_format("spncci", parser)

if (__name__ == "__main__"):

    # >>> import mfdnres
    # runaem0092-Z4-N3-Daejeon16-1-LGI0301-Nmax20.res
    # {'run': 'aem0092', 'descriptor': 'Z4-N3-Daejeon16-1-LGI0301-Nmax20', 'Z': 4, 'N': 3, 'interaction': 'Daejeon16', 'coulomb': 1, 'truncation_descriptor': 'LGI0301-Nmax20', 'filename': 'runaem0092-Z4-N3-Daejeon16-1-LGI0301-Nmax20.res', 'nuclide': (4, 3)}

    filename = r"runaem0092-Z4-N3-Daejeon16-1-LGI0301-Nmax20.res"
    info = input.parse_filename(filename, filename_format="spncci")
    print(filename)
    print(info)


    
