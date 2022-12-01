""" mfdn_format_c1.py -- declares descriptor parser

    Language: Python 3
    Mark A. Caprio
    Patrick J. Fasano
    University of Notre Dame

    02/06/17 (pjf): Initiated (based on format_6_ho.py).
    10/10/17 (mac): Generalize filename format to allow for variable code name
        ("mfdn", "mfdn15", etc.).
    04/27/18 (mac): Rename parameter Mj to M.
    09/06/18 (pjf): Allow negative M.
    09/06/18 (pjf): Allow hyphens in interaction name.

"""

import re

def parser(filename):
    """Parse results filename in format c1.

    NOTE: Currently restricted to the Nmax truncation special case.

    Args:
        filename (string) : filename (as basename)

    Returns:
        (dict) : info parsed from filename

    """

    regex = re.compile(
        # prolog
        r"run(?P<run>\w+)"
        r"\-(?P<code_name>mfdn[^\-]*)"
        r"\-(?P<descriptor>"
        # descriptor contents
        r"Z(?P<Z>\d+)\-N(?P<N>\d+)"
        r"\-Nmax(?P<Nmax>\d+)"
        r"(?P<mixed_parity_flag>x)?(?P<fci_flag>\-fci)?"
        r"\-Mj(?P<M>-?[\d\.]+)"
        # epilog
        r").res"
    )

    conversions = {
        "Z" : int,
        "N" : int,
        "Nmax" : int,
        "mixed_parity_flag" : (lambda s  :  (s=="x")),
        "fci_flag" : (lambda s  :  (s=="-fci")),
        "M" : float,
        }

    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for MFDn results filename: " + filename)
    info = match.groupdict()

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key])

    return info


if (__name__ == "__main__"):

    filename = r"runpjf0040-mfdn15-Z06-N08-Nmax00-Mj0.0.res"
    info = parser(filename)
    print(filename)
    print(info)
else:
    # intra-package references
    from .. import input
    input.register_filename_format("mfdn_format_c1", parser)
