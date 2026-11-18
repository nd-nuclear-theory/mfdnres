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
    11/18/26 (seb): Add support for num_segments and Nmax_orb parameters.

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
        # Nmax_orb (optional)
        r"(?:\-Nmax_orb(?P<Nmax_orb>\d+))?"
        # num_segments (optional)
        r"(?:\-n(?P<num_segments>\d+))?"
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
        "Nmax_orb": int,
        "num_segments": int,
        }

    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for MFDn results filename: " + filename)
    info = match.groupdict()

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key]) if (info[key] is not None) else None

    return info


if (__name__ == "__main__"):

    files = [
            r"runpjf0040-mfdn15-Z06-N08-Nmax00-Mj0.0.res",
            r"runmfdncounting02-mfdn15-Z03-N03-Nmax04-Mj1.0-Nmax_orb05-n001.res",
            ]
    for filename in files: 
        info = parser(filename)
        print(filename)
        print(info)
else:
    # intra-package references
    from .. import input
    input.register_filename_format("mfdn_format_c1", parser)
