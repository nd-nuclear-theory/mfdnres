""" mfdn_format_6_ho.py -- declares descriptor parser

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    07/03/15 (mac): Initiated (based on format_5_ho.py).
    04/27/18 (mac): Rename parameter Mj to M.

"""

import re

# intra-package references
from .. import descriptor

def parser(filename):
    """ Parses results filename in format 6, restricted to the ho basis special case.

    Args:
        filename (string) : filename (as basename)

    Returns:
        (dict) : info parsed from filename

    """

    regex = re.compile(
        # prolog
        r"run(?P<run>\w+)"
        r"\-mfdn"
        r"\-(?P<descriptor>"
        # descriptor conents
        r"Z(?P<Z>\d+)\-N(?P<N>\d+)"
        r"\-(?P<interaction>[^\-]+)\-(?P<coulomb>\d)"
        r"\-hw(?P<hw>[\d\.]+)"
        r"\-aL(?P<lawson>[\d\.]+)"
        r"\-Nmax(?P<Nmax>\d+)(?P<mixed_parity_flag>x?)(?P<fci_flag>\-fci)?"
        r"\-Mj(?P<M>[\d\.]+)"
        r"\-lan(?P<lanczos>\d+)"
        # epilog
        r").res"
    )

    conversions = {
        "Z" : int,
        "N" : int,
        "coulomb" : int,
        "hw" : float,
        "lawson" : float,
        "Nmax" : int,
        "mixed_parity_flag" : (lambda s  :  (s=="x")),
        "fci_flag" : (lambda s  :  (s=="-fci")),
        "M" : float,
        "lanczos" : int
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


descriptor.register_filename_format("mfdn_format_6_ho", parser)

if (__name__ == "__main__"):

    filename = r"run0376-mfdn-Z4-N6-JISP16-1-hw15.000-aL100-Nmax04-Mj0.0-lan1500.res"
    info = descriptor.parse_res_filename(filename, filename_format="mfdn_format_6_ho")
    print(filename)
    print(info)
