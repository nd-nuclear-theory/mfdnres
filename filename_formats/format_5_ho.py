""" format_5_ho.py -- declares descriptor parser

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    6/2/15 (mac): Initiated (as mfdn_descriptor.py).
    6/5/15 (mac): Restructure as subpackage.
    7/3/15 (mac): Rewrite as parser function.  Add fields fci_flag and Mj.
    Last modified 7/3/15.
    
"""

import re

# intra-package references
import mfdnres.descriptor

def parser(filename):
    """ Parses results filename in format 5, restricted to the ho basis special case.

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
        r"\-Nmax(?P<Nmax>\d+)(?P<mixed_parity_flag>[x]?)(?P<fci_flag>[-fci]?)"
        r"\-MM(?P<MM>\d+)"
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
        "MM" : int,
        "lanczos" : int
        }


    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for MFDn results filename")
    info = match.groupdict()

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key])

    # replace MM with Mj
    info["Mj"] = info["MM"]/2
    del info["MM"]

    return info

mfdnres.descriptor.register_filename_format("format_5_ho",parser)

if (__name__ == "__main__"):

    # run0352-mfdn-Z4-N5-JISP16-1-hw20.000-aL100-Nmax10-MM1-lan1000.res

    filename = r"run0352-mfdn-Z4-N5-JISP16-1-hw20.000-aL100-Nmax10-MM1-lan1000.res"
    parser(filename)
