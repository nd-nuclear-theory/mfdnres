"""descriptor.py

    Defines filename (or task descriptor) parsing for MFDn res files.


    A filename parsing function is assumed to provide the following
    mandatory fields:

        "run" (str): run name (may be null)
        "descriptor" (str): the part of the filename which
             describes the run parameters
        "Z", "N" (int): proton and neutron numbers

    E.g., under "format5ho", the filename

        "run0364-mfdn-Z4-N3-JISP16-1-hw20.000-aL100-Nmax10-MM1-lan1000.res"

    yields

        "run" : "0364"
        "descriptor" : "Z4-N3-JISP16-1-hw20.000-aL100-Nmax10-MM1-lan1000"
        "Z" : 4
        "N" : 3
        "interaction" : "JISP16"
        ...

    The wrapper parse_res_filename will add the field "nuclide" as a
    tuple of int, e.g.,

        "nuclide" : (4,3)

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    06/02/15 (mac): Initiated (as mfdn_descriptor.py).
    06/05/15 (mac): Restructure as part of package.
    09/06/18 (pjf): Include filename in returned info.

"""

import os
import re

################################################################
# parser registry
################################################################

# global registration variables
filename_format_parser = {}

def register_filename_format(format_name,parser):
    """Register information for parsing filename.

    Args:
        format_name (str): name for filename format
        parser (callable): function for parsing filename

    """

    filename_format_parser[format_name] = parser


################################################################
# wrapper function
################################################################

def parse_res_filename(filename,filename_format):
    """Parses mfdn results filename.

    Only the basename is considered, extracted via os.path.basename,
    i.e., any preceding path is ignored.

    Args:
        filename (str): filename to parse

    Returns: (dict) : dictionary with keys for parameters ("run",
        "descriptor", "Z", "N", ...) parsed from filename, plus
        "nuclide" as a tuple of int

    """

    # parse filename
    basename = os.path.basename(filename)
    parser = filename_format_parser[filename_format]
    info = parser(basename)

    # define nuclide tuple
    info["filename"] = filename
    info["nuclide"] = (info["Z"],info["N"])

    return info


if (__name__ == "__main__"):
    pass
