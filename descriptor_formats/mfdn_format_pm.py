"""mfdn_format_pm.py -- declares descriptor parser

PM typical format

    MFDn.res.Z4.N5.JISP16.Nmin1.Nm13.hw20.0.La500.St06.tol1e-6

Caveat: Nmin is only included if nonero.  Nasty double use of period as
decimal and as field separator.

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    7/26/15 (mac): Initiated (based on format_5_ho.py).
    Last modified 7/26/15.

"""

import re

# intra-package references
if (__name__ == "__main__"):
    import sys
    import os
    sys.path.append(os.path.join(sys.path[0],"..",".."))
import mfdnres.descriptor

def parser(filename):
    """ Parses results filename in Pieter's typical format, defined for ho basis.

    Note that no information on Coulomb is provided.

    Args:
        filename (string) : filename (as basename)

    Returns:
        (dict) : info parsed from filename

    """

    regex = re.compile(
        ## r"\-Mj(?P<Mj>[\d\.]+)"  # TODO Mj
        # prolog
        r"MFDn\.res"
        r"\.(?P<descriptor>"
        # descriptor conents
        r"Z(?P<Z>\d+)\.N(?P<N>\d+)"
        r"\.(?P<interaction>[^\.]+)"
        r"(\.Nmin1)?"  # ignore Nmin, though could be used in detecting mixed parity runs
        r"\.Nm(?P<Nmax>\d+)"
        r"\.hw(?P<hw>[\d]+.[\d]+)"
        r"\.La(?P<lanczos>\d+)"
        r"\.St(?P<states>\d+)"
        r"\.tol(?P<tolerance>[\de\+\-]+)"
        # epilog
        r")"
    )

    conversions = {
        "Z" : int,
        "N" : int,
        "hw" : float,
        "Nmax" : int,
        ##"Mj" : float,
        "lanczos" : int,
        "states" : int,
        "tolerance" : float,
        }

    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for MFDn results filename")
    info = match.groupdict()
    ## print(info)

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key])

    return info

mfdnres.descriptor.register_filename_format("mfdn_format_pm",parser)

if (__name__ == "__main__"):

    filename = r"MFDn.res.Z4.N5.JISP16.Nmin1.Nm13.hw20.0.La500.St06.tol1e-6"
    info = mfdnres.descriptor.parse_res_filename(filename,filename_format="mfdn_format_pm")
    print(filename)
    print(info)
