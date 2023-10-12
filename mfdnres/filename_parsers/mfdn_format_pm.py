"""mfdn_format_pm.py -- declares descriptor parser for Pieter Maris runs

    See mfdn_format_pm_test.py for parsing tests.

PM typical format

    MFDn.res.Z4.N5.JISP16.Nmin1.Nm13.hw20.0.La500.St06.tol1e-6

Caveat: Nmin is only included if nonero.  Nasty double use of period as
decimal and as field separator.

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    07/26/15 (mac): Initiated (based on format_5_ho.py).
    04/27/18 (mac): Rename parameter Mj to M.
    10/12/23 (mac): Provide support for LENPIC runs.

"""

import re

# intra-package references
from .. import input

def parser(filename):
    """ Parses results filename in Pieter's typical format, defined for ho basis.

    Note that no information on Coulomb is provided.

    Args:
        filename (string) : filename (as basename)

    Returns:
        (dict) : info parsed from filename

    """

    regex = re.compile(
        ## r"\-Mj(?P<M>[\d\.]+)"  # TODO Mj
        # prolog
        r"MFDn\.res"
        r"\.(?P<descriptor>"
        # descriptor contents
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
        "interaction" : str,
        "hw" : float,
        "Nmax" : int,
        ##"Mj" : float,
        "lanczos" : int,
        "states" : int,
        "tolerance" : float,
        }

    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for MFDn results filename: " + filename)
    info = match.groupdict()
    ## print(info)

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key])

    return info

input.register_filename_format("mfdn_format_pm", parser)

def parser_lenpic(filename):
    """ Parses results filename in Pieter's typical format, defined for ho basis, for LENPIC runs.

    Args:
        filename (string) : filename (as basename)

    Returns:
        (dict) : info parsed from filename

    """

    regex = re.compile(
        ## r"\-Mj(?P<M>[\d\.]+)"  # TODO Mj
        # prolog
        r"MFDn\.res"
        r"\.(?P<descriptor>"
        # descriptor contents
        r"A(?P<A>\d+)\.Z(?P<Z>\d+)\.N(?P<N>\d+)"
        r"\.Nm(?P<Nmax>\d+)"
        r"\.(?P<interaction>[^\.]+)"
        r"_HO(?P<hw>[\d]+)"
        # epilog
        r")"
    )

    conversions = {
        "A" : int,
        "Z" : int,
        "N" : int,
        "interaction" : str,
        "hw" : float,
        }

    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for MFDn results filename: " + filename)
    info = match.groupdict()
    ## print(info)

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key])

    return info

input.register_filename_format("mfdn_format_pm_lenpic", parser_lenpic)
