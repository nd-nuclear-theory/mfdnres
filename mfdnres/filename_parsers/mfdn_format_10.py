""" mfdn_format_10.py -- declares descriptor parser for traditional shell model runs

    Language: Python 3
    Mark A. Caprio
    Patrick J. Fasano
    University of Notre Dame

    08/01/24 (pjf): Initiated (based on format_7_ho.py).

"""

import re

# intra-package references
from .. import input

def parser(filename):
    """Parse results filename in format 10, for basic shell-model runs.

    Args:
        filename (string): filename (as basename)

    Returns:
        (dict): info parsed from filename

    """

    regex = re.compile(
        # prolog
        r"run(?P<run>\w+)"
        r"\-(?P<code_name>((mfdn)|(obscalc-ob)|(transitions-ob)|(transitions-tb))[^\-]*)"
        r"\-(?P<descriptor>"
        # descriptor contents
        r"Z(?P<Z>\d+)\-N(?P<N>\d+)"
        r"\-(?P<interaction>.+)"
        r"(\-Mj(?P<M>-?[\d\.]+))?"
        r"(\-lan(?P<lanczos>\d+))?"
        r"(\-tol(?P<tolerance>\d+\.\d+[eE][+-]\d+))?"
        # decomposition group (optional)
        r"("  # begin decomposition group
          r"\-J(?P<decomposition_J>[\d\.]+)\-g(?P<decomposition_g>[01])\-n(?P<decomposition_n>[\d]+)"
          # task_descriptor_decomposition_2 has no prefix before "decomposition_type"
          r"\-(?P<decomposition_type>.+)\-dlan(?P<decomposition_lanczos>\d+)"
          r"(?P<decomposition_flag>\-decomp)?"
        r")?"  # end decomposition group
        # subset index (optional)
        r"(\-subset(?P<subset_index>\d+))?"
        # epilog
        r").(?P<extension>((res)|(out)|(lanczos)))"
    )

    conversions = {
        "Z": int,
        "N": int,
        "interaction": str,
        "M": float,
        "lanczos": int,
        # decomposition group (optional)
        "decomposition_type": str,
        "decomposition_J": float, "decomposition_g": int, "decomposition_n": int,
        "decomposition_lanczos": int,
        # subset index (optional)
        "subset_index": int,
        }

    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for MFDn results filename: " + filename)
    info = match.groupdict()

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key]) if (info[key] is not None) else None

    # build decomposition state
    info["decomposition_state"] = (info.pop("decomposition_J"), info.pop("decomposition_g"), info.pop("decomposition_n"))

    return info

input.register_filename_format("mfdn_format_10", parser)
