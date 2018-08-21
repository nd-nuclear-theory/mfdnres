""" mfdn_format_8.py -- declares descriptor parser

    Language: Python 3
    Patrick J. Fasano
    University of Notre Dame

    3/18/18 (pjf): Initiated (based on mfdn_format_7_ho.py).

"""

import re

# intra-package references
from .. import descriptor

_parity_map = {"g0": +1, "g1": -1, "gx": 0}


def _truncation_parser(substr):
    regex = re.compile(
        r""
        )


def parser(filename):
    """Parse results filename in format 8.

    Args:
        filename (string) : filename (as basename)

    Returns:
        (dict) : info parsed from filename

    """

    regex = re.compile(
        # prolog
        r"run(?P<run>\w+)"
        r"\-(?P<code_name>[^\-]+)"
        r"\-(?P<descriptor>"
        # descriptor contents
        r"Z(?P<Z>\d+)\-N(?P<N>\d+)"
        r"\-(?P<interaction>[^\-]+)\-coul(?P<coulomb>\d)"
        r"\-hw(?P<hw>[\d\.]+)"
        r"\-a_cm(?P<lawson>[\d\.]+)"
        r"\-an(?P<n_coeff>\d+\.\d{3})"
        r"\-bl(?P<l_coeff>\d+\.\d{3})"
        r"\-spWTmax(?P<sp_weight_max>\d+\.\d{3})"
        r"\-((?P<fci_flag>FCI)|WTmax(?P<mb_weight_max>\d+\.\d{3}))"
        r"\-(?P<parity_indicator>g.)"
        r"\-Mj(?P<Mj>[\d\.]+)"
        r"\-its(?P<max_iterations>\d+)"
        r"\-tol(?P<tolerance>\d+\.\d+[eE][+-]\d+)"
        r"((?P<natural_orbital_flag>\-natorb)\-no(?P<natural_orbital_iteration>\d+))?"
        # epilog
        r").res"
    )

    conversions = {
        "Z": int,
        "N": int,
        "interaction" : str,
        "coulomb": int,
        "hw": float,
        "lawson": float,
        "n_coeff": float,
        "l_coeff": float,
        "sp_weight_max": float,
        "mb_weight_max": (lambda s: float(s) if (s is not None) else None),
        "fci_flag": (lambda s: (s == "FCI")),
        "parity_indicator": (lambda s: _parity_map[s]),
        "Mj": float,
        "max_iterations": int,
        "natural_orbital_flag": (lambda s: (s == "-natorb")),
        "natural_orbital_iteration": (lambda i: int(i) if (i is not None) else None)
        }

    match = regex.match(filename)
    if (match is None):
        raise ValueError("bad form for MFDn results filename: " + filename)
    info = match.groupdict()

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key])

    return info


descriptor.register_filename_format("mfdn_format_8", parser)

if (__name__ == "__main__"):

    filename = r"run0000-mfdn15-Z2-N6-Daejeon16-coul1-hw10.000-a_cm0-an01.500-bl01.000-spWTmax12.000-WTmax15.000-g0-Mj0.0-its200-tol1.0e-06.res"
    info = descriptor.parse_res_filename(filename, filename_format="mfdn_format_8")
    print(filename)
    print(info)
