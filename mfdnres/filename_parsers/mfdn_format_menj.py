""" mfdn_format_menj.py -- declares descriptor parser

    Language: Python 3
    Mark A. Caprio
    Patrick J. Fasano
    Zhou Zhou
    University of Notre Dame

    01/16/23 (zz): Initiated (based on format_7_ho.py).

"""

import re

# intra-package references
from .. import input

def parser(filename):
    """Parse results filename in format 7, restricted to the ho basis
    special case, but allowing for natural orbitals built on this basis.

    Args:
        filename (string): filename (as basename)

    Returns:
        (dict): info parsed from filename

    """

    regex = re.compile( # TODO: add 3b support
        # prolog
        r"run(?P<run>\w+)"
        r"\-(?P<code_name>((mfdn)|(obscalc-ob)|(transitions-ob)|(transitions-tb))[^\-]*)"
        r"\-(?P<descriptor>"
        # descriptor contents
        r"Z(?P<Z>\d+)\-N(?P<N>\d+)"
        r"\-(?P<MEID>.+)"
        # r"\-(?P<interaction>.+)\-coul(?P<coulomb>\d)"
        r"\-hw(?P<hw>[\d\.]+)"
        r"(\-a_cm(?P<lawson>[\d\.]+))?"
        r"\-Nmax(?P<Nmax>\d+)"
        # r"(\-Ncutob(?P<Ncut>\d+))?"
        # r"(?P<mixed_parity_flag>x)?(?P<fci_flag>\-fci)?"
        # r"(\-Mj(?P<M>-?[\d\.]+))?"
        # r"(\-lan(?P<lanczos>\d+))?"
        # r"(\-tol(?P<tolerance>\d+\.\d+[eE][+-]\d+))?"
        # # natorb group (optional)
        # r"("  # begin natorb group
        #   r"(?P<natural_orbital_flag>\-natorb)"
        #   r"(\-J(?P<natural_orbital_J>[\d\.]+)\-g(?P<natural_orbital_g>[01])\-n(?P<natural_orbital_n>[\d]+))?"
        #   r"(\-no(?P<natural_orbital_iteration>\d+))?"
        # r")?"  # end natorb group
        # # decomposition group (optional)
        # r"("  # begin decomposition group
        #   r"\-J(?P<decomposition_J>[\d\.]+)\-g(?P<decomposition_g>[01])\-n(?P<decomposition_n>[\d]+)"
        #   # task_descriptor_decomposition_1 has "op" prefix before "decomposition_operator_name"
        #   # task_descriptor_decomposition_2 has no prefix before "decomposition_type"
        #   r"\-(op)?(?P<decomposition_type>.+)\-dlan(?P<decomposition_lanczos>\d+)"
        #   r"(?P<decomposition_flag>\-decomp)?"
        # r")?"  # end decomposition group
        # # subset index (optional)
        # r"(\-subset(?P<subset_index>\d+))?"
        # epilog
        # r".*"
        r").(?P<extension>((res)|(out)|(lanczos)))"
    )

    flag_conversions = {
        # "mixed_parity_flag": (lambda s : (s=="x")),
        # "fci_flag": (lambda s : (s=="-fci")),
        # "natural_orbital_flag": (lambda s : (s=="-natorb")),
        # "decomposition_flag": (lambda s : (s=="-decomp")),
    }

    conversions = {
        "Z": int,
        "N": int,
        "MEID": str,
        # "interaction": str,
        # "coulomb": int,
        "hw": float,
        # "lawson": float,
        "Nmax": int,
        # "Ncut": int,
        # "M": float,
        # "lanczos": int,
        # # natorb group (optional)
        # "natural_orbital_J": float, "natural_orbital_g": int, "natural_orbital_n": int,
        # "natural_orbital_iteration": int,
        # # decomposition group (optional)
        # "decomposition_type": str,
        # "decomposition_J": float, "decomposition_g": int, "decomposition_n": int,
        # "decomposition_lanczos": int,
        # # subset index (optional)
        # "subset_index": int,
        }

    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for MFDn results filename: " + filename)
    info = match.groupdict()

    # convert fields
    for key, conversion in flag_conversions.items():
        info[key] = conversion(info[key])
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key]) if (info[key] is not None) else None

    return info

input.register_filename_format("mfdn_format_menj", parser)
