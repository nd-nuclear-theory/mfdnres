""" mfdn_format_7_ho.py -- declares descriptor parser

    See mfdn_format_7_ho_test.py for parsing tests.

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
    12/17/19 (mac): Support format_7_trans (i.e., make diagonalization-related fields optional).
    06/12/20 (mac): Support optional subset index from format_7_trans.
    09/17/20 (mac): Support code names for postprocessor output ("transitions-{ob,tb}").
    12/01/20 (pjf):
        + Support natural orbital base state information in descriptor.
        + Support "decomp" flag.
    12/06/20 (pjf): Add additional decomposition descriptor parsing support.
    10/12/23 (mac): Update decomposition descriptor parsing to support task_descriptor_decomposition_2.

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

    regex = re.compile(
        # prolog
        r"run(?P<run>\w+)"
        r"\-(?P<code_name>((mfdn)|(obscalc-ob)|(transitions-ob)|(transitions-tb))[^\-]*)"
        r"\-(?P<descriptor>"
        # descriptor contents
        r"Z(?P<Z>\d+)\-N(?P<N>\d+)"
        r"\-(?P<interaction>.+)\-coul(?P<coulomb>\d)"
        r"\-hw(?P<hw>[\d\.]+)"
        r"(\-a_cm(?P<lawson>[\d\.]+))?"
        r"\-Nmax(?P<Nmax>\d+)"
        r"(\-Ncutob(?P<Ncut>\d+))?"
        r"(?P<mixed_parity_flag>x)?(?P<fci_flag>\-fci)?"
        r"(\-Mj(?P<M>-?[\d\.]+))?"
        r"(\-lan(?P<lanczos>\d+))?"
        r"(\-tol(?P<tolerance>\d+\.\d+[eE][+-]\d+))?"
        # natorb group (optional)
        r"("  # begin natorb group
          r"(?P<natural_orbital_flag>\-natorb)"
          r"(\-J(?P<natural_orbital_J>[\d\.]+)\-g(?P<natural_orbital_g>[01])\-n(?P<natural_orbital_n>[\d]+))?"
          r"(\-no(?P<natural_orbital_iteration>\d+))?"
        r")?"  # end natorb group
        # decomposition group (optional)
        r"("  # begin decomposition group
          r"\-J(?P<decomposition_J>[\d\.]+)\-g(?P<decomposition_g>[01])\-n(?P<decomposition_n>[\d]+)"
          # task_descriptor_decomposition_1 has "op" prefix before "decomposition_operator_name"
          # task_descriptor_decomposition_2 has no prefix before "decomposition_type"
          r"\-(op)?(?P<decomposition_type>.+)\-dlan(?P<decomposition_lanczos>\d+)"
          r"(?P<decomposition_flag>\-decomp)?"
        r")?"  # end decomposition group
        # subset index (optional)
        r"(\-subset(?P<subset_index>\d+))?"
        # epilog
        r").(?P<extension>((res)|(out)|(lanczos)))"
    )

    flag_conversions = {
        "mixed_parity_flag": (lambda s : (s=="x")),
        "fci_flag": (lambda s  :  (s=="-fci")),
        "natural_orbital_flag": (lambda s  :  (s=="-natorb")),
        "decomposition_flag": (lambda s  :  (s=="-decomp")),
    }

    conversions = {
        "Z": int,
        "N": int,
        "interaction": str,
        "coulomb": int,
        "hw": float,
        "lawson": float,
        "Nmax": int,
        "Ncut": int,
        "M": float,
        "lanczos": int,
        # natorb group (optional)
        "natural_orbital_J": float, "natural_orbital_g": int, "natural_orbital_n": int,
        "natural_orbital_iteration": int,
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
    for key, conversion in flag_conversions.items():
        info[key] = conversion(info[key])
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key]) if (info[key] is not None) else None

    # set default natorb iteration
    if info.get("natural_orbital_iteration") is None:
        info["natural_orbital_iteration"] = 0

    # build natorb base state
    info["natorb_base_state"] = (info.pop("natural_orbital_J"), info.pop("natural_orbital_g"), info.pop("natural_orbital_n"))

    # build decomposition state
    info["decomposition_state"] = (info.pop("decomposition_J"), info.pop("decomposition_g"), info.pop("decomposition_n"))

    # provide legacy decomposition_operator field
    info["decomposition_operator"] = info.get("decomposition_type")
    
    return info

input.register_filename_format("mfdn_format_7_ho", parser)
