""" mfdn_format_7_ho.py -- declares descriptor parser

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
"""

import re

# intra-package references
from .. import input

def parser(filename):
    """Parse results filename in format 7, restricted to the ho basis
    special case, but allowing for natural orbitals built on this basis.

    Args:
        filename (string) : filename (as basename)

    Returns:
        (dict) : info parsed from filename

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
        r"("   # begin natorb group
          r"(?P<natural_orbital_flag>\-natorb)"
          r"(\-J(?P<J>[\d\.]+)\-g(?P<g>[01])\-n(?P<n>[\d]+))?"
          r"(\-no(?P<natural_orbital_iteration>\d+))?"
        r")?"  # end natorb group
        r"(?P<decomposition_flag>\-decomp)?"
        r"(\-subset(?P<subset_index>\d+))?"
        # epilog
        r").(?P<extension>((res)|(out)|(lanczos)))"
    )

    conversions = {
        "Z" : int,
        "N" : int,
        "interaction" : str,
        "coulomb" : int,
        "hw" : float,
        "lawson" : float,
        "Nmax" : int,
        "Ncut" : int,
        "mixed_parity_flag" : (lambda s  :  (s=="x")),
        "fci_flag" : (lambda s  :  (s=="-fci")),
        "M" : float,
        "lanczos" : int,
        "natural_orbital_flag" : (lambda s  :  (s=="-natorb")),
        "J": float, "g": int, "n": int,
        "decomposition_flag" : (lambda s  :  (s=="-decomp")),
        "natural_orbital_iteration" : int
        }

    match = regex.match(filename)
    if (match == None):
        raise ValueError("bad form for MFDn results filename: " + filename)
    info = match.groupdict()

    # convert fields
    for key in conversions:
        conversion = conversions[key]
        info[key] = conversion(info[key]) if (info[key] is not None) else None

    # build natorb base state
    info["natorb_base_state"] = (info.pop("J"), info.pop("g"), info.pop("n"))

    return info

input.register_filename_format("mfdn_format_7_ho", parser)

if (__name__ == "__main__"):

    filename = r"run0000-mfdn-Z2-N6-Daejeon16-coul1-hw05.000-a_cm20-Nmax02-Mj0.0-lan500-tol1.0e-06-natorb-no0.res"
    info = input.parse_filename(filename, filename_format="mfdn_format_7_ho")
    print(filename)
    print(info)

    filename = r"run0000-mfdn-Z2-N6-Daejeon16-coul1-hw05.000-a_cm20-Nmax02x-Mj0.0-lan500-tol1.0e-06.res"
    info = input.parse_filename(filename, filename_format="mfdn_format_7_ho")
    print(filename)
    print(info)

    filename = r"runpjf0015-mfdn15-Z3-N4-JISP16-coul1-hw20.000-a_cm40-Nmax02-Mj0.5-lan1000-tol1.0e-06-natorb-no0.res"
    info = input.parse_filename(filename, filename_format="mfdn_format_7_ho")
    print(filename)
    print(info)

    filename = r"runpjf0069-mfdn15-Z2-N1-Daejeon16-coul1-hw22.500-a_cm0-Nmax14-Mj0.5-lan400-tol1.0e-06-natorb-J00.5-g0-n01-no0.res"
    info = input.parse_filename(filename, filename_format="mfdn_format_7_ho")
    print(filename)
    print(info)
