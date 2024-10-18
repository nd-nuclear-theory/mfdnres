"""data.py

    NCCI observable tabulation and plotting.

    Mark A. Caprio
    University of Notre Dame

    - 03/16/21 (mac): Created, as ncci_plot.py.
    - 03/18/21 (mac):
        + Implement pandas-based tabulation as wrapper around mfdnres.analysis.make_<observable>_table
        + Implement pandas-based write_hw_scan_plot
        + Implement generic nuclide_observable support (including compount difference/ratio observables)
        + Implement hw scan plotting routines and canned write_hw_scan_plot
    - 04/02/21 (mac): Integrate into mfdnres as mfdnres.data.
    - 04/06/21 (mac): Support positive Nmax_relative in styling functions.
    - 04/08/21 (mac):
        + Support mesh selection and slicing in make_hw_scan_data.
        + Support styling overrides in Nmax_plot_style.
    - 04/13/21 (mac): Add add_expt_marker_band.
    - 04/15/20 (pjf): Add missing element symbols.
    - 04/24/21 (mac): Add unary "minus" compound observable.
    - 05/25/21 (mac): Add Nmax_label_text.  Rename make_qn_text to qn_text.
    - 06/14/21 (mac): Add options to qn_text to control subscript/superscript.
    - 07/14/21 (mac): Support asymmetric error in add_expt_marker_band and add add_data_marker.
    - 09/21/21 (mac): Refactor tabulation of observables to use extensible observable registry.
    - 10/14/21 (mac): Add nuclide and qn tuple manipulation tools (from emratio_obs.py).
    - 10/25/21 (mac/zz): Add "rme" observable and "fix-sign-to" compound observable.
    - 02/13/22 (mac): Extend isotope label formatting and add isotope_str.
    - 02/24/22 (zz): Add E1p,E1n,E1 support in rme and rtp funtions.
    - 03/14/22 (mac): Add break_label_at_symbol.
    - 04/08/22 (mac):
         + Add LevelSelector and associated machinery.
         + Hide "fix-sign-to" in observable descriptor.
    - 04/17/11 (mac): Propose ObservableExtractor interface.
    - 05/01/22 (mac): Add element_symbol labeling function.
    - 05/22/22 (mac): Add element_str labeling function.
    - 06/09/22 (mac/aem): Add LevelSelectorOverride.
    - 06/21/22 (mac): Add nuclide_str() and qn_str().
    - 07/24/22 (mac): Provide traceback diagnostics for failed level selection in resolve_qn().
    - 07/30/22 (mac): Provide support for observable.Observable objects in tabulation/plotting.
    - 07/31/22 (mac): Move LevelSelector out to submodule level.
    - 11/19/22 (mac): Overhaul handling of styling keyword arguments in plotting functions
    - 02/09/23 (mac): Provide observable_labelpad pass-through option to set_up_hw_scan_axes().
    - 02/27/23 (mac): Provide hw_labelpad pass-through option to set_up_hw_scan_axes().
    - 03/19/23 (mac):
      + Add add_hw_scan_plot_Nmax_labels().
      + Add tick specification support to set_up_hw_scan_axes().
      + Improve pass-through option handling in add_observable_panel_label().         
    - 03/21/23 (mac): Add set_up_hw_scan_secondary_axis().
    - 04/26/23 (mac): Support log axis in set_up_hw_scan_axes().
    - 05/24/23 (mac): Provide side option for add_hw_scan_plot_Nmax_labels().
    - 07/08/23 (mac): Provide legend_position option for add_hw_scan_plot_Nmax_labels().
    - 07/09/23 (mac): Support value None for option Nmax_label_tagged_index in
        add_hw_scan_plot_Nmax_labels().
    - 07/24/23 (mac): Simplify option names for add_hw_scan_plot_Nmax_labels().
    - 10/24/23 (mac): Provide observable_scale, tick_specifier, and labelpad options
        for set_up_Nmax_scan_axes().
    - 11/28/23 (mac): 
        + Provide label_text option for add_hw_scan_plot_Nmax_labels().
        + Change add_hw_scan_plot_Nmax_labels() option legend_index default to None.
    - 12/29/23 (mac): Provide zorder and linestyle options for add_expt_marker_band().
    - 09/01/24 (mac): Support (J,g, [n]) quantum numbers in qn_text(). 
    - 09/27/24 (mac):
        + Remove support for legacy "tuple" observables.
        + Add support for generic key tuples in make_hw_scan_data().
    - 10/17/24 (mac): Provide legend_xy option for add_hw_scan_plot_Nmax_labels().
"""

import collections
import os
import traceback

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

from . import (
    analysis,
    level,
    ## observable,
    ticks,
    tools
)

################################################################
# global plot styling
################################################################

# Options for mpl.rcParams:
#
#     >>> mpl.rcParams.update(mfdnres.data.SENSIBLE_PLOT_STYLE)

SENSIBLE_PLOT_STYLE = {
    # font
    "font.family": "serif",
    "mathtext.fontset": "dejavuserif",
    # lines
    "lines.linewidth": 1,
    "axes.prop_cycle": mpl.cycler(color=["black"]),
    # ticks
    "xtick.labelsize": "small",
    "xtick.direction": "in",
    "xtick.major.pad": 1,
    "ytick.labelsize": "small",
    "ytick.direction": "in",
    "ytick.major.pad": 1,
    # legend
    "legend.labelspacing": 0.2,  # compact spacing of entries
    "legend.frameon": False,  # no frame or background
    "legend.fancybox": False,  # no rounded corners (if turn frame back on)
    "legend.framealpha": None,  # no transparency (if turn frame back on)
    "legend.edgecolor": "black",
    "legend.fontsize": "small",
}

################################################################
# range utility
################################################################

# FigGeometry.m
#
# ExtendInterval[PRange:{_?NumericQ,_?NumericQ},PDiff:{_?NumericQ,_?NumericQ},Mode:(Abs|Absolute)]:=
#   PRange+PDiff*{-1,+1};
#
# ExtendInterval[PRange:{_?NumericQ,_?NumericQ},PFrac:{_?NumericQ,_?NumericQ},Mode:Scaled]:=
#   PRange+PFrac*{-1,+1}*-Subtract@@PRange;

def extend_interval(bounds,extensions):
    """ Extend bounds (2-vector) by extensions (scalar or 2-vector).
    """
    b = np.array(bounds)
    e = np.array(extensions)
    return b+e*np.array([-1,+1])

def extend_interval_relative(bounds,extensions):
    """ Extend bounds (2-vector) by relative extensions (scalar or 2-vector).
    """
    b = np.array(bounds)
    e = np.array(extensions)
    return b+e*np.array([-1,+1])*(b[1]-b[0])

################################################################
# partitioning iterable
################################################################

def partitions(iterable, r):
    """Generator function to decompose iterable into tuples of length r.

    The iteration may end with one final shorter tuple containing any residual
    elements.

    >>> list(partitions(range(10),3))

    [(0, 1, 2), (3, 4, 5), (6, 7, 8), (9,)]


    TODO 04/08/21 (mac): Reimplement more elegantly in terms of next()
    calls to iterable.

    NOTE 05/16/21 (mac): Compare note from Python Library help for zip():

    "The left-to-right evaluation order of the iterables is guaranteed. This
    makes possible an idiom for clustering a data series into n-length groups
    using zip(*[iter(s)]*n). This repeats the same iterator n times so that each
    output tuple has the result of n calls to the iterator. This has the effect
    of dividing the input into n-length chunks."

    Arguments:

        iterable (iterable): iterable from which to take items

        r (int): sublist length

    Yields:

        (tuple): group of (up to) r elements

    """
    pool = tuple(iterable)
    while len(pool) > 0:
        yield pool[:r]
        pool = pool[r:]

################################################################
# tools: nuclide and qn tuple manipulation
################################################################

def canonicalized_nuclide(nuclide):
    """ Transpose (Z,N) if needed to make Z<=N. """
    return tuple(sorted(nuclide))

def is_canonical_nuclide(nuclide):
    """ Nuclide (Z,N) with Z<=N. """
    return nuclide == canonicalized_nuclide(nuclide)

def conjugate_nuclide(nuclide):
    """ Transpose (Z,N) to (N,Z). """
    return tuple(reversed(nuclide))

def excited(qn):
    """ Excite (J,g,n) to (J,g,n+1). """
    (J,g,n) = qn
    return (J,g,n+1)

################################################################
# tools: label formatting
################################################################

def break_label_at_symbol(label, symbol, separator = "$\n$"):
    """Insert newline before given character in text label.

    Intended for to, e.g., split ratios at the solidus.

    By default, inserts not only a newline, but also a closing and reopening
    "$", on the assumption that the text is breaking in math mode.

    Example:

        >>> mfdnres.data.break_label_at_symbol(r"$E(4^+)/E(2^+)$", "/")

        '$E(4^+)$\n$/E(2^+)$'

    Arguments:

        label (str): Unbroken label

        symbol (str): Symbol at which to break

        separator (str, optional): String to introduce to generate linebreak

    Returns:

        (str): Broken label

    """

    return label.replace(symbol,separator+symbol)


################################################################
# basic text labels
################################################################

# FigText.m
# ElementAbbreviations=
# 	{"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P",
# 	 "S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu",
# 	 "Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc",
# 	 "Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La",
# 	 "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
# 	 "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",
# 	 "At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf",
# 	 "Es","Fm","Md","No","Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
# 	"Uun","Uuu","Uub"}; (* data circa Mathematica 5.0 *)
ELEMENT_SYMBOLS = [
    "n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P",
    "S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu",
    "Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc",
    "Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La",
    "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
    "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",
    "At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf",
    "Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn",
    "Nh","Fl","Mc","Lv","Ts","Og"
]

## def element_str(Z):
##     """Generate text label component for element symbol.
##
##     Arguments:
##
##         Z (int): Z for element
##
##     Returns:
##
##         label (str): label string, to be interpreted in math mode
##
##     """
##
##     element_symbol = ELEMENT_SYMBOLS[Z] if Z<len(ELEMENT_SYMBOLS) else str(Z)
##     label = r"\mathrm{{{}}}".format(element_symbol)
##     return label

def element_symbol(Z):
    """Generate text label component for element symbol.

    Arguments:

        Z (int): Z for element

    Returns:

        label (str): label string, to be interpreted in math mode

    """

    element_symbol = ELEMENT_SYMBOLS[Z] if Z<len(ELEMENT_SYMBOLS) else str(Z)
    label = r"\mathrm{{{}}}".format(element_symbol)
    return label

def element_str(Z, lower = False):
    """Generate simple string for element, for use in filenames, e.g., "Dy".

    Arguments:

        Z (int): Z for element

        lower (bool, optional): force lowercase

    Returns:

        (str): simple string representation of nuclide

    """
    element_symbol = ELEMENT_SYMBOLS[Z] if Z<len(ELEMENT_SYMBOLS) else str(Z)
    if lower:
        element_symbol = element_symbol.lower()
    return element_symbol

def isotope(nuclide, format = None, as_tuple = False):
    """Generate text label component for nuclide.
    
    Name inspired by eponymous commend from LaTeX isotope package.

    Arguments:

        nuclide (tuple): (Z,N)

        format (str, optional): format code for label ("AS"=A+Symbol,
        "AZSN"=A+Z+Symbol+N, "tuple"=(Z,N))

        as_tuple (bool, optional, deprecated): return (Z,N) label rather than
        standard nuclide symbol (deprecated in favore of using format option)

    Returns:

        label (str): label string, to be interpreted in math mode

    """

    (Z,N) = nuclide
    element_symbol = ELEMENT_SYMBOLS[Z] if Z<len(ELEMENT_SYMBOLS) else str(Z)
    A = sum(nuclide)
    if as_tuple:
        format = "tuple"  # legacy

    if (format == None) or (format == "AS"):
        # standard "^A Symbol"
        label = r"^{{{}}}\mathrm{{{}}}".format(A,element_symbol)
    elif format == "AZSN":
        # "^A _Z Symbol _N"
        label = r"^{{{}}}_{{{}}}\mathrm{{{}}}^{{}}_{{{}}}".format(A,Z,element_symbol,N)  # dummy superscript is for subscript alignment
    elif format == "tuple":
        label = r"({nuclide[0]:d},{nuclide[1]:d})".format(nuclide=nuclide)
    else:
        raise ValueError("unrecognized format option".format(format))
    
    return label


def nuclide_str(nuclide):
    """Generate simple string for nuclide code, for use in filenames, e.g., "Z03-N03".

    Implements special case of isotope_str for format="ZN".

    Arguments:

        nuclide (tuple): (Z,N)

    Returns:

        (str): simple string representation of nuclide

    """
    (Z,N) = nuclide
    label = "Z{nuclide[0]:02d}-N{nuclide[1]:02d}".format(nuclide=nuclide)
    return label

def isotope_str(nuclide, format = None, lower = False):
    """Generate simple string for isotope symbol, for use in filenames, e.g., "156Dy".

    Arguments:

        nuclide (tuple): (Z,N)

        format (str, optional): format code for label ("AS"=A+Symbol,
        "As"=A+symbol, "ZN"=Zxx-Nxx); if None, defaults to "AS"

        lower (bool, optional, deprecated): force lowercase (redundant to format
        option value "As")

    Returns:

        (str): simple string representation of isotope symbol

    """
    (Z,N) = nuclide
    if (format == None) or (format == "AS") or (format == "As"):
        element_symbol = ELEMENT_SYMBOLS[Z] if Z<len(ELEMENT_SYMBOLS) else str(Z)
        if (format == "As") or lower:
            element_symbol = element_symbol.lower()
        A = sum(nuclide)
        label = "{}{}".format(A, element_symbol)
    elif format == "ZN":
        label = nuclide_str(nuclide)
    else:
        raise ValueError("unrecognized format option".format(format))
    return label

ISOTOPE_STR_PATTERN = re.compile(r"^(?P<A>\d+)(?P<sym>[A-Za-z]+)$")

def parse_isotope_str(label):
    """Parse simple string for nuclide, e.g. "156Dy", into (Z,N) tuple.

    This is the inverse of isotope_str(). Note that this function is *almost*
    case insensitive; "n" (neutron) is distinguished from "N" (nitrogen). This
    is the only instance where the case is relevant.

    Arguments:
        label (str): simple string representation of nuclide

    Returns:
        (tuple of int): (Z,N)
    """
    # extract string components
    match = ISOTOPE_STR_PATTERN.match(label)
    if not match:
        raise ValueError("invalid isotope string: {}".format(label))
    A = int(match["A"])
    element_symbol = match["sym"]

    # canonicalize to capitalized symbol if not "neutron"
    if element_symbol != "n":
        element_symbol = element_symbol.capitalize()

    # determine atomic number
    Z = ELEMENT_SYMBOLS.index(element_symbol)
    if Z > A:
        raise ValueError("invalid nuclide: A={:d}, Z={:d}".format(A,Z))

    nuclide = (Z,A-Z)
    return nuclide

def qn_text(qn,show_parity=True,show_index=True):
    """ Generate text label component for quantum numbers.

    Arguments:

        qn (tuple): (J,g) or (J,g,n) quantum numbers

        show_parity (bool, optional): whether or not to show parity (subscript)

        show_index (bool, optional): whether or not to show index (subscript)

    Returns:

        label (str): label string, to be interpreted in math mode
    """

    if len(qn)==3:
        J, g, n = qn
    elif len(qn)==2:
        J, g = qn
        n = None
    else:
        raise ValueError("Unexpected form for quantum numbers: {}".format(qn))

    # am.HalfInt.Str() is missing in Python
    ## J = am.HalfInt(int(2*qn[0]),2)
    twice_J=int(2*J)
    J_str = "{}/2".format(twice_J) if twice_J % 2 else twice_J//2
    if show_parity:
        P_str = "+" if g==0 else "-"
    else:
        P_str = ""
    if show_index and (n is not None):
        n_str = "{:d}".format(n)
    else:
        n_str = ""

    label = r"{{{}}}^{{{}}}_{{{}}}".format(J_str,P_str,n_str)
    return label

def qn_str(qn):
    """Generate simple string for (J,g,n) quantum numbers, for use in filenames, e.g., "00.0-0-1".

    Arguments:

        qn (tuple): (J,g,n) quantum numbers

    Returns:

        (str): simple string representation of qn

    """
    label = "{:04.1f}-{:1d}-{:02d}".format(*qn)
    return label

HW_AXIS_LABEL_TEXT = r"$\hbar\omega~(\mathrm{MeV})$"
NMAX_AXIS_LABEL_TEXT = r"$N_{\mathrm{max}}$"

################################################################
# tools for dispatching level.Level or traditional quantum numbers
################################################################

def resolve_qn(results_data, level_selector, verbose=False):
    """Resolve (J,g,n) quantum number tuple or level selector to quantum number tuple.

    Arguments:

        results_data (mfdnres.ResultsData): Results data

        level_selector (tuple): (J,g,n) tuple or level.Level object

    Returns:

        resolved_qn_list (list): (J,g,n) tuple or None

    """
    if isinstance(level_selector, tuple):
        resolved_qn = level_selector
    elif isinstance(level_selector, level.Level):
        try:
            resolved_qn = level_selector.select_level(results_data)
        except Exception as err:
            print("level_selector.select_level failed with exception: {}".format(err))
            traceback.print_exception(type(err), value=err, tb=err.__traceback__)
            ##traceback.print_tb(err.__traceback__)
            raise
    else:
        raise(TypeError("Unexpected value for level selector"))
    
    return resolved_qn


def resolve_qn_text(level_selector):
    """Resolve LaTeX label text for (J,g,n) quantum number tuple or level selector.

    Arguments:

        results_data (mfdnres.ResultsData): Results data

        level_selector (list): (J,g,n) or level.Level object

    Returns:

        resolved_qn_text (str): LaTeX label text

    """

    if isinstance(level_selector,tuple):
        label = qn_text(level_selector)
    elif isinstance(level_selector,level.Level):
        label = level_selector.label_text
    else:
        raise(TypeError("Unexpected value for level selector"))

    return label


def make_observable_axis_label_text(nuclide_observable):
    """ Generate axis label (with units) for observable, given nuclide_observable.

    Arguments:

        nuclide_observable (tuple): standard nuclide/observable pair or compound

    Returns:

        label (str): label string, to be interpreted in math mode
    """

    observable_str, units_str = nuclide_observable.axis_label_text
        
    if units_str is None:
        label = observable_str
    else:
        label = r"{}~({})".format(observable_str,units_str)

    return label


def make_interaction_text(interaction_coulomb):
    """ Make interaction text, given interaction_coulomb.

    Arguments:

        interaction_coulomb (tuple): interaction/coulomb specifier

    Returns:

        label (str): label string, to be interpreted in math mode
    """
    label = r"\mathrm{{{}}}".format(interaction_coulomb[0])
    return label


def Nmax_label_text(Nmax_highlight,Nmax_max=None):
    """ Generate Nmax label of form Nmax=* or Nmax=*(*).

    Arguments:

        Nmax_highlight (int): last "full scan" Nmax

        Nmax_max (int, optional): last Nmax

    Returns:

        label (str): label string, to be interpreted in math mode

    """
    if Nmax_max == None:
        Nmax_max = Nmax_highlight

    if Nmax_highlight==Nmax_max:
        Nmax_combo_text = "{:d}".format(Nmax_max)
    else:
        Nmax_combo_text = "{:d}({:d})".format(Nmax_highlight,Nmax_max)
    label = "N_{{\mathrm{{max}}}}={}".format(Nmax_combo_text)
    return label

################################################################
# Nmax line styling
################################################################

def Nmax_dashing_emratio(Nmax_relative,base_length=8,exponent_scale=8):
    """Provide dashing pattern based on relative Nmax.

    Inspired by Zhou's dashing:

        (4*r, 4*(1-r))
        r=nmax/max_Nmax  (goes from 1 at highest Nmax, to 0 at floor Nmax0)

    But here we use logarithmic scale to maximum Nmax need not be specified, and
    line for Nmax0 would not vanish.

    Used in, e.g., emratio [Phys. Rev. C 104, 034319 (2021)].

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:

        (tuple): dashing directive

    """
    if Nmax_relative >= 0:
        return (None,None)
    else:  # Nmax_relative < 0
        r = 2**(Nmax_relative/exponent_scale)
        return (base_length*r,base_length*(1-r))

    raise ValueError("invalid Nmax_relative {}".format(Nmax_relative))


def Nmax_dashing_scidraw(Nmax_relative):
    """Provide dashing pattern based on relative Nmax.

    Follows RotationPlots.m "NmaxDashing" dashing convention, based on
    Mathematica code for plots for, e.g., [Eur. Phys. J. A 56, 120 (2020)].

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:

        (tuple): dashing directive

    """

    length_by_Nmax_relative = {
        -2: 6,
        -4: 2,
    }

    if Nmax_relative >= 0:
        return (None,None)
    elif Nmax_relative in length_by_Nmax_relative:
        r = length_by_Nmax_relative[Nmax_relative]
        return (r,r)
    else:  # Nmax_relative < -4:
        return (1,2)

    raise ValueError("invalid Nmax_relative {}".format(Nmax_relative))


def Nmax_dashing_natorb(Nmax_relative):
    """Provide dashing pattern based on relative Nmax.

    Inspired by RotationPlots.m "NmaxDashing" dashing convention, as adapted by
    pjf for plots for natorb [Phys. Rev. C 105, 054301 (2022)].

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:

        (tuple): dashing directive

    """
    
    length_by_Nmax_relative = {
        -2: 10,
        -4: 8,
        -6: 6,
        -8: 4,
        -10: 3,
        -12: 2,
    }

    if Nmax_relative >= 0:
        return (None,None)
    elif Nmax_relative in length_by_Nmax_relative:
        r = length_by_Nmax_relative[Nmax_relative]
        return (r,r)
    else:  # Nmax_relative < -12:
        return (1,1)

    raise ValueError("invalid Nmax_relative {}".format(Nmax_relative))


def Nmax_color(Nmax_relative):
    """Provide color based on relative Nmax.

    Positive Nmax_relative is permitted to allow addition of ad hoc high-Nmax
    data points.

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:

        (str): color directive

    """
    if Nmax_relative >= 0:
        return "firebrick"
    elif Nmax_relative < 0:
        return "black"

    raise ValueError("invalid Nmax_relative {}".format(Nmax_relative))


def Nmax_marker_face_color(Nmax_relative):
    """Provide face color based on relative Nmax.

    Positive Nmax_relative is permitted to allow addition of ad hoc high-Nmax
    data points.

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:

        (str): color directive

    """
    if Nmax_relative > 0:
        return "white"
    else:
        return None

    raise ValueError("invalid Nmax_relative {}".format(Nmax_relative))


def Nmax_symbol_scale(Nmax_relative,exponent_scale=8):
    """ Provide symbol scale based on relative Nmax.

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:

        (float): symbol scale factor (max 1.0)
    """

    return 2**(Nmax_relative/exponent_scale)

    
def Nmax_plot_style(
        Nmax_relative,
        marker_size=6,
        Nmax_symbol_scale=Nmax_symbol_scale,
        Nmax_marker_face_color=None,
        Nmax_dashing=Nmax_dashing_emratio,
        Nmax_color=Nmax_color
):
    """Provide plot styling kwargs based on relative Nmax.

    Styling is meant for curve of fixed Nmax in hw scan plot.

    The various optional arguments control how styling details are determined.

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

        marker_size (float,optional): base marker size for maximal Nmax

        Nmax_symbol_scale, ... (callable, optional): functions to provide
            specific styling parameters as function of relative Nmax

    Returns:

        (dict): kwargs for ax.plot

    """

    return dict(
        markersize=marker_size*Nmax_symbol_scale(Nmax_relative),
        markerfacecolor=(None if Nmax_marker_face_color is None else Nmax_marker_face_color(Nmax_relative)),
        ## linewidth=1,
        dashes=Nmax_dashing(Nmax_relative),
        color=Nmax_color(Nmax_relative),
        )


def hw_plot_style(
        hw,
        ## marker_size=6,
        ## Nmax_symbol_scale=Nmax_symbol_scale,
        ## Nmax_marker_face_color=None,
        ## Nmax_dashing=Nmax_dashing,
        ## Nmax_color=Nmax_color
):
    """Provide plot styling kwargs based on hw.  (WIP)

    Styling is meant for curve of fixed hw in Nmax scan plot.

    The various optional arguments control how styling details are determined.

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

        marker_size (float,optional): base marker size for maximal Nmax

        Nmax_symbol_scale, ... (callable, optional): functions to provide
            specific styling parameters as function of relative Nmax

    Returns:

        (dict): kwargs for ax.plot

    """

    return dict(
        ##markersize=marker_size,
        ##markersize=marker_size*Nmax_symbol_scale(Nmax_relative),
        ##markerfacecolor=(None if Nmax_marker_face_color is None else Nmax_marker_face_color(Nmax_relative)),
        ## linewidth=1,
        ##dashes=Nmax_dashing(Nmax_relative),
        ##color=Nmax_color(Nmax_relative),
        )


################################################################
# (Nmax,hw) multi-indexed data ("hw scan")
################################################################


def hw_scan_descriptor(interaction_coulomb, nuclide_observable, verbose=False):
    """Generate standard descriptor string for a (nuclide,observable) pair.

    Arguments:

        interaction_coulomb (tuple or str): tuple of (interaction, use_coulomb),
        or may simply be given as an interaction string (in which case
        use_coulomb defaults to False)

        nuclide_observable (tuple): standard nuclide/observable pair or compound

    Returns:
        descriptor (str): descriptor string

    """

    if verbose:
        print("Generating hw_scan_descriptor: {} {}".format(interaction_coulomb,nuclide_observable))

    # trap interaction only (no use_coulomb)
    if isinstance(interaction_coulomb, tuple):
        interaction, use_coulomb = interaction_coulomb
    else:
        interaction, use_coulomb = interaction_coulomb, False
        
    observable_descriptor = nuclide_observable.descriptor_str
    
    descriptor="hw-scan_{interaction:s}-{use_coulomb:1d}_{observable_descriptor}".format(
        interaction=interaction, use_coulomb=use_coulomb,
        observable_descriptor=observable_descriptor,
    )

    return descriptor


def hw_scan_drop_nan(observable_data):
    """ Drop nan values from (Nmax,hw) observable data.

    Arguments:
        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

    Returns:
        (pd.DataFrame): data multi-indexed by (Nmax,hw)
    """

    clean_data = observable_data[observable_data["value"].notna()]
    return clean_data


KEY_DESCRIPTOR_NMAX_HW = (("Nmax", int), ("hw", float))

def make_hw_scan_data(
        mesh_data, observable, *,
        selector=None,
        key_descriptor=KEY_DESCRIPTOR_NMAX_HW,
        Nmax_range=None, hw_range=None,
        mesh_ranges=None,
        verbose=False):
    """Tabulate generic observable vs. (Nmax,hw), for scan plots.

    Despite the name, this data tabulation function serves to generate data for
    either hw scans (curves representing fixed Nmax) or Nmax scans (curves
    representing fixed hw).

    Tabulation format (for default key descriptor):

        Nmax hw value

    Arguments:

        mesh_data (list of ResultsData): data set to include

        observable (mfdnres.observable.Observable): observable object

        selector (dict, optional): parameter-value pairs for selection using analysis.selected_mesh_data,
            e.g., {"interaction":interaction,"coulomb":coulomb}

        key_descriptor (tuple of tuple, optional): dtype descriptor for key

        Nmax_range (tuple of float, optional): range to which to limit first
        mesh parameter, which is by default Nmax

        hw_range (tuple of float, optional): range to which to limit second mesh
        parameter, which is by default hw

        mesh_ranges (tuple of tuple of float, optional): range to which to limit
        mesh parameters, which have meanings as specified by key_descriptor

    Returns:
        observable_data (np.array): scan data, with rows (Nmax,hw,value)

    """

    

    ## if not isinstance(observable, observable.Observable):
    ##     raise ValueError("Invalid observable {}".format(observable))
    
    if selector is None:
        selector = {}
    mesh_data_selected = analysis.selected_mesh_data(mesh_data,selector,verbose=verbose)
    observable_data = observable.data(mesh_data_selected, key_descriptor, verbose=verbose)

    # drop nan values
    if verbose:
        print("Observable data before purging NaNs")
        print(observable_data)
    observable_data = hw_scan_drop_nan(observable_data)

    # slice on (Nmax,hw)
    Nmax_slice = slice(None) if Nmax_range is None else slice(*Nmax_range)
    hw_slice = slice(None) if hw_range is None else slice(*hw_range)
    observable_data = observable_data.loc[(Nmax_slice, hw_slice), :]  # pandas MultiIndex slicing
    if mesh_ranges is not None:
        mesh_slices = tuple(map(slice,mesh_ranges))
        observable_data = observable_data.loc[mesh_slices, :]  # pandas MultiIndex slicing

    if verbose:
        print(observable_data)

    return observable_data

def write_hw_scan_data(descriptor,observable_data,directory="data",format_str_obs="9.5f"):
    """
    Write hw scan table for one or more observables.

    Arguments:

        descriptor (str): descriptor string to use in filename

        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)
    """
    # generate table
    table = observable_data.reset_index()
    output_str=table.to_string(
        index=False,header=False,
        formatters={
            "Nmax":(lambda x:format(x,"2d")),
            "hw":(lambda x:format(x,"7.3f")),
        },
        float_format=(lambda x:format(x,format_str_obs)),
    )

    # write table
    output_file_name=os.path.join(directory,"{}_table.dat").format(descriptor)
    with open(output_file_name, 'wt') as out_file:
        out_file.write(output_str)

def set_up_hw_scan_axes(
        ax, nuclide_observable, hw_range, observable_range,
        hw_range_extension=(0.05,0.05), observable_range_extension=(0.05,0.05),
        observable_scale=None,
        observable_axis_label_text=None,
        hw_labelpad=None,
        observable_labelpad=None,
        hw_tick_specifier=None,
        observable_tick_specifier=None,
):
    """ Set up axis ranges, labels, and ticks for hw scan plot.

    Arguments:

        ax (mpl.axes.Axes): axes object

        nuclide_observable (tuple): standard nuclide/observable pair or compound

        hw_range (tuple of float): x range, before extension

        observable_range (tuple of float): y range, or None for matplotlib auto

        observable_scale (str): y scale ("linear" or "log")

        observable_axis_label_text (str, optional): override for observable axis label text

        hw_range_extension (tuple of float, optional): x range relative extension

        observable_range_extension (tuple of float, optional): y range relative extension

        hw_labelpad (scalar, optional): pass-though labelpad option for xlabel

        observable_labelpad (scalar, optional): pass-though labelpad option for ylabel

        hw_tick_specifier (tuple, optional): tick specification (min,max,step,num_subdivision) for hw ticks

        observable_tick_specifier (tuple, optional): tick specification (min,max,step,num_subdivision) for observable ticks

    """

    # set ticks
    #
    # Note that the tick specification must come *before* setting range limits,
    # since the range limits automatically readjust when you set the ticks.
    if hw_tick_specifier is not None:
        x_ticks = ticks.linear_ticks(*hw_tick_specifier)
        ticks.set_ticks(ax,"x",x_ticks)
    if observable_tick_specifier is not None:
        y_ticks = ticks.linear_ticks(*observable_tick_specifier)
        ticks.set_ticks(ax,"y",y_ticks)

    # set limits
    ax.set_xlim(*extend_interval_relative(hw_range,hw_range_extension))
    if observable_scale=="log":
        # Note: Override any range extension for log scale
        observable_range_extension=(0.,0.)
        ax.set_yscale("log")
    if (observable_range is not None) and np.isfinite(observable_range[0]).all():
        ax.set_ylim(*extend_interval_relative(observable_range,observable_range_extension))
        
    # set axis labels
    ax.set_xlabel(HW_AXIS_LABEL_TEXT, labelpad=hw_labelpad)
    if observable_axis_label_text is None:
        observable_axis_label_text = make_observable_axis_label_text(nuclide_observable)
    ax.set_ylabel(
        r"${}$".format(observable_axis_label_text),
        labelpad=observable_labelpad,
    )


def set_up_hw_scan_secondary_axis(
        ax, observable, 
        observable_norm_scale,
        observable_norm_labelpad=None,
        observable_norm_axis_label_text=None,
        observable_norm_tick_specifier=None,
):
    """Set up axis ranges, labels, and ticks for hw scan plot.

    Arguments:

        ax (mpl.axes.Axes): axes object

        observable (observable.Observable): observable object (must provide
            secondary_axis_label_text property)

        observable_norm_scale (float): normalization scale (i.e., denominator of
            ration, such as r^2 or r^4), or None to suppress axis

        observable_norm_labelpad (scalar, optional): pass-though labelpad option for ylabel

        observable_norm_axis_label_text (str, optional): override for observable axis label text

        observable_norm_tick_specifier (tuple, optional): tick specification (min,max,step,num_subdivision) for observable ticks

    Returns:

        ax_secondary_y (mpl.axes.Axes): axes object

    """

    if observable_norm_scale is None:
        return
    
    # create secondary axis
    def scale_functions(secondary_scale):
        """Generate scaling functions for secondary scale.

        DEBUGGING: If secondary_scale is a loop variable, such a wrapping
        function as this one (which rebinds the current value to the local
        argument variable) is needed.  Direct use of lambdas in

            functions=((lambda x: x*secondary_scale), (lambda x: x/secondary_scale))

        leaves variable secondary_scale inside lambda bound to the variable
        secondary_scale outside the lambda, and thus mutable.  All secondary
        scales then are determined by the single final value of the loop
        variable secondary_scale, which is the value in effect when the axes are
        rendered.
        """
        return ((lambda x: x*secondary_scale), (lambda x: x/secondary_scale))
    ax_secondary_y = ax.secondary_yaxis(
        'right',
        functions=scale_functions(observable_norm_scale)
    )
    
    # set ticks
    if observable_norm_tick_specifier is not None:
        y_ticks = ticks.linear_ticks(*observable_norm_tick_specifier)
        ticks.set_ticks(ax_secondary_y,"y",y_ticks)

    # set axis label
    if observable_norm_axis_label_text is None:
        observable_norm_axis_label_text = observable.secondary_axis_label_text
    ax_secondary_y.set_ylabel(
        r"${}$".format(observable_norm_axis_label_text),
        labelpad=observable_norm_labelpad,
    )

    return ax_secondary_y
    
def set_up_Nmax_scan_axes(
        ax, nuclide_observable, Nmax_range, observable_range,
        Nmax_range_extension=(0.05,0.05),
        observable_range_extension=(0.05,0.05),
        observable_scale=None,
        observable_axis_label_text=None,
        Nmax_labelpad=None,
        observable_labelpad=None,
        Nmax_tick_specifier=None,
        observable_tick_specifier=None,
):
    """ Set up axes.

    Arguments:

        ax (mpl.axes.Axes): axes object

        nuclide_observable (tuple): standard nuclide/observable pair or compound

        Nmax_range (tuple of int): x range, before extension

        observable_range (tuple of float): y range, or None for matplotlib auto

        observable_scale (str): y scale ("linear" or "log")

        observable_axis_label_text (str, optional): override for observable axis label text

        Nmax_range_extension (tuple of float, optional): x range relative extension

        observable_range_extension (tuple of float, optional): y range relative extension

        Nmax_labelpad (scalar, optional): pass-though labelpad option for xlabel

        observable_labelpad (scalar, optional): pass-though labelpad option for ylabel

        Nmax_tick_specifier (tuple, optional): tick specification (min,max,step,num_subdivision) for Nmax ticks

        observable_tick_specifier (tuple, optional): tick specification (min,max,step,num_subdivision) for observable ticks

    """

    # set ticks
    #
    # Note that the tick specification must come *before* setting range limits,
    # since the range limits automatically readjust when you set the ticks.
    if Nmax_tick_specifier is not None:
        x_ticks = ticks.linear_ticks(*Nmax_tick_specifier)
        ticks.set_ticks(ax,"x",x_ticks)
    if observable_tick_specifier is not None:
        y_ticks = ticks.linear_ticks(*observable_tick_specifier)
        ticks.set_ticks(ax,"y",y_ticks)

    # set limits
    ax.set_xlim(*extend_interval_relative(Nmax_range,Nmax_range_extension))
    if observable_scale=="log":
        # Note: Override any range extension for log scale
        observable_range_extension=(0.,0.)
        ax.set_yscale("log")
    if (observable_range is not None) and np.isfinite(observable_range[0]).all():
        ax.set_ylim(*extend_interval_relative(observable_range,observable_range_extension))
        
    # set axis labels
    ax.set_xlabel(NMAX_AXIS_LABEL_TEXT, labelpad=Nmax_labelpad)
    if observable_axis_label_text is None:
        observable_axis_label_text = make_observable_axis_label_text(nuclide_observable)
    ax.set_ylabel(
        r"${}$".format(observable_axis_label_text),
        labelpad=observable_labelpad,
    )
        

def add_observable_panel_label(ax,interaction_coulomb,nuclide_observable,**kwargs):
    """ Add observable panel label to plot.

    Standardized label provides: nuclide, observable, interaction

    Arguments:

        ax (mpl.axes.Axes): axes object

        interaction_coulomb (tuple): interaction/coulomb specifier

        nuclide_observable (tuple): standard nuclide/observable pair or compound

        **kwargs: pass-through keyword arguments to ax.annotate

    """

    # panel label
    if isinstance(nuclide_observable, tuple):
        nuclide_text = make_nuclide_text(nuclide_observable)
        observable_text = make_observable_text(nuclide_observable)
    else:
        ##if isinstance(nuclide_observable, observable.Observable):
        nuclide_text = nuclide_observable.nuclide_label_text
        observable_text = nuclide_observable.observable_label_text

    interaction_text = make_interaction_text(interaction_coulomb)

    # combine styling options (last takes precedence)
    kw_defaults=dict(
        xy=(0.95,0.05),
        xycoords="axes fraction",
        multialignment="left",
        horizontalalignment="right",
        verticalalignment="bottom",
        bbox=dict(boxstyle="round",facecolor="white"),
    )
    
    kw_full = {
        **kw_defaults,
        **kwargs,
    }
    
    ax.annotate(
        "${}$ ${}$\n$^{}$".format(nuclide_text,observable_text,interaction_text),
        **kw_full
    )

def add_hw_scan_plot(
        ax,observable_data,Nmax_max,
        Nmax_plot_style=Nmax_plot_style,
        Nmax_plot_style_kw={},
        **kwargs,
):
    """Add hw scan plot to axes.

    Arguments:

        ax (mpl.axes.Axes): axes object

        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

        Nmax_max (int): highest Nmax for styling purposes

        Nmax_plot_style (callable, optional): function to provide styling kwargs
            as function of Nmax_relative (rarely needed, since normally instead
            can simply use Nmax_plot_style_kw with the default styling function
            Nmax_plot_style)

        Nmax_plot_style_kw (dict, optional): styling options to pass through to
        plot styling function defined in Nmax_plot_style argument

        kwargs (Line2D properties, optional): kwargs are used to specify line
        properties not otherwise fixed by the prior arguments

    Returns:

       Nmax_groups (pd.DataFrameGroupBy): curve data grouped by Nmax (for
           possible use in subsequent calls to labeling functions)

    """

    kw_defaults = {
        "markersize": 6,
        "marker": ".",
    }

    Nmax_groups = observable_data.reset_index().groupby("Nmax")
    for Nmax, group in Nmax_groups:

        # combine styling options (last takes precedence)
        kw_full = {
            **kw_defaults,
            **Nmax_plot_style(Nmax-Nmax_max,**Nmax_plot_style_kw),
            **kwargs,
        }
        
        ax.plot(
            group["hw"],group["value"],
            **kw_full,
        )

    return Nmax_groups

def add_hw_scan_plot_Nmax_labels(
        ax, Nmax_groups, label_list,
        legend_index=None,
        side="right",
        data_point_index=None,
        text_displacement=None,
        legend_position="bottom",
        legend_xy=None,
        label_text=None,
):
    """Add Nmax curve labels to previously drawn hw scan plot.

    Arguments:

        ax (mpl.axes.Axes): axes object

        Nmax_groups (pd.DataFrameGroupBy): curve data grouped by Nmax (as
            returned by add_hw_scan_plot())

        label_list (list of int): list of Nmax values for labels [formerly
        Nmax_label_list]

        legend_index (int, optional): index within label_list for Nmax
        label to which to attach the legend "Nmax" (or None to omit legend)
        [formerly Nmax_label_tagged_index]

        side (str, optional): side of curve for label "left" or "right"

        data_point_index (int, optional): index of data point within curve for
        labeling (0 for "left" end of curve, -1 for "right" end of curve); or
        None for default based on side

        text_displacement (tuple, optional): xy displacement in points of text relative to curve point; or
        None for default based on side

        legend_position (str, optional): position of Nmax
        legend label relative to Nmax labels ("bottom" or "top")

        legend_xy (tuple of float, optional): explicit position of Nmax legend
        label relative to the Nmax label to which it is attached, for manual
        fine-tuning (defaults legend_xy=(1,0) for legend_position="bottom",
        legend_xy=(1,1) for or legend_position="top")

        label_text (str, optional): template string for label (applied as format
        string to the value of Nmax), default "{}"; may be used to provide an arbitrary text label for hw
        scan curves

    """

    # TODO (mac, 03/19/23): add in generalizations for call-out lines

    if side=="right":
        if data_point_index is None:
            data_point_index=-1
        if text_displacement is None:
            text_displacement=(+12,+0)
    elif side=="left":
        if data_point_index is None:
            data_point_index=0
        if text_displacement is None:
            text_displacement=(-2,+0)
        
    for Nmax, group in Nmax_groups:

        # extract curve endpoint
        curve_points = group[["hw", "value"]].to_numpy()
        endpoint = curve_points[data_point_index]
        
        # generate Nmax label
        if Nmax in label_list:
            if label_text is None:
                label_text = r"{}"
            resolved_label_text = label_text.format(Nmax)
            Nmax_label = ax.annotate(
                r"${}$".format(resolved_label_text),
                ##r"${}$".format(Nmax),
                xy=endpoint, xycoords="data",
                xytext=text_displacement, textcoords="offset points",
                fontsize="x-small",
                horizontalalignment="right", verticalalignment="center",
                ##arrowprops=dict(arrowstyle="-", linewidth=0.5, shrinkA=1, shrinkB=3),
                ##bbox=dict(boxstyle="square", visible=False, pad=0.),  # to clip call-out line under text
            )

        # add "Nmax" legend label
        if (
                len(label_list)>0
                and legend_index is not None
                and Nmax == label_list[legend_index]
        ):
            if legend_position=="bottom":
                if legend_xy is None:
                    legend_xy = (1,0)
                verticalalignment = "top"
            elif legend_position=="top":
                if legend_xy is None:
                    legend_xy = (1,1)
                verticalalignment = "bottom"
            ax.annotate(
                r"$N_{\mathrm{max}}$",
                xy=legend_xy, xycoords=Nmax_label,
                fontsize="x-small",
                horizontalalignment="right", verticalalignment=verticalalignment,
            )

            
def add_Nmax_scan_plot(
        ax,observable_data,
        hw_plot_style=hw_plot_style,
        hw_plot_style_kw={},
        verbose=False,
        **kwargs,
):
    """Add Nmax scan plot to axes.

    Arguments:

        ax (mpl.axes.Axes): axes object

        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

        hw_plot_style (callable, optional): function to provide styling kwargs
            as function of hw

        hw_plot_style_kw (dict, optional): styling options to pass through to
        plot styling function defined in hw_plot_style argument

        kwargs (Line2D properties, optional): kwargs are used to specify plot
        properties (e.g., marker, markersize) not otherwise fixed by the prior arguments

    """

    kw_defaults = {
        "markersize": 6,
        "marker": ".",
    }
    
    for hw, group in observable_data.reset_index().groupby("hw"):
        if verbose:
            print("hw {}\n {}".format(hw,group))

        # combine styling options (last takes precedence)
        kw_full = {
            **kw_defaults,
            **hw_plot_style(hw,**hw_plot_style_kw),
            **kwargs,
        }
        ## print(kw_defaults,hw_plot_style(hw,**hw_plot_style_kw),kwargs,kw_full)
        
        ax.plot(
            group["Nmax"],group["value"],
            **kw_full,
        )

def write_hw_scan_plot(
        descriptor,
        interaction_coulomb,nuclide_observable,
        observable_data,
        hw_range,observable_range,Nmax_max,
        hw_range_extension=(0.02,0.02),
        observable_range_extension=(0.02,0.02),
        figsize=(6,4),
        directory=".",
        panel_label_kwargs = {},
        verbose=False
):
    """Generate full "canned" hw scan plot.

    This is for "quick and dirty" plots for analysis, and to serve as a model
    for the basic elements of creating a figure containing an hw scan plot.  For
    "production" figures for publication, where you need more control, and may
    combine several plots in different panels, you would generally write control
    code along these lines yourself, instead of using this canned routine.

    Arguments:

        descriptor (str): descriptor string to use in filename

        interaction_coulomb (tuple): interaction/coulomb specifier

        nuclide_observable (tuple): standard nuclide/observable pair or compound

        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

        hw_range (tuple): horizontal axis range

        observable_range (tuple): vertical axis range

        Nmax_max (int): highest Nmax assumed for styling purposes

        hw_range_extension (tuple, optional): fractional range extensions for horizontal axis

        observable_range_extension (tuple, optional): fractional range extensions for vertical axis

        figsize (tuple, optional): mpl figure size argument

        directory (str, optional): output directory

        panel_label_kwargs (dict, optional): keyword arguments to pass through to add_observable_panel_label()

    """
    # initialize plot
    fig, ax = plt.subplots(figsize=figsize)

    # provide axis labeling
    set_up_hw_scan_axes(
        ax,
        nuclide_observable,
        hw_range,
        observable_range,
        hw_range_extension=hw_range_extension,
        observable_range_extension=observable_range_extension,
    )

    # make panel label
    add_observable_panel_label(ax,interaction_coulomb,nuclide_observable,**panel_label_kwargs)

    # make Nmax label
    ax.annotate(
        "$N_{{\mathrm{{max}}}}={}$".format(Nmax_max),
        xy=(0.95,0.97),
        xycoords="axes fraction",
        multialignment="left",
        horizontalalignment="right",
        verticalalignment="top",
        fontsize="x-small",
    )
    
    # generate plot
    add_hw_scan_plot(ax,observable_data,Nmax_max)

    # finalize plot
    figure_file_name = os.path.join(
        directory,
        "{}_plot.pdf".format(descriptor)
        )
    if verbose:
        print(figure_file_name)
    plt.savefig(figure_file_name)
    plt.close()

################################################################
# plot annotations
################################################################

def add_expt_marker_band(
        ax,x_range, y_with_error,
        error_full_rectangle=False,
        color="black", linewidth=1, linestyle="solid",
        error_facecolor="lightgray",
        error_edgecolor="black", error_linewidth=0.5, error_linestyle="solid",
        zorder=None,
):
    """Add marker indicating value with error band (rectangle) and central value (line).

    Limitation: Giving a value for error_linestyle other than "solid" seems to
    break the rendering of the whole axis.

    Arguments:

        ax (mpl.axes.Axes): axes object

        x_range (tuple of float): (x1,x2) range

        y_with_error (float or tuple): (y,dy) given as y, (y,None), (y,dy), or
           (y,(dy_plus,dy_minus)), where all errors should have *positive*
           values (in keeping with the conventions of Axes.errorbar); if None,
           no marker is drawn

        error_full_rectangle (bool, optional): draw full rectangle for error band, instead of just top and bottom lines

        color, linewidth (optional): styling parameters for central value

        error_facecolor, error_edgecolor, error_linewidth, error_linestyle (optional): styling parameters error band

        zorder (optional): styling parameters for whole object


    """

    (x0,x1) = x_range
    if y_with_error is None:
        return
    if type(y_with_error)==tuple:
        (y,y_error) = y_with_error
    else:
        y = y_with_error
        y_error = None
    if type(y_error)==tuple:
        dy_plus, dy_minus = y_error
    else:
        dy_plus = y_error
        dy_minus = y_error

    # support zorder overrides
    if zorder is None:
        zorder_options = {}
    else:
        zorder_options = dict(zorder=zorder)
        
    # error band
    if y_error is not None:
        y0 = y-dy_minus
        y1 = y+dy_plus
        if error_full_rectangle:
            ax.fill(
                [x0,x1,x1,x0],
                [y0,y0,y1,y1],
                edgecolor=error_edgecolor, linewidth=error_linewidth, linestyle=error_linestyle,
                facecolor=error_facecolor,
                **zorder_options,
            )
        else:
            ax.fill(
                [x0,x1,x1,x0],
                [y0,y0,y1,y1],
                edgecolor=error_edgecolor, linewidth=0, linestyle=error_linestyle,
                facecolor=error_facecolor,
                **zorder_options,
            )
            ax.hlines(
                [y0,y1], *x_range,
                color=error_edgecolor, linewidth=error_linewidth,
                **zorder_options,
            )

    # central value
    ax.hlines(
        y, *x_range,
        color=color, linewidth=linewidth, linestyle=linestyle,
        **zorder_options,
    )

def add_data_marker(ax,x,y_with_error,errorbar_kw=dict()):
    """Add marker indicating value with error band (rectangle) and central value (line).

    Arguments:

        ax (mpl.axes.Axes): axes object

        x (float): x coordinate

        y_with_error (float or tuple): (y,dy) given as y, (y,None), (y,dy), or
           (y,(dy_plus,dy_minus)), where all errors should have *positive*
           values (in keeping with the conventions of Axes.errorbar); if None,
           no marker is drawn

        errorbar_kw (dict, optional): options to Axes.errorbar

    """

    if y_with_error is None:
        return
    if type(y_with_error)==tuple:
        (y,y_error) = y_with_error
    else:
        y = y_with_error
        y_error = None
    if type(y_error)==tuple:
        dy_plus, dy_minus = y_error
    else:
        dy_plus = y_error
        dy_minus = y_error
    if y_error is None:
        yerr = None
    else:
        ## yerr = [y_error]  # for symmetric error bars
        yerr = [[dy_plus],[dy_minus]]
    ax.errorbar([x],[y],yerr=yerr,**errorbar_kw)
