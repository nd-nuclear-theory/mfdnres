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

"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import (
    analysis
)

################################################################
# global plot styling
################################################################

## def init_plt_rc():
##     """ Set global plot styling parameters.
##     """
##     plt.rc('font', family='serif')
##     plt.rc('mathtext', fontset='dejavuserif')
##     plt.rc('xtick', labelsize='small')
##     plt.rc('ytick', labelsize='small')

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
    b = np.array(bounds)
    e = np.array(extensions)
    return b+e*np.array([-1,+1])

def extend_interval_relative(bounds,extensions):
    b = np.array(bounds)
    e = np.array(extensions)
    return b+e*np.array([-1,+1])*(b[1]-b[0])

################################################################
# observable parsing
################################################################

def unpack_observable(observable):
    """Unpack standard observable tuple into type, operator, and qn list.

    Arguments

        observable (tuple): standard observable specifier tuple

    Returns:

        observable_type (str): "energy", ...
    
        observable_operator (str): "M1", ..., or None for energy

        observable_qn_list (list of tuple): [(J1,g1,n1),...]

    """
    observable_type = observable[0]
    if observable_type in {"energy","isospin"}:
        observable_operator = None
        observable_qn_list = observable[1:]
    else:
        observable_operator = observable[1]
        observable_qn_list = observable[2:]

    return (observable_type,observable_operator,observable_qn_list)

################################################################
# descriptor string construction
################################################################

def nuclide_observable_descriptor(nuclide_observable):
    """ Generate standard descriptor string for a (nuclide,observable) pair.

    Arguments:
        nuclide_observable (tuple): standard nuclide/observable pair or compound

    Returns:
        descriptor (str): descriptor string
    """
    # trap compound observable
    if nuclide_observable[0] in {"diff","ratio"}:
        (arithmetic_operation,nuclide_observable1,nuclide_observable2) = nuclide_observable
        return r"{}_{}_{}".format(  # "{}-{}-{}"
            arithmetic_operation,
            nuclide_observable_descriptor(nuclide_observable1),
            nuclide_observable_descriptor(nuclide_observable2)
        )

    # unpack arguments
    (nuclide,observable) = nuclide_observable
    (observable_type,observable_operator,observable_qn_list) = unpack_observable(observable)

    qn_list_str = "-".join([
        "{:04.1f}-{:1d}-{:02d}".format(*observable_qn)
        for observable_qn in observable_qn_list
    ])

    if observable_operator is None:
        descriptor_template = "Z{nuclide[0]:02d}-N{nuclide[1]:02d}-{observable_type}-{qn_list_str}"
    else:
        descriptor_template = "Z{nuclide[0]:02d}-N{nuclide[1]:02d}-{observable_type}-{observable_operator}-{qn_list_str}"
        
    descriptor=descriptor_template.format(
            nuclide=nuclide,
            observable_type=observable_type,
            observable_operator=observable_operator,
            observable=observable,
            qn_list_str=qn_list_str
        )

    return descriptor
    
################################################################
# text labels
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
    "Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt",
    "Uun","Uuu","Uub"
]

def isotope_symbol(nuclide,as_tuple=False):
    """Generate text label component for nuclide.

    Arguments:

        nuclide (tuple): (Z,N)

        as_tuple (bool, optional): return (Z,N) label rather than standard nuclide symbol

    Returns:
    
        label (str): label string, to be interpreted in math mode

    """

    (Z,N) = nuclide
    element_symbol = ELEMENT_SYMBOLS[Z] if Z<len(ELEMENT_SYMBOLS) else str(Z)
    A = sum(nuclide)
    if as_tuple:
        label = r"({nuclide[0]:d},{nuclide[1]:d})".format(nuclide=nuclide)
    else:
        label = r"^{{{}}}\mathrm{{{}}}".format(A,element_symbol)
    return label
    
def make_nuclide_text(nuclide_observable,as_tuple=False):
    """Generate text label component for nuclide.

    Arguments:

        nuclide_observable (tuple): standard nuclide/observable pair or compound

        as_tuple (bool, optional): return (Z,N) label rather than standard nuclide symbol

    Returns:
    
        label (str): label string, to be interpreted in math mode

    """
    # trap compound observable
    if nuclide_observable[0] in {"diff","ratio"}:
        (arithmetic_operation,nuclide_observable1,nuclide_observable2) = nuclide_observable
        same_nuclide = (nuclide_observable1[0]==nuclide_observable2[0])
        if same_nuclide:
            return make_nuclide_text(nuclide_observable1)
        else:
            return r"{}/{}".format(
                make_nuclide_text(nuclide_observable1,as_tuple=as_tuple),
                make_nuclide_text(nuclide_observable2,as_tuple=as_tuple)
            )

    (nuclide,observable) = nuclide_observable

    return isotope_symbol(nuclide,as_tuple=as_tuple)

def make_qn_text(qn):
    """ Generate text label component for quantum numbers.

    Arguments:

        qn (tuple): (J,g,n) quantum numbers

    Returns:
    
        label (str): label string, to be interpreted in math mode
    """

    # am.HalfInt.Str() is missing in Python
    ## J = am.HalfInt(int(2*qn[0]),2)
    twice_J=int(2*qn[0])
    J_str = "{}/2".format(twice_J) if twice_J % 2 else twice_J//2
    P_str = "+" if qn[1]==0 else "-"
    n = qn[2]
    label = r"{}^{}_{}".format(J_str,P_str,n)
    return label

def make_observable_text(nuclide_observable):
    """ Generate text label for observable.

    Arguments:

        nuclide_observable (tuple): standard nuclide/observable pair or compound

    Returns:
    
        label (str): label string, to be interpreted in math mode
    """
   
    # trap compound observable
    if nuclide_observable[0] in {"diff","ratio"}:
        (arithmetic_operation,nuclide_observable1,nuclide_observable2) = nuclide_observable
        if arithmetic_operation == "diff":
            arithmetic_symbol = "-"
        elif arithmetic_operation == "ratio":
            arithmetic_symbol = "/"
        return r"{}{}{}".format(
            make_observable_text(nuclide_observable1),
            arithmetic_symbol,
            make_observable_text(nuclide_observable2)
        )
    
    # unpack arguments
    (nuclide,observable) = nuclide_observable
    (observable_type,observable_operator,observable_qn_list) = unpack_observable(observable)

    # construct label
    if observable_type == "energy":
        observable_str = r"E"
        qn_str = make_qn_text(observable_qn_list[0])
        label = r"{}({})".format(observable_str,qn_str)
    elif observable_type == "isospin":
        observable_str = r"\bar{T}"
        qn_str = make_qn_text(observable_qn_list[0])
        label = r"{}({})".format(observable_str,qn_str)
    elif observable_type == "radius":
        pass  # TODO
    elif observable_type in {"moment","momentsqr"}:
        if observable_operator == "M1":
            observable_str = r"\mu"
        elif observable_operator in {"Dlp","Dln","Dsp","Dsn","Dl0","Dl1","Ds0","Ds1"}:
            observable_str = r"\mu_{{{}}}".format(observable_operator[1:])
        elif observable_operator in {"E2p","E2n","E20","E21"}:
            observable_str = r"eQ_{{{}}}".format(observable_operator[2:])
        qn_str = make_qn_text(observable_qn_list[0])
        label = r"{}({})".format(observable_str,qn_str)
        if observable_type == "momentsqr":
            label = r"[{}]^2".format(label)
    elif observable_type == "rtp":
        if observable_operator == "M1":
            observable_str = r"M1"
        elif observable_operator in {"Dlp","Dln","Dsp","Dsn","Dl0","Dl1","Ds0","Ds1"}:
            observable_str = r"M1_{{{}}}".format(observable_operator[1:])
        elif observable_operator in {"E2p","E2n","E20","E21"}:
            observable_str = r"E2_{{{}}}".format(observable_operator[2:])
        qn_str_1 = make_qn_text(observable_qn_list[0])
        qn_str_2 = make_qn_text(observable_qn_list[1])
        label = r"B({};{}\rightarrow{})".format(observable_str,qn_str_2,qn_str_1)  # <1|O|2> = 2->1
    else:
        raise(ValueError("unrecognized observable type {}".format(observable_type)))
                                                
    return label

HW_AXIS_LABEL_TEXT = r"$\hbar\omega~(\mathrm{MeV})$"

def make_observable_axis_label_text(nuclide_observable,bare_math_mode_text=False):
    """ Generate axis label (with units) for observable.

    Arguments:

        nuclide_observable (tuple): standard nuclide/observable pair or compound

    Returns:
    
        label (str): label string, to be interpreted as mpl text
    """
    # trap compound observable
    if nuclide_observable[0] in {"diff","ratio"}:
        (arithmetic_operation,nuclide_observable1,nuclide_observable2) = nuclide_observable
        if arithmetic_operation == "diff":
            return r"$\Delta {}$".format(make_observable_axis_label_text(nuclide_observable1,bare_math_mode_text=True))
        elif arithmetic_operation == "ratio":
            return r"Ratio"

    # unpack arguments
    (nuclide,observable) = nuclide_observable
    (observable_type,observable_operator,observable_qn_list) = unpack_observable(observable)

    # construct label
    if observable_type == "energy":
        observable_str = r"E"
        units_str = r"\mathrm{MeV}"
    elif observable_type == "isospin":
        observable_str = r"\bar{T}"
        units_str = None
    elif observable_type == "radius":
        observable_str = r"r"
        units_str = r"\mathrm{fm}"
    elif observable_type in {"moment","momentsqr"}:
        if observable_operator in {"M1","Dlp","Dln","Dsp","Dsn","Dl0","Dl1","Ds0","Ds1"}:
            observable_str = r"\mu"
            units_str = r"\mu_N"
        elif observable_operator in {"E2p","E2n","E20","E21"}:
            observable_str = r"eQ"
            units_str = r"e\,\mathrm{fm}^{2}"
        if observable_type == "momentsqr":
            observable_str = r"({})^2".format(observable_str)
    elif observable_type == "rtp":
        if observable_operator in {"M1","Dlp","Dln","Dsp","Dsn","Dl0","Dl1","Ds0","Ds1"}:
            observable_str = r"B(M1)"
            units_str = r"\mu_N"
        elif observable_operator in {"E2p","E2n","E20","E21"}:
            observable_str = r"B(E2)"
            units_str = r"e^2\,\mathrm{fm}^{4}"
    else:
        raise(ValueError("unrecognized observable type {}".format(observable_type)))

    if units_str is not None:
        label = r"{}~({})".format(observable_str,units_str)
    else:
        label = observable_str
    if not bare_math_mode_text:
        label = r"${}$".format(label)

    return label

def make_interaction_text(interaction_coulomb):
    """ Interaction text.

    Arguments:
    
        interaction_coulomb (tuple): interaction/coulomb specifier

    Returns:
    
        label (str): label string, to be interpreted in math mode
    """
    label = r"\mathrm{{{}}}".format(interaction_coulomb[0])
    return label

################################################################
# Nmax line styling
################################################################

def Nmax_dashing(Nmax_relative,base_length=8,exponent_scale=8):
    """Provide dashing pattern based on relative Nmax.

    Inspired by Zhou's dashing:

        (4*r, 4*(1-r))
        r=nmax/max_Nmax  (goes from 1 at highest Nmax, to 0 at floor Nmax0)

    But here we use logarithmic scale to maximum Nmax need not be specified, and
    line for Nmax0 would not vanish.

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:

        (tuple): dashing directive

    """
    if Nmax_relative == 0:
        return (None,None)
    elif Nmax_relative < 0:
        r = 2**(Nmax_relative/exponent_scale)
        return (base_length*r,base_length*(1-r))

    elif Nmax_relative == -4:
        return (2,2)
    elif Nmax_relative < -4:
        return (1,2)

    raise ValueError("invalid Nmax_relative {}".format(Nmax_relative))

def Nmax_dashing_scidraw(Nmax_relative):
    """ Provide dashing pattern based on relative Nmax.

    Follows RotationPlots.m "NmaxDashing" dashing convention.

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:

        (tuple): dashing directive
    """
    if Nmax_relative == 0:
        return (None,None)
    elif Nmax_relative == -2:
        return (6,6)
    elif Nmax_relative == -4:
        return (2,2)
    elif Nmax_relative < -4:
        return (1,2)

    raise ValueError("invalid Nmax_relative {}".format(Nmax_relative))

def Nmax_color(Nmax_relative):
    """ Provide color based on relative Nmax.

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:
    
        (str): color directive
    """
    if Nmax_relative == 0:
        return "firebrick"
    elif Nmax_relative < 0:
        return "black"

    raise ValueError("invalid Nmax_relative {}".format(Nmax_relative))

def Nmax_symbol_scale(Nmax_relative,exponent_scale=8):
    """ Provide symbol scale based on relative Nmax.

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

    Returns:

        (float): symbol scale factor (max 1.0)
    """

    if Nmax_relative <= 0:
        return 2**(Nmax_relative/exponent_scale)

    raise ValueError("invalid Nmax_relative {}".format(Nmax_relative))

def Nmax_plot_style(Nmax_relative,marker_size=6):
    """

    Arguments:

        Nmax_relative (int): Nmax-Nmax_max

        marker_size (float,optional): base marker size for maximal Nmax

    Returns:

        (dict): kwargs for ax.plot
    """
    return dict(
        markersize=marker_size*Nmax_symbol_scale(Nmax_relative),
        linewidth=1,
        dashes=Nmax_dashing(Nmax_relative),
        color=Nmax_color(Nmax_relative),
        )

################################################################
# hw scan
################################################################

def hw_scan_descriptor(interaction_coulomb,nuclide_observable):
    """ Generate standard descriptor string for a (nuclide,observable) pair.

    Arguments:
        nuclide_observable (tuple): standard nuclide/observable pair or compound

    Returns:
        descriptor (str): descriptor string
    """
    descriptor="hw-scan_{interaction_coulomb[0]:s}-{interaction_coulomb[1]:1d}_{nuclide_observable_descriptor}".format(
        interaction_coulomb=interaction_coulomb,
        nuclide_observable_descriptor=nuclide_observable_descriptor(nuclide_observable)
    )

    return descriptor

def make_hw_scan_data(mesh_data,nuclide_observable,verbose=False):
    """ Tabulate generic observable for hw scan.

    Tabulation format:
        Nmax hw value

    Arguments:
        mesh_data (list of ResultsData): data set to include
        ##interaction_coulomb (tuple): (interaction,use_coulomb)
        nuclide_observable (tuple): (nuclide,observable)
          observable: (observable_type,observable_operator,(Jf,gf,nf),...)

    Returns:
        observable_data (np.array): scan data, with rows (Nmax,hw,value)
    """

    # trap compound observable
    if nuclide_observable[0] in {"diff","ratio"}:
        (arithmetic_operation,nuclide_observable1,nuclide_observable2) = nuclide_observable
        if arithmetic_operation == "diff":
            return make_hw_scan_data(mesh_data,nuclide_observable1)-make_hw_scan_data(mesh_data,nuclide_observable2)
        elif arithmetic_operation == "ratio":
            return make_hw_scan_data(mesh_data,nuclide_observable1)/make_hw_scan_data(mesh_data,nuclide_observable2)

    # unpack arguments
    (nuclide,observable) = nuclide_observable
    (observable_type,observable_operator,observable_qn_list) = unpack_observable(observable)

    # select nuclide
    mesh_data_selected = analysis.selected_mesh_data(mesh_data,{"nuclide":nuclide})
    analysis.mesh_key_listing(mesh_data_selected,("nuclide","interaction","coulomb","hw","Nmax","parity"),verbose=verbose)

    # generate table

    # NOTE: May ultimately supplant analysis.make_obs_table ndarray step with
    # direct construction of pandas data frame.
    
    KEY_DESCRIPTOR_NMAX_HW = (("Nmax",int),("hw",float))
    if observable_type == "energy":
        table = analysis.make_obs_table(
            mesh_data_selected,KEY_DESCRIPTOR_NMAX_HW,
            lambda results_data : results_data.get_energy(*observable_qn_list)
        )
    elif observable_type == "isospin":
        table = analysis.make_obs_table(
            mesh_data_selected,KEY_DESCRIPTOR_NMAX_HW,
            lambda results_data : results_data.get_isospin(*observable_qn_list)
        )
    elif observable_type == "radius":
        table = analysis.make_obs_table(
            mesh_data_selected,KEY_DESCRIPTOR_NMAX_HW,
            lambda results_data : results_data.get_radius(observable_operator,*observable_qn_list)
        )
    elif observable_type == "moment":
        table = analysis.make_obs_table(
            mesh_data_selected,KEY_DESCRIPTOR_NMAX_HW,
            lambda results_data : results_data.get_moment(observable_operator,*observable_qn_list)
        )
    elif observable_type == "momentsqr":
        table = analysis.make_obs_table(
            mesh_data_selected,KEY_DESCRIPTOR_NMAX_HW,
            lambda results_data : results_data.get_moment(observable_operator,*observable_qn_list)**2
        )
    elif observable_type == "rtp":
        table = analysis.make_obs_table(
            mesh_data_selected,KEY_DESCRIPTOR_NMAX_HW,
            lambda results_data : results_data.get_rtp(observable_operator,tuple(observable_qn_list))
        )
    else:
        raise(ValueError("unrecognized observable type {}".format(observable_type)))

    # convert to DataFrame
    observable_data = pd.DataFrame(table).set_index(["Nmax","hw"])

    # drop nan values
    observable_data = observable_data[observable_data["value"].notna()]
    
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
        ax,nuclide_observable,hw_range,observable_range,
        hw_range_extension=(0.05,0.05),observable_range_extension=(0.05,0.05)
):
    """ Set up axes.

    Arguments:
    
        ax (mpl): descriptor string to use in filename

        nuclide_observable (tuple): standard nuclide/observable pair or compound

        hw_range (tuple of float): x range, before extension

        observable_range (tuple of float): y range, or None for matplotlib auto

        hw_range_extension (tuple of float, optional): x range relative extension

    """
    ax.set_xlabel(HW_AXIS_LABEL_TEXT)
    ax.set_xlim(*extend_interval_relative(hw_range,hw_range_extension))
    ax.set_ylabel(make_observable_axis_label_text(nuclide_observable))
    if (observable_range is not None) and np.isfinite(observable_range[0]).all():
        ax.set_ylim(*extend_interval_relative(observable_range,observable_range_extension))

def add_observable_panel_label(ax,interaction_coulomb,nuclide_observable,**kwargs):
    """ Add observable panel label to plot.

    Arguments:
    
        ax (mpl): axes object

        interaction_coulomb (tuple): interaction/coulomb specifier

        nuclide_observable (tuple): standard nuclide/observable pair or compound

    """

    # panel label
    ax.annotate(
        "${}$ ${}$\n$^{}$".format(
            make_nuclide_text(nuclide_observable),
            make_observable_text(nuclide_observable),
            make_interaction_text(interaction_coulomb)
        ),
        xy=(0.95,0.05),xycoords="axes fraction",
        multialignment="left",
        horizontalalignment="right",
        verticalalignment="bottom",
        bbox=dict(boxstyle="round",facecolor="white")
    )
        
def add_hw_scan_plot(ax,observable_data,Nmax_max):
    """ Add hw scan plot to axes.

    Arguments:
    
        ax (mpl): descriptor string to use in filename

        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

        Nmax_max (int): highest Nmax for styling purposes
    """

    for Nmax, group in observable_data.reset_index().groupby("Nmax"):
        ax.plot(
            group["hw"],group["value"],
            marker=".",
            **Nmax_plot_style(Nmax-Nmax_max)
        )

def write_hw_scan_plot(
        descriptor,
        interaction_coulomb,nuclide_observable,
        observable_data,
        hw_range,observable_range,Nmax_max,
        hw_range_extension=(0.02,0.02),
        observable_range_extension=(0.02,0.02),
        figsize=(6,4),
        directory="data"
):
    """ Generate full "canned" hw scan plot.

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
    add_observable_panel_label(ax,interaction_coulomb,nuclide_observable)

    # generate plot
    add_hw_scan_plot(ax,observable_data,Nmax_max)

    # finalize plot
    figure_file_name = os.path.join(
        directory,
        "{}_plot.pdf".format(descriptor)
        )
    plt.savefig(figure_file_name)
    plt.close()
        
