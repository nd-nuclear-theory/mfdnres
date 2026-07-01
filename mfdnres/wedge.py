"""wedge.py

    Plotting functions for beta-gamma deformation wedge plots.

    Mark A. Caprio
    Scott R. Carmichael
    University of Notre Dame

    - 02/20/24 (mac): Created, from code in c_su3.py (src).
    - 04/29/25 (mac): Add option share_probability to add_decomposition_wedge_plot.
    - 07/01/26 (mac): Absorb into mfdnres.
"""

import itertools
import sys

import matplotlib as mpl
## import matplotlib.pyplot as plt
import numpy as np

import mfdnres.data
import mfdnres.decomposition


################################################################
# helper functions
################################################################

cast_to_su3 = mfdnres.decomposition.labels_subsetting_function(mfdnres.decomposition.SU3Labels)


################################################################
# axis setup
################################################################

def set_up_wedge_plot_axes(
        ax,
        *,
        beta_max=1,
        linewidth = 1,
        fontsize="large",
        show_beta_label=True,
        show_gamma_label=True,
        show_axes=False,
):
    """ Set up axis ranges, wedge, and wedge labels for wedge plot.

    Arguments:

        ax (mpl.axes.Axes): Axes object.

        beta_max (float, optional): Radial axis range.

        linewidth (float, optional): Line width for wedge.

        fontsize (str, optional): Font size for beta and gamma labels.

        show_beta_label (bool, optional): Whether or not to show beta label.

        show_gamma_label (bool, optional): Whether or not to show gamma label.

        show_axes (bool, optional): Whether or not to show axis frame (for
            debugging purposes).

    """

    # draw wedge
    #
    max_angle_deg = 60
    mid_angle_deg = 30
    ax.add_artist(
        mpl.patches.Wedge(
            (0,0), beta_max, 0, max_angle_deg,
            facecolor="none", edgecolor="black",
            linewidth=linewidth,
        ),
    )
    ax.add_line(
        mpl.lines.Line2D(
            [0, beta_max*np.cos(mid_angle_deg*np.pi/180)],
            [0, beta_max*np.sin(mid_angle_deg*np.pi/180)],
            linestyle=":", linewidth=linewidth, color="dimgray",
            zorder=0,  # to lie behind wedge and behind scatter plots
        ),
    )
    
    # set frame options
    if not show_axes:
        ax.set_axis_off()
    ax.set(
        xlim=(-0.01,1.01*beta_max),
        ylim=(-0.15,1.01*np.sin(60*np.pi/180)*beta_max),
    )
    ax.set_aspect('equal', adjustable='box')

    # add beta/gamma labels
    if show_beta_label:
        beta_label_text = r"\beta"
        ## ax.set_xlabel(r"${}$".format(x_axis_label_text), fontsize=20)
        ax.annotate(
            r"${}$".format(beta_label_text), 
            xy = (0.5*beta_max, -0.01*beta_max), 
            verticalalignment="top",
            fontsize=fontsize,
        )
    if show_gamma_label:
        gamma_label_text = r"\gamma"
        ax.annotate(
            r"${}$".format(gamma_label_text), 
            xy = ((np.sqrt(3)/2)*beta_max, 0.5*beta_max), 
            rotation = 30,
            fontsize=fontsize,
        )


################################################################
# data extraction
################################################################

def extract_leading_irrep_evolution(
        leading_irreps_by_nuclide_Nex,
        irrep_provenance_by_u3s_labels_by_nuclide_Nex,
        nuclide_list, Nex,
        ):
    """Extract SU(3) leading irreps (with provenance) for chain of nuclides.

    Arguments:

        ...

    Returns:

        list[tuple]: List of (nuclide, proton_irreps, neutron_irreps, irreps).
            Due to possible ambiguous provenance or degenerate leading irrep,
            irreps are given as tuples containing groups of either one or two
            irreps.

    """

    evolution_data = []
    for nuclide in nuclide_list:
        leading_irreps_u3s = leading_irreps_by_nuclide_Nex[(nuclide, Nex)]
        proton_irreps = set()
        neutron_irreps = set()
        irreps = set()
        for irrep_u3s in leading_irreps_u3s:
            provenance = irrep_provenance_by_u3s_labels_by_nuclide_Nex[(nuclide, Nex)][irrep_u3s]
            for proton_irrep_u3s, neutron_irrep_u3s in provenance:
                proton_irreps.add(cast_to_su3(proton_irrep_u3s))
                neutron_irreps.add(cast_to_su3(neutron_irrep_u3s))
            irreps.add(cast_to_su3(irrep_u3s))
        proton_irreps = tuple(sorted(proton_irreps))
        neutron_irreps = tuple(sorted(neutron_irreps))
        irreps = tuple(sorted(irreps))
        assert(len(proton_irreps)<=2 and len(neutron_irreps)<=2 and len(irreps)<=2)
        
        evolution_data.append((nuclide, proton_irreps, neutron_irreps, irreps))

    return evolution_data


################################################################
# plotting
################################################################

def irrep_beta_gamma(irrep):
    """For a given irrep, calculate beta and gamma.

    Here beta is calculated sans the dimensionful, A-dependent scale factor.

    Arguments:

        irrep (collections.namedtuple): Irrep with lambda_omega and mu_omega (or
            lambda_sigma and mu_sigma) fields.

    Returns:

        tuple[float]: beta, gamma

    """

    if "lambda_omega" in irrep._fields:
        lam = irrep.lambda_omega
        mu = irrep.mu_omega
    elif "lambda_sigma" in irrep._fields:
        lam = irrep.lambda_sigma
        mu = irrep.mu_sigma
    else:
        raise ValueError("Irrep does not contain suitable lambda and mu fields.")
    
    beta = np.sqrt(lam**2 + lam*mu + mu**2 + 3*lam + 3*mu + 3)
    gamma = np.arctan(np.sqrt(3)*(mu + 1) / (2*lam + mu + 3))
    
    return beta, gamma


COLOR_BY_SPECIES = {
    "m": "blue",
    "p": "red",
    "n": "dimgray",
}
TEXT_BY_SPECIES = {
    "m": "Matter",
    "p": "Proton",
    "n": "Neutron",
}
LEGEND_HANDLES_FOR_SPECIES = [
    mpl.lines.Line2D(
        [0], [0],
        marker="o",
        color=COLOR_BY_SPECIES[species],
        ## color="white",  # kludge to suppress line
        ## markeredgecolor="darkgreen",markerfacecolor="darkgreen",
        markersize=8,
        label=TEXT_BY_SPECIES[species],
    )
    for species in ["m", "p", "n"]
]


def add_evolution_wedge_plot(
        ax,
        evolution_data,
        *,
        species="m",
        color=None,
        **kwargs,
):
    """Draw (beta,gamma) evolution wedge plot.

    Arguments:

        ax (mpl.axes.Axes): Axes object.

        evolution_data (list[tuple]): Evolution data as returned by
            extract_leading_irrep_evolution.
 
        species (str, optional): Species as "m" (matter), "p" (proton), or "n" (neutron).

        color (str, optional): Color for plot symbols.  Default determined from species.

        kwargs (Line2D properties, optional): kwargs are used to specify line
        properties not otherwise fixed by the prior arguments.

    Returns:

        xy_by_nuclide_by_branch (dict): Database of plotted points.  Dictionary
            signature is:

            branch (0 or -1) -> nuclide (N, Z) -> coordinates (x, y)

        TODO (mac): Restructure return database to better match xyr_by_irrep,
        for use by annotation functions?

    """

    xy_by_nuclide_by_branch = dict()
    for branch_index in [0, -1]:  # pick member of each group (singlet or pair) of alternative irreps

        # accumulate points
        x_list = []
        y_list = []
        ## xy_by_nuclide = xy_by_nuclide_by_branch.get(branch_index, dict())
        xy_by_nuclide_by_branch.setdefault(branch_index, dict())
        xy_by_nuclide = xy_by_nuclide_by_branch[branch_index]
        for nuclide, proton_irreps, neutron_irreps, irreps in evolution_data:
            if species == "m":
                irrep = irreps[branch_index]
            elif species == "p":
                irrep = proton_irreps[branch_index]
            elif species == "n":
                irrep = neutron_irreps[branch_index]
            beta, gamma = irrep_beta_gamma(irrep)
            x, y = beta*np.cos(gamma), beta*np.sin(gamma)
            x_list.append(x)
            y_list.append(y)
            xy_by_nuclide[nuclide] = (x, y)
            
        # draw points
        if color is None:
            color = COLOR_BY_SPECIES[species]
        ax.plot(
            x_list, y_list,
            color=color, marker='o',
            **kwargs,
        )

    return xy_by_nuclide_by_branch


def add_nuclide_labels(
        ax,
        xy_by_nuclide_by_branch,
        *,
        label_style="isotope",
        species="m",
        color=None,
        nuclide_list=None,
        label_displacement = (0, 0),
        **kwargs,
):
    """Add nuclide labels to leading irrep evolution wedge plot.

    Arguments:

        ax (mpl.axes.Axes): Axes object.

        xy_by_nuclide_by_branch (dict): Database of plotted points as returned by
            add_evolution_wedge_plot.
 
        species (str, optional): Species as "m" (matter), "p" (proton), or "n" (neutron).
           Used here only for color selection.

        color (str, optional): Color for plot labels.  Default determined from species.

        nuclide_list (list[tuple], optional): List of nuclides to label.
        Defaults to all.

        label_displacement (tuple or dict): Displacement of label (in
        point) from symbol center, or dict of such displacements by nuclide.

        kwargs (Line2D properties, optional): kwargs are used to specify line
        properties not otherwise fixed by the prior arguments.

    """

    for branch_index in [0, -1]:  # pick member of each group (singlet or pair) of alternative irreps
        xy_by_nuclide = xy_by_nuclide_by_branch[branch_index]

        if color is None:
            color = COLOR_BY_SPECIES[species]

        if nuclide_list is None:
            nuclide_list = list(xy_by_nuclide.keys())

        # annotate points with nuclide labels
        for nuclide in nuclide_list:
            xy = xy_by_nuclide[nuclide]
            if label_style=="isotope":
                label_text = mfdnres.data.isotope(nuclide)
            elif label_style=="A":
                label_text = str(sum(nuclide))
            else:
                raise(ValueError("unrecognized label style {}".format(label_style)))
            if type(label_displacement) is tuple:
                the_label_displacement = label_displacement
            elif type(label_displacement) is dict:
                the_label_displacement = label_displacement[nuclide]
            else:
                raise(ValueError("unrecognized label_displacement {}".format(label_displacement)))
            ax.annotate(
                r"${}$".format(label_text),
                xy=xy,
                xytext=the_label_displacement,
                textcoords="offset points",
                color=color,
                **kwargs,
                )


def add_decomposition_wedge_plot(
        ax,
        decomposition,
        *,
        scale=1000,
        connector_kwargs={},
        share_probability=True,
        **kwargs,
):
    """Draw (beta,gamma) evolution wedge plot.

    Arguments:

        ax (mpl.axes.Axes): Axes object.

        decomposition (dict): SU3 decomposition, with signature labels
        (mfdnres.decomposition.SU3Labels) -> probability (float).
 
        scale (float, optional): Scale factor from probability to symbol area.

        share_probability (bool, optional): Whether probability should be shared
        equally over members of degenerate irrep group, or total probability
        should be shown for each (so displayed probabilities no longer sum to
        total probability).

        kwargs (Line2D properties, optional): kwargs are used to specify line
        properties not otherwise fixed by the prior arguments.

    Returns:

        xyr_by_irrep (dict): Database of plotted points, as irrep
        (mfdnres.decomposition.SU3Labels) -> (xy, radius).

    """

    xyr_by_irrep = dict()
    x_list = []
    y_list = []
    s_list = []
    for irrep_group, probability in decomposition.items():
        # extract data for irreps within group
        x_list_group = []
        y_list_group = []
        s_list_group = []
        if share_probability:
            probability /= len(irrep_group)
        for irrep in irrep_group:
            beta, gamma = irrep_beta_gamma(irrep)
            x, y = beta*np.cos(gamma), beta*np.sin(gamma)
            x_list_group.append(x)
            y_list_group.append(y)
            area = scale*probability
            s_list_group.append(area)
            irrep = cast_to_su3(irrep)  # ensure key is SU3 (may be, e.g., U3, in given decomposition)
            xyr_by_irrep[irrep] = ((x, y), np.sqrt(area/np.pi))

        # draw connector line for degenerate irreps
        if len (irrep_group)>1:
            ax.plot(
                x_list_group, y_list_group,
                marker="",
                zorder=0,
                **connector_kwargs,
            )

        # accumulate group data to master list
        x_list += x_list_group
        y_list += y_list_group
        s_list += s_list_group
            
    # draw points
    ax.scatter(
        x_list, y_list, s=s_list,
        marker="o",
        **kwargs,
    )

    return xyr_by_irrep


def add_dimension_wedge_plot(
        ax,
        dimension_data,
        *,
        scale=np.pi*4**2,
        unit_dimension=False,
        **kwargs,
):
    """Draw (beta,gamma) basis dimension (or irreps) wedge plot.

    The task is essentially a special case of that performed by
    add_decomposition_wedge_plot(), except that each irrep label is unique, not
    part of a degenerate tuple, and degeneracies need not be taken into account
    in plotting.

    Also, the symbols can be forced to a uniform size, to show only the lattice
    of irreps, not their dimensions.

    Arguments:

        ax (mpl.axes.Axes): Axes object.

        dimensions (dict): SU3 (or U3) dimensions, with signature labels
        (mfdnres.decomposition.SU3Labels) -> probability (float).
 
        scale (float, optional): Scale factor from dimension to symbol area.

        unit_dimension (bool, optional): Whether to show actual dimension
        (False) or to show all irreps with unit dimension (True), for equal
        sized symbols.

        kwargs (Line2D properties, optional): kwargs are used to specify line
        properties not otherwise fixed by the prior arguments.

    Returns:

        xyr_by_irrep (dict): Database of plotted points, as nuclide->(xy, radius).

    """

    xyr_by_irrep = dict()
    x_list = []
    y_list = []
    s_list = []
    for irrep, dimension in dimension_data.items():
        beta, gamma = irrep_beta_gamma(irrep)
        x, y = beta*np.cos(gamma), beta*np.sin(gamma)
        x_list.append(x)
        y_list.append(y)
        if unit_dimension:
            dimension = 1
        area = scale*dimension
        s_list.append(area)
        irrep = cast_to_su3(irrep)  # ensure key is SU3 (may be, e.g., U3, in given dimensions)
        xyr_by_irrep[irrep] = ((x, y), np.sqrt(area/np.pi))
            
    # draw points
    ax.scatter(
        x_list, y_list, s=s_list,
        marker="o",
        **kwargs,
    )

    return xyr_by_irrep


def add_irrep_labels(
        ax,
        xyr_by_irrep,
        *,
        irrep_list=None,
        label_displacement=None,
        label_padding=2,
        verbose=False,
        **kwargs,
):
    """Add irrep labels to wedge plot.

    Arguments:

        ax (mpl.axes.Axes): Axes object.

        xyr_by_irrep (dict): Database of plotted points, as nuclide->(xy, radius).
 
        irrep_list (list[tuple] or list[(mfdnres.decomposition.SU3Labels],
        optional): List of irreps to label.  Defaults to all.  Beware that
        values of quantum numbers may need to be float, even if conceptually
        they are integers, to successfully match those in the underlying
        xyr_by_irrep database.

        label_displacement (tuple or dict, optional): Displacement of label (in
        point) from symbol center, or dict of such displacements by nuclide.
        Default (or a value of None) placed label to right of irrep at position
        given by radius of symbol.

        label_padding (float, optional): Padding of label (in point) from nominal
        position, when calculated automatically.

        kwargs (Line2D properties, optional): kwargs are used to specify line
        properties not otherwise fixed by the prior arguments.  The following
        defaults are provided: horizontalalignment="left",
        verticalalignment="center_baseline".

    """

    if irrep_list is None:
        irrep_list = list(xyr_by_irrep.keys())
    irrep_list = list(map(cast_to_su3, irrep_list))  # cast irreps to SU3
    if verbose:
        print("  Available irreps: {}".format(list(xyr_by_irrep.keys())))
        print("  Labeling irreps: {}".format(irrep_list))

    default_kwargs = dict(horizontalalignment="left", verticalalignment="center_baseline")
    full_kwargs = default_kwargs | kwargs

    # annotate points with irrep labels
    for irrep in irrep_list:
        irrep = cast_to_su3(irrep)
        if irrep not in xyr_by_irrep:
            continue
        xy, r = xyr_by_irrep[irrep]
        if type(label_displacement) is tuple or label_displacement is None:
            the_label_displacement = label_displacement
        elif type(label_displacement) is dict:
            the_label_displacement = label_displacement[nuclide]
        else:
            raise(ValueError("unrecognized label_displacement {}".format(label_displacement)))
        if the_label_displacement is None:
            the_label_displacement = (r+label_padding, 0)
        label_text = str(cast_to_su3(irrep))
        ax.annotate(
            r"${:s}$".format(label_text),
            xy=xy,
            xytext=the_label_displacement,
            textcoords="offset points",
            **full_kwargs,
            )


def redraw_irreps(
        ax,
        xyr_by_irrep,
        *,
        irrep_list=None,
        verbose=False,
        **kwargs,
):
    """Highlights previously drawn irreps by redrawing them with given options.

    Accepts the same options as the original scatter plot.

    Arguments:

        ax (mpl.axes.Axes): Axes object.

        xyr_by_irrep (dict): Database of plotted points, as nuclide->(xy, radius).
 
        irrep_list (list[tuple] or list[SU3Labels], optional): List of irreps to
        highlight.  Defaults to all.

        kwargs (Collection properties, optional): kwargs are used to specify
        scatter plot properties not otherwise fixed by the prior arguments.

    """

    if irrep_list is None:
        irrep_list = list(xyr_by_irrep.keys())
    irrep_list = list(map(cast_to_su3, irrep_list))  # cast irreps to SU3
    if verbose:
        print("  Available irreps: {}".format(list(xyr_by_irrep.keys())))
        print("  Redrawing irreps: {}".format(irrep_list))
    
    # redraw points
    x_list = []
    y_list = []
    s_list = []
    for irrep in irrep_list:
        (x, y), r = xyr_by_irrep[irrep]
        x_list.append(x)
        y_list.append(y)
        area = np.pi*r**2
        s_list.append(area)
        
    ax.scatter(
        x_list, y_list, s=s_list,
        marker="o",
        **kwargs,
    )


################################################################
# main
################################################################

def main():
    """
    """
    beta, gamma = irrep_beta_gamma(mfdnres.decomposition.U3Labels(0, 2, 2))
    print("{} {}".format(beta, gamma*180/np.pi))
    
if __name__ == "__main__":
    main()

