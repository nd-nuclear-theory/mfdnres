"""network.py

    Network plotting machinery.

    Mark A. Caprio
    University of Notre Dame


    - 04/27/23 (mac): Created, incorporating code extracted from 12be-shape testbed.
    - 05/01/23 (mac): Add level selection, transition initial level selection, and valence marker line.

"""

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from . import (
    data,
    level,
    ticks,
)

################################################################
# level tabulation selector functions
################################################################

def level_test_0hw(results_data, qn, verbose=False):
    """Test if level is 0hw.

    Arguments:

        results_data (mfdnres.ResultsData): Results data

        qn (tuple): (J,g,n) quantum nubmers

    Returns:

        (bool): Result of test

    """

    Nex_decomposition = results_data.get_decomposition("Nex", qn)
    # Note: Provide short-circuit test in Nmax=0 space, where all states are
    # deemed 0hw, to avoid out-of-range index on Nex_decomposition.
    selected = (len(Nex_decomposition)==1) or (Nex_decomposition[0]>=Nex_decomposition[1])
    return selected

def level_test_T_excited(results_data,verbose=False):
    """Test if level has excited isospin, i.e., a deduced T > Tz.

    The threshold on deduced T is Tz+0.5.

    Arguments:

        results_data (mfdnres.ResultsData): Results data

        qn (tuple): (J,g,n) quantum nubmers

    Returns:

        (bool): Result of test

    """

    nuclide = results_data.params["nuclide"]
    Tz=1/2*(nuclide[0]-nuclide[1])
    T = results_data.get_isospin(qn)
    selected = (T-abs(Tz))>0.5
    return selected

################################################################
# network plotting
################################################################

def jj1(J):
    """Helper function providing J->J*(J+1)."""
    return J*(J+1)

def set_up_network_axes(
        ax, J_max, E_range,
        J_range_extension=(0.05,0.05), E_range_extension=(0.05,0.05),
        E_tick_specifier=None,
        show_frame=True,
):
    """Set up axis ranges, labels, and ticks for network plot.

    Arguments:

        ax (mpl.axes.Axes): axes object

        J_max (float): highest J for horizontal axis

        E_range (tuple of float): (E_min, E_max) for vertical axies

        J_range_extension (tuple of float, optional): horizontal axis range
        relative extension; this is applied after J(J+1) scaling

        E_range_extension (tuple of float, optional): vertical axis range
        relative extension

        E_tick_specifier (tuple, optional): tick specification
        (min,max,step,num_subdivision) for vertical ticks

        show_frame (bool, optional): whether or not to draw frame

    """

    if show_frame:
        # set ticks
        #
        # Note that the tick specification must come *before* setting range limits,
        # since the range limits automatically readjust when you set the ticks.
        # 
        # TODO (mac): implement tick_post_transformation and tick_label_function in
        # linear_ticks
        x_ticks = [J for J in range(0,int(J_max)+1,1)]  # TODO (mac): currently only provides integer ticks
        x_tick_values = [x*(x+1) for x in x_ticks]
        ax.xaxis.set_major_formatter(ticks.HalfIntFormatter())
        ## x_tick_labels=[mfdnres.ticks.half_int_str(x) for x in x_ticks]
        ax.set_xticks(x_tick_values, x_ticks)
        if E_tick_specifier is not None:
            y_ticks = ticks.linear_ticks(*E_tick_specifier)
            ticks.set_ticks(ax,"y",y_ticks)
        
        # set axis labels
        #
        # TODO (mac): support alternative labels, e.g., absolute "E"
        ax.set_xlabel(r"$J$",labelpad=1)
        ax.set_ylabel("$E_x~(\mathrm{MeV})$",labelpad=1)
    else:
        # disable frame and ticks
        ax.set_xticks([], [])
        ax.set_yticks([], [])
        ax.set_frame_on(False)

    # set limits
    ax.set_xlim(*data.extend_interval_relative((0, jj1(J_max)), J_range_extension))
    ax.set_ylim(*data.extend_interval_relative(E_range, E_range_extension))

    
def select_network_levels(
        results_data,
        reference_energy=None,
        reference_results_data=None,
        J_max=None,
        E_max=None,
        n_max_for_J=None,
        verbose=False,
):
    """Select levels for network plot (with plotting data).

    Arguments:

        results_data (mfdnres.ResultsData): Results data

        reference_energy (float or level.Level, optional): energy, or level providing zero
        point for energy, or None to use raw energies

        reference_results_data (mfdnres.ResultsData, optional): Results data for
        use in extracting reference energy (e.g., may contain spectrum of
        different parity or of isobar)

        J_max (float, optional): cap on J

        E_max (float, optional): cap on E 

        n_max_for_J (dict): mapping J->n_max giving cap on number of states for
        each J

    Returns:

        network_levels (dict): Plotting data for levels, stored as qn->E

    """

    # find reference energy
    if reference_results_data is None:
        reference_results_data = results_data
    if reference_energy is None:
        E0 = 0
    elif isinstance(reference_energy, level.Level):
        E0 = reference_results_data.get_energy(reference_energy.select_level(reference_results_data))
    else:
        E0 = reference_energy
    if verbose:
        print("Reference energy {}".format(E0))

    # collect levels to include in network
    #
    network_levels = {}
    for qn in results_data.levels:

        # retrieve level energy
        J, g, n = qn
        E = results_data.get_energy(qn) - E0

        # filter level
        if not abs(J%0.5)==0.:  # filter by (half) integrity (to suppress some unconverged levels)
            continue
        if J_max is not None and J > J_max:  # filter by J limit
            continue
        if E_max is not None and E>E_max:  # by max energy (relative to reference)
            continue
        if n_max_for_J is not None:  # by max n for given J
            n_max = n_max_for_J.get(J)
            if n_max is not None and n>n_max:
                continue
            
        # save level
        network_levels[qn] = E
    
    if verbose:
        print("Network levels: {}".format(network_levels))
    return network_levels


def filter_network_levels(
        results_data,
        network_levels,
        test,
        verbose=False,
):
    """Filter levels for network plot (with plotting data).

    Arguments:

        results_data (mfdnres.ResultsData): Results data

        network_levels (dict): Plotting data for levels, stored as qn->E

        test (callable): Test function with signature test(results_data,qn)->bool

    Returns:

        filtered_network_levels (dict): Plotting data for levels, stored as qn->E
    """

    filtered_network_levels = dict()
    for qn, E in network_levels.items():
        if test(results_data, qn):
            filtered_network_levels[qn] = E
    return filtered_network_levels

def draw_network_valence_marker(
        ax,Jv,
        plot_kw={},
):
    """ Draw vertical marker line for maximal valence angular momentum.

    Arguments:

        ax (mpl.axes.Axes): axes object

        Jv (float): maximal valence angular momentum

        plot_kw (dict, optional): Line2D properties properties overriding
        the default line style
    """
    # draw levels
    plot_kw_defaults = dict(
        color="dimgray",
        linestyle="dashed",
        linewidth=0.5,
        marker=None,
    )
    plot_kw_full = {
        **plot_kw_defaults,
        **plot_kw,
    }
    ax.axvline(jj1(Jv), **plot_kw_full)
    
def draw_network_levels(
        ax,
        network_levels,
        plot_kw={},
        verbose=False,
):
    """Draw levels for network plot.

    Arguments:

        ax (mpl.axes.Axes): axes object

        network_levels (dict): Plotting data for levels, stored as qn->E

        plot_kw (dict, optional): Line2D properties properties overriding
        the default marker style (e.g., marker="s", markerfacecolor = "white",
        markeredgecolor="black", markersize=7)

    """

    # extract coordinates for plotting
    #
    # Levels are ordered by (J, E) for plotting.
    JJ1_values = [
        jj1(J)
        for J, _, _ in network_levels.keys()
    ]
    E_values = [
        E
        for E in network_levels.values()
    ]
    if verbose:
        print("JJ1 {} E {}".format(JJ1_values,E_values))
    
    # marker face coloring
    ## if marker_face_coloring is not None:
    ##     print(marker_face_coloring)
    ##     test, true_color, false_color = marker_face_coloring
    ##     color_values = [
    ##         true_color if test(results_data, qn) else false_color
    ##         for qn in network_levels.keys()
    ##     ]
    ##     print("color_values: {}".format(color_values))
        
    # draw levels
    plot_kw_defaults = dict(
        linestyle="None",
        marker="s",
        markerfacecolor = "white",
        markeredgecolor="black",
        markersize=4,
    )
    plot_kw_full = {
        **plot_kw_defaults,
        **plot_kw,
    }
    
    ax.plot(JJ1_values, E_values, **plot_kw_full)


def draw_expt_levels(
        ax,
        expt_energies,
        plot_kw={},
        verbose=False,
):
    """Draw levels for network plot.

    Arguments:

        ax (mpl.axes.Axes): axes object

        expt_energies (list): Plotting data for levels, stored as [J, E]

        plot_kw (dict, optional): Line2D properties properties overriding
        the default marker style

    """

    # extract coordinates for plotting
    #
    # Levels are ordered by (J, E) for plotting.
    print(expt_energies)
    ## level_coordinates = np.array([
    ##     [J*(J+1), E]
    ##     for J, E in expt_energies
    ## ])
    level_coordinates = np.array(expt_energies)
    level_coordinates[:,0] = level_coordinates[:,0]*(level_coordinates[:,0]+1)

    # draw levels
    plot_kw_defaults = dict(
        linestyle="None",
        marker="_",
        color="darkgreen",
        ##markerfacecolor = "darkgreen",
        ##markeredgecolor="darkgreen",
        markersize=7,
        markeredgewidth=2,
    )
    plot_kw_full = {
        **plot_kw_defaults,
        **plot_kw,
    }

    
    ax.plot(level_coordinates[:,0], level_coordinates[:,1], **plot_kw_full)
    

def band_fit(
        results_data,
        network_levels,
        band_members,
        J_values_for_fit=None,
        with_coriolis=False,
        verbose=False
):
    """Obtain band rotational energy parameters, by fitting band energies.

    The results_data is used for level selection, but energies are taken from
    network_levels, to ensure they are with respect to the proper reference
    energy.

    Arguments:

        results_data (mfdnres.ResultsData): Results data

        network_levels (dict): Plotting data for levels, stored as qn->E

        band_members (list of level.Level): Band members (may include some
        beyond those used in fit)

        J_values_for_fit (list of float, optional): J values for band members to use in fit

        with_coriolis (bool, optional): Whether or not to allot Coriolis contribution
        (for K=1/2)

    Returns:
        parameters (tuple of float): band energy parameters (E0,A,a)

    """

    # extract data for fit points
    selected_levels = [
        band_member.select_level(results_data)
        for band_member in band_members
        ]
    if J_values_for_fit is not None:
        selected_levels = [
            qn
            for qn in selected_levels
            if qn[0] in J_values_for_fit
            ]
    J_list = [
        qn[0]
        for qn in selected_levels
    ]
    E_list = [
        network_levels[qn]
        for qn in selected_levels
    ]
    if (verbose):
        print("Band members for fit: {}".format(selected_levels))
        ## print("J: {}".format(J_list))
        print("Energies: {}".format(E_list))

    # construct coefficient matrices
    b = np.array(E_list)
    if with_coriolis:
        A = np.array([
            [1,J*(J+1),(-1)**(J+1/2)*(J+1/2)]
            for J in J_list
        ])
    else:
        A = np.array([
            [1,J*(J+1)]
            for J in J_list
        ])

    # solve system
    parameters = np.linalg.lstsq(A,b,rcond=None)[0]

    # postprocess Coriolis parameter
    if with_coriolis:
        # convert last coefficient to a = c3/c2
        parameters[2] /= parameters[1]
    else:
        # zero pad for parameter a if not already present
        parameters = np.append(parameters,[0],axis=0)
    if (verbose):
        print("Parameters: {}".format(parameters))

    return tuple(parameters)


def draw_band_fit(
        ax,
        results_data,
        network_levels,
        band_members,
        J_values_for_fit=None,
        with_coriolis=False,
        J_range=None,
        plot_kw={},
        verbose=False
):
    """Obtain band rotational energy parameters, by fitting band energies.

    The results_data is used for level selection, but energies are taken from
    network_levels, to ensure they are with respect to the proper reference
    energy.

    Arguments:

        results_data (mfdnres.ResultsData): Results data

        network_levels (dict): Plotting data for levels, stored as qn->E

        band_members (list of level.Level): Band members (may include some
        beyond those used in fit)

        J_values_for_fit (list of float, optional): J values for band members to use in fit

        with_coriolis (bool, optional): Whether or not to allot Coriolis contribution
        (for K=1/2)

        J_range (tuple of float, optional): Range of J values to plot; currently
        must be specified but might someday default to J range for specified
        band members

        plot_kw (dict, optional): Line2D properties properties overriding
        the default line

    """

    if J_range is None:
        raise ValueError("parameter J_range not specified")
    
    # generate fit points
    parameters = band_fit(
        results_data, network_levels, band_members,
        J_values_for_fit=J_values_for_fit,
        with_coriolis=with_coriolis, verbose=verbose,
    )
    plot_points = np.array([
            [
                J*(J+1),
                parameters[0]*1+parameters[1]*J*(J+1)+
                (parameters[2]*(-1)**(J+1/2)*(J+1/2) if with_coriolis else 0)  # suppress noninteger J+1/2
            ]
            for J in np.linspace(*J_range,num=int(J_range[1]-J_range[0])+1)
    ])
    if verbose:
        print("Parameters: {}".format(parameters))
        print("Plot points: {}".format(plot_points))

    # plot curve
    plot_kw_defaults = dict(
        linestyle="-",
        marker=None,
    )
    plot_kw_full = {
        **plot_kw_defaults,
        **plot_kw,
    }
    ax.plot(plot_points[:,0], plot_points[:,1], **plot_kw_full)


def select_network_transitions(
        results_data,
        network_levels,
        operator,
        initial_levels=None,
        reference_strength=None,
        reference_results_data=None,
        strength_threshold=None,
        strength_mode="rtp",
        verbose=False,
):
    """Select transitions for network plot (with plotting data).

    Arguments:

        results_data (mfdnres.ResultsData): Results data

        network_levels (dict): Plotting data for levels, stored as qn->(J*(J+1), Ex)

        initial_levels (list): Levels to take for initial state of (canonical)
        transition, as qn tuples

        reference_strength (float or tuple of level.Level, optional):
        Transition strength, or level providing self-transition strength, or
        tuple of two levels (level_f,level_i) providing transition strength, or None to use
        largest transition strength
        each J

        reference_results_data (mfdnres.ResultsData, optional): Results data for
        use in extracting reference strength
       
        strength_threshold (float, optional): Threshold for transition strength
        (after scaling to reference strength)

        strength_mode (str, optional): Mode ("rtp", "rme", "me") specifying which accessor
        to use for the transition strength

    Returns:

        network_transitions (list): Plotting data for transitions, stored as
        ((qn_f,qn_i),strength)

    """

    # collect transitions
    levels = sorted(list(network_levels.keys()))
    canonical_qn_pairs = [
        (levels[bra_index], levels[ket_index])
        for ket_index in range(len(levels))
        for bra_index in range(ket_index)
    ]
    transition_strengths = []
    if strength_mode == "rme":
        strength_accessor=results_data.get_rme
    elif strength_mode == "rtp":
        strength_accessor=results_data.get_rtp
    elif strength_mode == "me":
        strength_accessor=results_data.get_me
    for qn_pair in canonical_qn_pairs:
        strength = abs(strength_accessor(operator,qn_pair))
        if not np.isnan(strength):
            transition_strengths.append((qn_pair, strength))
    
    # filter and rescale strengths
    if reference_strength is None:
        s0 = max(map(lambda x:x[-1], transition_strengths))
    elif isinstance(reference_strength, tuple):
        if strength_mode == "rme":
            reference_strength_accessor=reference_results_data.get_rme
        elif strength_mode == "rtp":
            reference_strength_accessor=reference_results_data.get_rtp
        elif strength_mode == "me":
            reference_strength_accessor=reference_results_data.get_me
        if reference_results_data is None:
            reference_results_data = results_data
        level_f, level_i = reference_strength
        qn_pair = level_f.select_level(reference_results_data), level_i.select_level(reference_results_data)
        s0 = abs(reference_strength_accessor(operator,qn_pair))
    else:
        s0 = reference_strength
    if verbose:
        print("Reference stregth {}".format(s0))
    network_transitions = [
        (qn_pair, strength/s0)
        for qn_pair, strength in transition_strengths
        if strength_threshold is None or strength/s0 >= strength_threshold
        if initial_levels is None or qn_pair[1] in initial_levels
    ]
        
    if verbose:
        print("Network transitions: {}".format(network_transitions))
    return network_transitions

def draw_network_transitions(
        ax,
        network_levels,
        network_transitions,
        width_scale=1,
        color_scale=1,
        colormap=mpl.cm.get_cmap(name="Greys"),
        plot_kw={},
        verbose=False,
):
    """Draw transitions for network plot.

    Arguments:

        ax (mpl.axes.Axes): axes object

        network_levels (dict): Plotting data for levels, stored as qn->(J*(J+1), E)

        network_transitions (list): Plotting data for transitions, stored as
        ((qn_f,qn_i),strength)

        width_scale (float, optional): Width (in pt) for transition of unit strength

        color_scale (float, optional): Scale factor for strength before applying colormap

        colormap (matplotlib.colors.Colormap, optional): Colormap for color as
        function of scaled transition strength

        plot_kw (dict, optional): Line2D properties properties overriding
        the default line style (but not linewidth or color)

    """

    # sort transitions so weakest is rendered backmost
    sorted_network_transitions = sorted(network_transitions, key = lambda x: x[1])

    plot_kw_defaults = dict(
        solid_capstyle="butt", marker=None,
    )
    plot_kw_full = {
        **plot_kw_defaults,
        **plot_kw,
    }
    
    # draw transitions
    for qn_pair, strength in sorted_network_transitions:
        qnf, qni = qn_pair
        Jf, _, _ = qnf
        Ji, _, _ = qni
        Ef = network_levels[qnf]
        Ei = network_levels[qni]
        ax.plot(
            [jj1(Ji),jj1(Jf)], [Ei,Ef],
            linewidth=width_scale*strength,
            color=colormap(color_scale*strength),
            **plot_kw_full
        )
    
