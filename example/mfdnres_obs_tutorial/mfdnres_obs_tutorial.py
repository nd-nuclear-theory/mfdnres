"""mfdnres_obs_tutorial.py

Analysis and plotting examples, drawing on the data for Fig. 6 of "bebands":

    M. A. Caprio et al., EPJA 56, 120 (2020),
    https://doi.org/10.1140/epja/s10050-020-00112-0

Note for user

    Comments starting with "TUTORIAL" are commentary meant for your reading
    pleasure in this tutorial.  They also flag optional, illustrative "print"
    statements which you can turn on to get a better feel for what the data look
    like.  They are not meant as comments for actual production code, where they
    would just be clutter.  So, if you use this example code as a template for
    your analysis code, please, please, please, delete the TUTORIAL comments.
    Don't be that person whose LaTeX manuscript files all look like the RevTeX
    sample file, complete with the chatty "helpful hint" comments...

Setup

    Make sure you have downloaded the results from the following runs:

        "mac0455"
        "mac0468"
        "mac0543"

    To do so:

    1) Set up your local results directory.  Although you can set directory
    names as you wish, our default tree, assumed below, mirrors the structure
    shown above for m2032:
    
       ${GROUP_HOME}/results/<user>/<code>
    
    So, under your personal account, I recommend setting the environment variable
    
       setenv GROUP_HOME $HOME
       mkdir -p ${GROUP_HOME}/results/mcaprio/mfdn
    
    2) Download the results file archives into your local
    ${GROUP_HOME}/results/mcaprio/mfdn.  These can be found as
    
       /global/cfs/cdirs/m2032/results/mcaprio/mfdn/run<run>-archive-<date>-res.tgz
    
    or, for older runs,
    
       /global/cfs/cdirs/m2032/results/mcaprio/mfdn/run<run>-archive-<date>.tgz
    
    3) Untar them!  E.g.,
    
       tar xvf runmac0455-archive-210214-res.tgz

  Mark A. Caprio
  Department of Physics
  University of Notre Dame

    - 04/07/21 (mac): Created, drawing on tabulate_obs_8li-coulex.py and emratio_obs.py.
    - 04/08/21 (mac): Clean up mfdnres.data interface.  Add survey and multipanel examples.
    - 04/10/21 (mac): Use os.makedirs.  Add experimental use of mfdnres.ticks.

"""

import itertools
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pgf
import numpy as np
import pandas as pd

import mfdnres
import mfdnres.ticks
import mfdnres.data  # TODO 04/08/21 (mac): Rename mfdnres.data, once new name decided...
import mfdnres.ncci

################################################################
# global plot styling
################################################################

# TUTORIAL: Define some sensible plot styling defaults for matplotlib...

mpl.rcParams.update(mfdnres.data.SENSIBLE_PLOT_STYLE)

################################################################
# global mesh and styling definitions
################################################################

# TUTORIAL: While normally global definitions are to be avoided, sometimes
# it is desirable to have a uniform set of definitions throughout an
# analysis.  Here we make some definitions which define itneractions and
# parameter meshes used throughout the analysis examples below.  To avoid
# confusion with local variables, we will follow Google Python style
# convention and use all caps.
    
# list of interactions to use: interaction name and Coulomb on/off flag
INTERACTION_COULOMB_LIST = [
    ("Daejeon16",1),
    ("JISP16",1),
    ("LENPICchi2bi2C",0)
]

# hw range for each interaction
HW_RANGE_BY_INTERACTION_COULOMB = {
    ("Daejeon16",1): (5.,25.),
    ("JISP16",1): (10.,30.),
    ("LENPICchi2bi2C",0): (15.,35.),
}

# special preferred hw for each interaction
HW_BY_INTERACTION_COULOMB = {
    ("Daejeon16",1): 15.,
    ("JISP16",1): 20.,
    ("LENPICchi2bi2C",0): 25.,
}

# matplotlib marker style for each interaction
MARKER_BY_INTERACTION_COULOMB = {
    ("Daejeon16",1): "o",
    ("JISP16",1): "s",
    ("LENPICchi2bi2C",0): "D",
}

# lower Nmax cutoff
NMAX_MIN = 4

# highest Nmax to use for each nuclide
NMAX_MAX_BY_NUCLIDE = {
    (4,5): 10
}

################################################################
# data input
################################################################

def read_data():
    """"""

    # define run list for input

    # TUTORIAL: Here we build a list of directories from which to "slurp" up
    # data files.  In this example, we read all data files ("*.res"), but we can
    # be more selective using glob patterns on the filenames.
    
    run_list = [

        # 9Be (4,5)
        "mac0455",  # 9Be-
        "mac0468",  # 9Be- SUPPLEMENT
        "mac0543",  # 9Be- (+TRANS) SUPPLEMENT Daejeon16 integer mesh at low hw
        
        # 9B (5,4)
        ## "mac0477",  # 9B-
        ##"Z5.N4",
    ]
    data_dir_list = [
        mfdnres.res.res_file_directory("mcaprio","mfdn",run)
        for run in run_list
    ]

    # TUTORIAL: To see what directory names we generated, enable the following...
    if False:
        print(data_dir_list)
        exit()
    
    # read runs
    res_format = "mfdn_v15"
    filename_format="mfdn_format_7_ho"
    mesh_data = mfdnres.res.slurp_res_files(
        data_dir_list,res_format,filename_format,
        glob_pattern="*.res",
        verbose=True
    )

    # TUTORIAL: To see what our initial, raw mesh of MFDnResultsData objects looks like...
    if False:
        print(mesh_data)
    if False:
        for results_data in mesh_data:
            print(results_data.params)
        exit()
    
    # merge results data

    # TUTORIAL: Here is where we merge data from from different calculations at
    # the same "mesh point", e.g., calculations with different M, or from mfdn
    # and obscalc-ob.
    
    mesh_data = mfdnres.analysis.merged_mesh(
        mesh_data,
        ("nuclide","interaction","coulomb","hw","Nmax"),
        postprocessor=mfdnres.ncci.augment_params_with_parity,
        verbose=False
    )

    # filter for cases of interest

    # TUTORIAL: Here we can optionally prune out data we want to completely
    # exclude from the analysis.  E.g., in this example, we exclude data
    # calculated at the "nonstandard" mesh point hw=7.5.
    
    mesh_data = [
        results_data
        for results_data in mesh_data
        if results_data.params["hw"]!=7.5
    ]
    
    # TUTORIAL: To see what our mesh of MFDnResultsData objects looks like
    # now...
    if False:
        print("Merged mesh (keys)")
        mfdnres.analysis.mesh_key_listing(
            mesh_data,
            ("nuclide","interaction","coulomb","hw","Nmax","parity"),
            verbose=True
        )
        exit()

    return mesh_data

################################################################
# Example: basic single plot
################################################################

def make_basic_plot(mesh_data):
    """ Provides a bare-bones example of making a single analysis plot.

    We take the ground state energy from bebands Fig. 6(a) for this example.

    """

    # TUTORIAL: Here is where the output will go.
    
    plot_directory="plots/basic"
    os.makedirs(plot_directory, exist_ok=True)

    # mesh parameters
    interaction_coulomb = INTERACTION_COULOMB_LIST[0]  # ("Daejeon16",1)
    hw_range = HW_RANGE_BY_INTERACTION_COULOMB[interaction_coulomb]  # (5.,25.)
    nuclide = (4,5)
    Nmax_max = NMAX_MAX_BY_NUCLIDE[nuclide]  # 10

    # define observable

    # TUTORIAL: It turns out to be most convenient to group our nuclide and
    # observble together in a tuple.  We usually need to know both together.
    # See the docstring to data.make_hw_scan_data for a summary of the syntax
    # for observables.
    
    nuclide_observable = ((4,5), ("energy", (1.5,1,1)))
        
    # generate descriptor

    # TUTORIAL: The "descriptor" is a standard string representation of the
    # observable, which we will use in output filenames.
    
    descriptor = mfdnres.data.hw_scan_descriptor(interaction_coulomb,nuclide_observable)

    # TUTORIAL: To see the descriptors...
    if False:
        print(descriptor)

    # tabulate

    # TUTORIAL: Here we digest the observable data into a pandas.DataFrame
    # object.  This is basically a "spreadsheet" with just one column,
    # containing the observable values.  The (Nmax,hw) for each value is given
    # as a compound row label (technically, a pandas MultiIndex).

    # TUTORIAL: The tabulation routine mfdnres.data.make_hw_scan_data will take
    # care of selecting the nuclide. But our mesh contains results from many
    # interactions.  We therefore have to explicitly filter for the interaction
    # of interest, using the optional argument
    #
    #    selector =  {"interaction": interaction, "coulomb": coulomb}
    #
    # This is equivalent to manually filtering the mesh:
    #
    #     mesh_data_selected = mfdnres.analysis.selected_mesh_data(
    #         mesh_data,
    #         {"interaction":interaction,"coulomb":coulomb}
    #     )

    # TUTORIAL: In general, we may want to prune to a given window in hw and
    # Nmax.  E.g., this gets rid of stray low-Nmax or high-Nmax test runs, or
    # extremes in hw beyond the desired plotting range.  We could do this cut
    # "before" tabulation, by triming the data mesh.  However, instead, we opt
    # to do this cut "afterwards", by slicing the pandas DataFrame on its
    # (Nmax,hw) MultiIndex.  The ranges are specified using the optional
    # arguments
    #
    #     Nmax_range, hw_range
    #
    # This is equivalent to
    #
    #     observable_data = observable_data.loc[(slice(*Nmax_range),slice(*hw_range)),:]

    # TUTORIAL: To see selected mesh and final sliced DataFrame, set
    # verbose=True.
    
    interaction_coulomb = INTERACTION_COULOMB_LIST[0]  # ("Daejeon16",1)
    (interaction,coulomb) = interaction_coulomb
    observable_data = mfdnres.data.make_hw_scan_data(
        mesh_data,nuclide_observable,
        selector =  {"interaction": interaction, "coulomb": coulomb},
        Nmax_range = (NMAX_MIN,Nmax_max), hw_range = hw_range,
        verbose = True
        )

    # write data

    # TUTORIAL: Right here and now we write out the tabular data.  This way we
    # can go back and check numerical values, and we have it for posterity, to
    # submit as suppmemental data with a publication, re-plot it years from now
    # in other software, etc.
    
    mfdnres.data.write_hw_scan_data(
        descriptor,observable_data,
        directory=plot_directory
    )

    # make plot

    # TUTORIAL: This is the "canned" routine to set up a single-panel figure,
    # label it, and plot the data.  We will see later how to break out the
    # different parts of the task, if we want more control.
    #
    # Also, here we manually set the plot range for the observable.  We'll see
    # ways to automate this later.  If observable_range=None, then matplotlib
    # will be allowed to pick the plot range according to its defaults.
    
    observable_range = (-60.,-40.)
    ## print(observable_range)
    mfdnres.data.write_hw_scan_plot(
        descriptor,
        interaction_coulomb,nuclide_observable,
        observable_data,
        hw_range=hw_range,
        observable_range=observable_range,
        Nmax_max=Nmax_max,
        directory=plot_directory,
        verbose=True
    )

################################################################
# Example: series of individual plots
################################################################

def make_plot_series(mesh_data):
    """ Make a series of standalone plots for analysis purposes.

    We take the observables from bebands Fig. 6 for this example, but with
    mutltiple interactions.
    """

    plot_directory="plots/series"
    os.makedirs(plot_directory, exist_ok=True)

    # observable definitions

    # TUTORIAL: This was originally just a list of (nuclide,observable) pairs.
    # But we've turned it into a dictionary, which we use to manually specify a
    # plot range, as well.
    #
    # Note the use of "compound" observables, automatically constructed as
    # differences or ratios of other (simple) observables.  The descriptor
    # string, plot labels, etc., will also be automatically constructed out of
    # those for the individual simple observables, as well.

    nuclide = (4,5)  # for choosing Nmax_max
    nuclide_observable_list = {

        # energies
        ((4,5), ("energy", (1.5,1,1))): (-60,-40),
        ## ((4,5), ("energy", (2.5,1,1))): (-60,-40),

        # excitation energy
        ("diff", ((4,5), ("energy", (2.5,1,1))), ((4,5), ("energy", (1.5,1,1)))): (0.,5.),

        # B(E2)s
        ((4,5), ("rtp", "E2p", (1.5,1,1), (2.5,1,1))): (0.,30.),
        ## ((4,5), ("rtp", "E2p", (1.5,1,1), (3.5,1,1))): (0.,30.),

        # B(E2) ratio
        ("ratio", ((4,5), ("rtp", "E2p", (1.5,1,1), (2.5,1,1))), ((4,5), ("rtp", "E2p", (1.5,1,1), (3.5,1,1)))): (0.,5.),
    }

    # tabulate and plot
    for interaction_coulomb in INTERACTION_COULOMB_LIST:
        for nuclide_observable in nuclide_observable_list:

            (interaction,coulomb) = interaction_coulomb

            # generate descriptor
            descriptor=mfdnres.data.hw_scan_descriptor(interaction_coulomb,nuclide_observable)

            # tabulate
            hw_range = HW_RANGE_BY_INTERACTION_COULOMB[interaction_coulomb]
            Nmax_max = NMAX_MAX_BY_NUCLIDE[nuclide]
            observable_data = mfdnres.data.make_hw_scan_data(
                mesh_data,nuclide_observable,
                selector =  {"interaction": interaction, "coulomb": coulomb},
                Nmax_range = (NMAX_MIN,Nmax_max), hw_range = hw_range
            )
            
            # TUTORIAL: To see the pandas data frames...
            if False:
                print(observable_data)
        
            # write data
            mfdnres.data.write_hw_scan_data(
                descriptor,observable_data,
                directory=plot_directory
            )

            # make plot

            observable_range = nuclide_observable_list[nuclide_observable]

            # TUTORIAL: Alternatively, if you wanted to to choose a range
            # including zero and the extremes of data, ...
            ## observable_range = (min(observable_data.min()["value"],0),max(observable_data.max()["value"],0))
            
            mfdnres.data.write_hw_scan_plot(
                descriptor,
                interaction_coulomb,nuclide_observable,
                observable_data,
                hw_range=hw_range,
                observable_range=observable_range,
                Nmax_max=Nmax_max,
                directory=plot_directory,
                verbose=True
            )

################################################################
# Example: make multipage survey plot
################################################################

def make_survey_plot(mesh_data):
    """Make a multipanel survey of observables for multiple interactions.

    We take a sampling of excitation energies in the same nuclide as our
    observable set, but these observables could also typically represent a
    survey across multiple nuclides.

    """

    plot_directory="plots/survey"
    os.makedirs(plot_directory, exist_ok=True)

    # figure layout parameters
    num_rows = 3
    num_cols = len(INTERACTION_COULOMB_LIST)
    figsize = np.array((3,2))*np.array((num_cols,num_rows))

    # figure contents
    nuclide = (4,5)
    qn_ref = (1.5,1,1)  # reference state for Ex
    qn_list = [
        # yrast band
        ## (1.5,1,1),
        (2.5,1,1),
        (3.5,1,1),
        (4.5,1,1),
        # yrare band
        (0.5,1,1),
        (1.5,1,2),
        (2.5,1,2),
        (3.5,1,2),
    ]
        
    # initialize multipage pdf file
    pdf_file_name = os.path.join(
        plot_directory,
        "survey.pdf"
    )
    print(pdf_file_name)
    pdf = mpl.backends.backend_pgf.PdfPages(pdf_file_name)

    # for each page of observables
    for qn_sublist in mfdnres.data.partitions(qn_list,num_rows):

        # initialize figure for page
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(nrows=num_rows, ncols=num_cols, hspace=0., wspace=0.)

        # for each row
        for (row,qn) in enumerate(qn_sublist):
        
            observable_data_by_col = {}
            
            # define observable
            nuclide_observable = ("diff", (nuclide, ("energy", qn)), (nuclide, ("energy", qn_ref)))

            # first pass over columns -- tabulate
            for (col,interaction_coulomb) in enumerate(INTERACTION_COULOMB_LIST):

                (interaction,coulomb) = interaction_coulomb

                # generate descriptor
                descriptor=mfdnres.data.hw_scan_descriptor(interaction_coulomb,nuclide_observable)

                # tabulate
                hw_range = HW_RANGE_BY_INTERACTION_COULOMB[interaction_coulomb]
                Nmax_max = NMAX_MAX_BY_NUCLIDE[nuclide]
                observable_data = mfdnres.data.make_hw_scan_data(
                    mesh_data,nuclide_observable,
                    selector =  {"interaction": interaction, "coulomb": coulomb},
                    Nmax_range = (NMAX_MIN,Nmax_max), hw_range = hw_range
                )
                
                # write data
                mfdnres.data.write_hw_scan_data(
                    descriptor,observable_data,
                    directory=plot_directory
                )

                # save data for plotting
                observable_data_by_col[col] = observable_data
                
            # second pass over columns -- find common y-axis range
            observable_range = (0,0)
            for (col,interaction_coulomb) in enumerate(INTERACTION_COULOMB_LIST):
                observable_data = observable_data_by_col[col]
                observable_range = (min(observable_data.min()["value"],observable_range[0]),max(observable_data.max()["value"],observable_range[1]))

            # third pass over columns -- plot
            for (col,interaction_coulomb) in enumerate(INTERACTION_COULOMB_LIST):
                observable_data = observable_data_by_col[col]

                # construct axes
                ax = fig.add_subplot(gs[row, col])
            
                # draw axes
                hw_range = HW_RANGE_BY_INTERACTION_COULOMB[interaction_coulomb]
                mfdnres.data.set_up_hw_scan_axes(
                    ax,
                    nuclide_observable,
                    hw_range,
                    observable_range,
                    observable_range_extension=(0.05,0.10)
                )

                # eliminate labels from interior panel edges
                if not ax.is_last_row():
                    ax.set_xlabel(None)
                    ax.set_xticklabels([])
                if not ax.is_first_col():
                    ax.set_ylabel(None)
                    ax.set_yticklabels([])

                # make panel label
                mfdnres.data.add_observable_panel_label(
                    ax,
                    interaction_coulomb,
                    nuclide_observable
                )

                # make plot
                mfdnres.data.add_hw_scan_plot(ax,observable_data,Nmax_max)
                
        # finalize figure for page
        pdf.savefig()
        plt.close()
                    
    # finalize multipage pdf
    pdf.close()
    
################################################################
# Example: heterogeneous multipanel plot
################################################################

def make_multipanel_plot(mesh_data):
    """Make a custom multipanel plot (with some observables overlaid).

    We roughly reproduce bebands Fig. 6 for this example, sans some annotations.

    """

    plot_directory="plots/multipanel"
    os.makedirs(plot_directory, exist_ok=True)

    # plot contents
    interaction_coulomb = INTERACTION_COULOMB_LIST[0]
    nuclide = (4,5)
    Nmax_max = NMAX_MAX_BY_NUCLIDE[nuclide]
    panel_letter_by_panel = {
        # TODO 04/08/21 (mac): Implement generic panel letter generator function
        # a la SciDraw Multipanel.
        (0,0): "a",
        (1,0): "b",
        (0,1): "c",
        (1,1): "d",
        }
    panel_label_text_by_panel = {
        (0,0): "$E$",
        (1,0): "$\Delta E$",
        (0,1): "$B(E2)$",
        (1,1): "$B(E2)$ ratio",
        }
    observable_range_by_panel = {
        (0,0): (-60.,-43.),
        (1,0): (0.,5.),
        (0,1): (0.,25.),
        (1,1): (0.,5.),
        }
    nuclide_observable_list_by_panel = {

        (0,0): [
            # energies
            ((4,5), ("energy", (1.5,1,1))),
            ((4,5), ("energy", (2.5,1,1))),
        ],

        (1,0): [
        # excitation energy
        ("diff", ((4,5), ("energy", (2.5,1,1))), ((4,5), ("energy", (1.5,1,1)))),
        ],

        (0,1): [
            # B(E2)s
            ((4,5), ("rtp", "E2p", (1.5,1,1), (2.5,1,1))),
            ((4,5), ("rtp", "E2p", (1.5,1,1), (3.5,1,1))),
        ],

        (1,1): [
            # B(E2) ratio
            (
                "ratio",
                ((4,5), ("rtp", "E2p", (1.5,1,1), (2.5,1,1))),
                ((4,5), ("rtp", "E2p", (1.5,1,1), (3.5,1,1)))
            ),
        ],
    }

    # figure layout parameters
    panel_size = (3.,2.)
    dimensions = (2,2)
    figsize = np.array(panel_size)*np.array(dimensions)

    # initialize figure
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(*dimensions, hspace=0., wspace=0.25)

    for (row, col) in itertools.product(range(2),range(2)):

        nuclide_observable_list = nuclide_observable_list_by_panel[(row, col)]

        # find range parameters
        hw_range = HW_RANGE_BY_INTERACTION_COULOMB[interaction_coulomb]
        observable_range = observable_range_by_panel[(row, col)]

        # construct axes
        ax = fig.add_subplot(gs[row, col])
        
        # set manual ticks (experimental interface)

        # TUTORIAL: I am in the process of implementing manual control over tick
        # mark intervals, a la CustomTicks.  The interface is still under
        # development.  But you can see that matplotlib's default hw ticks, for
        # the plot range we set below, is atrocious (major tick marks at 10 and
        # 20, with no minor ticks at all).  Compare bebands Fig. 6, where major
        # ticks are in steps of 5, with minor ticks in steps of 2.5 (i.e., 2
        # subintervals).  To override the default ticks on the hw axis, enable
        # the following (experimental) code...

        if False:
            major_ticks, minor_ticks = mfdnres.ticks.linear_tick_locations(0,50,5,2)
            ax.set_xticks(major_ticks)
            ax.set_xticks(minor_ticks, minor=True)

        # draw axes

        # TUTORIAL: All observables in each panel are of the same type (e.g.,
        # all are energies), so we can simply pick the first observable in the list, to use as a
        # "representative" observable for defining the y axis label.

        # TUTORIAL: We also use the optional arguments
        #
        #     hw_range_extension=(0.05,0.05),
        #     observable_range_extension=(0.05,0.05),
        #
        # to set bigger (5% or 10%) margins on the specified x and y axis ranges.

        mfdnres.data.set_up_hw_scan_axes(
            ax,
            nuclide_observable_list[0],
            hw_range,
            observable_range,
            hw_range_extension=(0.15,0.15),
            observable_range_extension=(0.05,0.05),
        )
        
        # eliminate labels from interior panel edges
        if not ax.is_last_row():
            ##ax.get_xaxis().get_label().set_visible(False)
            ax.set_xlabel(None)
            ax.set_xticklabels([])

        # panel letter
        ax.annotate(
            "({})".format(panel_letter_by_panel[(row, col)]),
            xy=(0.05,0.95),xycoords="axes fraction",
            horizontalalignment="left",
            verticalalignment="top",
            fontsize="small"
        )

        # panel label
        ax.annotate(
            panel_label_text_by_panel[(row, col)],
            xy=(0.90,0.90),xycoords="axes fraction",
            multialignment="left",
            horizontalalignment="right",
            verticalalignment="top",
            bbox=dict(boxstyle="round",facecolor="white")
        )
            
        # tabulate and plot each observable
        for plot_index, nuclide_observable in enumerate(nuclide_observable_list):

            # generate descriptor
            descriptor=mfdnres.data.hw_scan_descriptor(interaction_coulomb,nuclide_observable)

            # tabulate
            (interaction,coulomb) = interaction_coulomb
            observable_data = mfdnres.data.make_hw_scan_data(
                mesh_data,nuclide_observable,
                selector =  {"interaction": interaction, "coulomb": coulomb},
                Nmax_range = (NMAX_MIN,Nmax_max), hw_range = hw_range
            )
        
            # write data
            mfdnres.data.write_hw_scan_data(
                descriptor,observable_data,
                directory=plot_directory
            )

            # add plot

            # TUTORIAL: To override the default styling for the curves, we use
            # the optional argument Nmax_plot_style_kw.  This provides "hooks"
            # to control various plot attributes as functions of Nmax (see
            # mfdnres.data.Nmax_plot_style).  In particular, to replicate the
            # bebands figure, we override the default dashing, so that it no
            # longer reflects Nmax, but rather is used to distinguish plots of
            # different observables.
            
            if plot_index==0:
                dashing = (None,None)
            else:
                dashing = (1,1)
            mfdnres.data.add_hw_scan_plot(
                ax,
                observable_data,
                Nmax_max=Nmax_max,
                Nmax_plot_style_kw = dict(marker_size=4, Nmax_dashing=(lambda Nmax_relative : dashing))
            )
            
    # finalize plot
    figure_file_name = os.path.join(
        plot_directory,
        "multipanel.pdf"
        )
    print(figure_file_name)
    plt.savefig(figure_file_name)
    plt.close()

################################################################
# teardrop plot
################################################################

EXPT_M1_MOMENT_BY_NUCLIDE = {
    # M1 moment (with uncertainty); from Stone 2005
    (4,5): (-1.177,None),
    }

def make_teardrop_plot(mesh_data):
    """ Generate "teardrop" plot showing convergence of observables with Nmax at fixed hw.
    """

    # plotting parameters
    plot_directory="plots/teardrop"
    os.makedirs(plot_directory, exist_ok=True)

    # plot contents
    nuclide = (4,5)
    nuclide_observable_list = [
        ((4,5), ("moment", "M1", (1.5,1,1))),
        ((4,5), ("moment", "Dlp", (1.5,1,1))),
        ((4,5), ("moment", "Dln", (1.5,1,1))),
        ((4,5), ("moment", "Dsp", (1.5,1,1))),
        ((4,5), ("moment", "Dsn", (1.5,1,1))),
        ]

    # initialize plot
    figsize = (6,4)
    observable_range = (-1.75,1)
    observable_range_extension = (0.02,0.02)
    marker_size = 6
    num_panels = len(nuclide_observable_list)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=1, ncols=num_panels, hspace=0., wspace=0.)

    # tabulate observables
    for (observable_index,nuclide_observable) in enumerate(nuclide_observable_list):

        # construct axes
        row = 0
        col = observable_index
        ax = fig.add_subplot(gs[row, col])
        ax.set_xlim(0,1)
        ax.set_xticks([])
        ax.set_axisbelow(b=True)
        ax.set_ylabel(mfdnres.data.make_observable_axis_label_text(nuclide_observable))
        ax.set_ylim(*mfdnres.data.extend_interval_relative(observable_range,observable_range_extension))
        ax.grid(axis="y",linewidth=0.5,linestyle=":",color="gray")
        mfdnres.data.suppress_interior_labels(ax)
        if not ax.is_first_col():
            ax.tick_params(axis="y", length=0)
        
        # label for observable
        ax.annotate(
            r"${}$".format(mfdnres.data.make_observable_text(nuclide_observable)),
            xy=(0.5,0),
            xycoords=("axes fraction","axes fraction"),
            multialignment="left",
            horizontalalignment="center",
            verticalalignment="bottom"
        )

        # experimental value marker
        if observable_index==0:
            mfdnres.data.add_expt_marker_band(ax,(0.1,0.9),EXPT_M1_MOMENT_BY_NUCLIDE[nuclide])
        
        for (interaction_index,interaction_coulomb) in enumerate(INTERACTION_COULOMB_LIST):

            (interaction,coulomb) = interaction_coulomb

            hw = HW_BY_INTERACTION_COULOMB[interaction_coulomb]
            Nmax_max = NMAX_MAX_BY_NUCLIDE[nuclide]
            observable_data = mfdnres.data.make_hw_scan_data(
                mesh_data,nuclide_observable,
                selector =  {"interaction": interaction, "coulomb": coulomb},
                Nmax_range = (NMAX_MIN,Nmax_max), hw_range = (hw,hw)
            )
            
            # plot data
            plot_data = observable_data.reset_index()
            plot_data["x"] = 0.3+0.2*interaction_index
            plot_data["s"] = (marker_size*((plot_data["Nmax"]-Nmax_max).transform(mfdnres.data.Nmax_symbol_scale)))**2
            plot_data["c"] = (plot_data["Nmax"]-Nmax_max).transform(mfdnres.data.Nmax_color)
            ##print(plot_data)
            ax.plot(
                "x","value",
                data=plot_data,
                marker=None,
                zorder=2.1
            )
            ax.scatter(
                "x","value",
                s="s",c="c",
                ##c=plot_data["c"],s=plot_data["s"],
                data=plot_data,
                edgecolors="black",
                marker=MARKER_BY_INTERACTION_COULOMB[interaction_coulomb],
                zorder=2.2
            )


    # finalize plot
    figure_file_name = os.path.join(
        plot_directory,
        "teardrop.pdf"
        )
    print(figure_file_name)
    plt.savefig(figure_file_name)
    plt.close()
    
################################################################
# main
################################################################

def main():    

    mesh_data=read_data()

    make_basic_plot(mesh_data)
    ## make_plot_series(mesh_data)
    ## make_survey_plot(mesh_data)
    ## make_multipanel_plot(mesh_data)
    ## make_teardrop_plot(mesh_data)

if __name__ == "__main__":
    main()