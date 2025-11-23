"""Analysis and plotting examples.

These plots make use of data from the calculations for Fig. 6 of "bebands":

    M. A. Caprio et al., EPJA 56, 120 (2020),
    https://doi.org/10.1140/epja/s10050-020-00112-0

****************************************************************

WARNING: This file demonstrated the use of the legacy "keyword tuple" notation
for specifying observables for plotting.  New code should use
observable.Observable instead.  This legacy tutorial code is only intended for
backward compatibility tests.

****************************************************************

Contents

    * global plot styling

    * global mesh and styling definitions

    * data input

    * Example: basic single plot ("canned" version)

    * Example: series of basic single plots ("canned" version)

    * Example: multipage survey plot

    * Example: heterogeneous multipanel plot

    * Example: teardrop plot

    * main

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
    - 11/01/21 (mac,zz): Update examples to use newer mfdnres.ticks and mfdnres.multipanel tools.
    - 05/18/22 (mac): Add tick styling to multipanel example.

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
import mfdnres.data
import mfdnres.multipanel
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
# analysis.  Here we make some definitions which define interactions and
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
    #
    # Note: A typical results_data object will contain a params dictionary like the following:
    #
    #   {'Version': 15, 'Revision': 'v15b01-34-gd651299', 'ndiags': 9, 'MPIranks': 45,
    #   'OMPthreads': 64, 'Nprotons': 4, 'Nneutrons': 5, 'TwoMj': 1, 'parity': -1,
    #   'Nmin': 0, 'Nmax': 8, 'DeltaN': 2, 'WTmax': 13.1, 'M': 0.5, 'nuclide': (4, 5),
    #   'A': 9, 'dimension': 63003395, 'numnonzero': 45472165505, 'Hrank': 2, 'hbomeg':
    #   8.0, 'fmass': 938.92, 'TBMEfile': ['tbme-rrel2.bin', 'tbme-Ncm.bin'], 'hw': 8.0,
    #   'tbo_names': ['rrel2', 'Ncm'], 'numTBops': 2, 'run': 'mac0543', 'code_name':
    #   'mfdn15', 'descriptor':
    #   'Z4-N5-Daejeon16-coul1-hw08.000-a_cm50-Nmax08-Mj0.5-lan1500-tol1.0e-04', 'Z': 4,
    #   'N': 5, 'interaction': 'Daejeon16', 'coulomb': 1, 'lawson': 50.0, 'Ncut': None,
    #   'mixed_parity_flag': False, 'fci_flag': False, 'lanczos': 1500, 'tolerance':
    #   '1.0e-04', 'natural_orbital_flag': False, 'natural_orbital_iteration': 0,
    #   'decomposition_operator': None, 'decomposition_lanczos': None,
    #   'decomposition_flag': False, 'subset_index': None, 'extension': 'res',
    #   'natorb_base_state': (None, None, None), 'decomposition_state': (None, None,
    #   None), 'filename':
    #   '/home/mcaprio/results/mcaprio/mfdn/runmac0543/results/res/runmac0543-mfdn15-Z4-N5-Daejeon16-coul1-hw08.000-a_cm50-Nmax08-Mj0.5-lan1500-tol1.0e-04.res'}

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
# Example: basic single plot ("canned" version)
################################################################

def make_basic_plot(mesh_data):
    """ Provides a bare-bones example of making a single analysis plot.

    We take the ground state energy from bebands Fig. 6(a) for this example.

    """

    # plot directory

    # TUTORIAL: Here is where the output files will go.
    
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

    interaction, coulomb = interaction_coulomb
    observable_data = mfdnres.data.make_hw_scan_data(
        mesh_data,nuclide_observable,
        selector =  {"interaction": interaction, "coulomb": coulomb},
        Nmax_range = (NMAX_MIN,Nmax_max), hw_range = hw_range,
        verbose = False
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

    # TUTORIAL: To get you started quickly, this example uses a "canned" routine
    #
    #     mfdnres.data.write_hw_scan_plot()
    #
    # to make a single-panel figure containing the plots for a single
    # observable, and some standardized labels.  This is convenient for a
    # "quick-and-dirty" check.  But often you will want to have more control
    # over your figure, and will need to take over these tasks (generating the
    # axes, adding the plots, adding labels) yourself.  You can look in data.py
    # to see the code for this function.  Its main ingredients are:
    #
    #     # initialize plot
    #     fig, ax = plt.subplots(...)
    #
    #     # provide axis labeling
    #     mfdnres.data.set_up_hw_scan_axes(ax,...)
    #     
    #     # make panel label
    #     mfdnres.data.add_observable_panel_label(ax,...)
    #     
    #     # generate plot
    #     mfdnres.data.add_hw_scan_plot(ax,observable_data,...)
    #     
    #     # finalize plot
    #     plt.savefig(figure_file_name)
    #     plt.close()
    #
    # Then make_survey_plot() below provides an example where we take over these
    # tasks ourselves.
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
# Example: series of basic single plots ("canned" version)
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

            interaction, coulomb = interaction_coulomb

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
# Example: multipage survey plot
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
        for row, qn in enumerate(qn_sublist):
        
            observable_data_by_col = {}
            
            # define observable
            nuclide_observable = ("diff", (nuclide, ("energy", qn)), (nuclide, ("energy", qn_ref)))

            # first pass over columns -- tabulate
            for col, interaction_coulomb in enumerate(INTERACTION_COULOMB_LIST):

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
            for col, interaction_coulomb in enumerate(INTERACTION_COULOMB_LIST):
                observable_data = observable_data_by_col[col]
                observable_range = (min(observable_data.min()["value"],observable_range[0]),max(observable_data.max()["value"],observable_range[1]))

            # third pass over columns -- plot
            for col, interaction_coulomb in enumerate(INTERACTION_COULOMB_LIST):
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
    #
    # TUTORIAL: The basic idea here is that we will loop over all panel indices
    # below, to create the panels.  The code inside the loop is "generic".  It
    # is the same for all panels.  It will know to choose different contents to
    # include in each panel, based on the panel row and column indices, by using
    # the information provided in the following dictionaries.
    interaction_coulomb = INTERACTION_COULOMB_LIST[0]
    nuclide = (4,5)
    Nmax_max = NMAX_MAX_BY_NUCLIDE[nuclide]
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
    observable_tick_specifier_by_panel = {
        (0,0): (-70,-30,5,5),
        (1,0): (-1,6,1,5),
        (0,1): (-1,30,5,5),
        (1,1): (-1,6,1,5),
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

    # initialize figure
    dimensions=(2,2)
    panel_size=(2.2,2.0)
    fig, gs = mfdnres.multipanel.multipanel_fig_gs(
        dimensions=dimensions,
        panel_size=panel_size,
        x_panel_gaps=0.25,
        y_panel_gaps=0.02,
    )
    
    for panel_indices in mfdnres.multipanel.grid_iterator(dimensions):

        nuclide_observable_list = nuclide_observable_list_by_panel[panel_indices]

        # find range parameters
        hw_range = HW_RANGE_BY_INTERACTION_COULOMB[interaction_coulomb]
        observable_range = observable_range_by_panel[panel_indices]

        # construct axes
        ax = fig.add_subplot(gs[panel_indices[0], panel_indices[1]])
        
        # set manual ticks

        # TUTORIAL: You can see that matplotlib's default hw ticks, for the plot
        # range we set below, is atrocious (major tick marks at 10 and 20, with
        # no minor ticks at allT).  Compare bebands Fig. 6, where major ticks
        # are in steps of 5, with minor ticks in steps of 2.5 (i.e., 2
        # subintervals).  In mfdnres.ticks, we find some functions to provide
        # manual control over tick mark intervals, a la my Mathematica package
        # CustomTicks.  To override the default ticks on the hw axis, enable the
        # following code...

        if True:
            x_ticks = mfdnres.ticks.linear_ticks(0,50,5,2)
            mfdnres.ticks.set_ticks(ax,"x",x_ticks)
            y_ticks = mfdnres.ticks.linear_ticks(*observable_tick_specifier_by_panel[panel_indices])
            mfdnres.ticks.set_ticks(ax,"y",y_ticks)

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
        mfdnres.multipanel.suppress_interior_labels(ax,axis="x")

        # format ticks
        ax.tick_params(which="both", top=True, right=True, labelsize="x-small")
        
        # panel letter
        ax.annotate(
            mfdnres.multipanel.panel_letter(dimensions,panel_indices,direction="vertical"),
            xy=(0.05,0.95),xycoords="axes fraction",
            horizontalalignment="left",
            verticalalignment="top",
            fontsize="small",
        )

        # panel label
        ax.annotate(
            panel_label_text_by_panel[panel_indices],
            xy=(0.90,0.90),xycoords="axes fraction",
            multialignment="left",
            horizontalalignment="right",
            verticalalignment="top",
            ##bbox=dict(boxstyle="round",facecolor="white"),
        )

        # Nmax label (on just one panel)
        if panel_indices==(0,0):
            ax.annotate(
                r"${}$".format(mfdnres.data.Nmax_label_text(Nmax_max)),
                xy=(0.05,0.035),
                xycoords=("axes fraction","axes fraction"),
                multialignment="left",
                horizontalalignment="left",
                verticalalignment="bottom",
                fontsize="x-small",
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
                Nmax_range = (NMAX_MIN,Nmax_max), hw_range = hw_range,
            )
        
            # write data
            mfdnres.data.write_hw_scan_data(
                descriptor,observable_data,
                directory=plot_directory,
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
                Nmax_plot_style_kw = dict(marker_size=4, Nmax_dashing=(lambda Nmax_relative : dashing)),
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
# Example: teardrop plot
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
    dimensions = (1,len(nuclide_observable_list))
    fig, gs = mfdnres.multipanel.multipanel_fig_gs(dimensions=dimensions,panel_size=(0.75,3.))

    # tabulate observables
    for panel_indices in mfdnres.multipanel.grid_iterator(dimensions):
        observable_index = mfdnres.multipanel.panel_index(dimensions,panel_indices)
        nuclide_observable = nuclide_observable_list[observable_index]

        # construct axes
        row, col = panel_indices
        ax = fig.add_subplot(gs[row, col])
        
        ax.set_xlim(0,1)
        ax.set_xticks([])
        ax.set_axisbelow(b=True)
        
        ax.set_ylabel(r"${}$".format(mfdnres.data.make_observable_axis_label_text(nuclide_observable)))
        ax.set_ylim(*mfdnres.data.extend_interval_relative(observable_range,observable_range_extension))
        ax.grid(axis="y",linewidth=0.5,linestyle=":",color="gray")
        mfdnres.multipanel.suppress_interior_labels(ax)
        ## if not ax.is_first_col():
        ##     ax.tick_params(axis="y", length=0)

        # set manual ticks
        y_ticks = mfdnres.ticks.linear_ticks(-2,2,0.5,5)
        mfdnres.ticks.set_ticks(ax,"y",y_ticks)
            
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
        
        for interaction_index, interaction_coulomb in enumerate(INTERACTION_COULOMB_LIST):

            interaction, coulomb = interaction_coulomb

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
    make_plot_series(mesh_data)
    make_survey_plot(mesh_data)
    make_multipanel_plot(mesh_data)
    make_teardrop_plot(mesh_data)

if __name__ == "__main__":
    main()
