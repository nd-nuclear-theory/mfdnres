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

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mfdnres
import mfdnres.data
import mfdnres.ncci

################################################################
# global plot styling
################################################################

# TUTORIAL: Define some sensible defaults for matplotlib...

mpl.rcParams["font.family"] = "serif"
mpl.rcParams["mathtext.fontset"] = "dejavuserif"
mpl.rcParams["xtick.labelsize"] = "small"
mpl.rcParams["ytick.labelsize"] = "small"
mpl.rcParams["lines.linewidth"] = 1
mpl.rcParams["axes.prop_cycle"] = mpl.cycler(color=["black"])
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["xtick.major.pad"] = 1
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["ytick.major.pad"] = 1

################################################################
# control code
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

    # select cases of interest

    # TUTORIAL: Here we can optionally prune out data we don't want to include,
    # e.g., here we exclude any data with Nmax<4 and data calculated at the
    # "nonstandard" mesh point hw=7.5.
    
    mesh_data = [
        results_data
        for results_data in mesh_data
        if results_data.params["Nmax"]>=4 and results_data.params["hw"]!=7.5
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

def make_tabulations_be2_ratio(mesh_data):
    """ Generate ratios of "upward" transition B(E2) to "moment" B(E2).

    Legacy: tabulate_obs -> ObservablePlots
    """

    # content lists
    interaction_coulomb_list = [("Daejeon16",1),("JISP16",1),("LENPICchi2bi2C",0)]
    nuclide_transition_pair_list = [
        
        # 7Be/7Li transition to moment ratio
        (((3,4),("E2p",(0.5,1,1),(1.5,1,1))),((3,4),("E2p",(1.5,1,1),(1.5,1,1)))),
        (((4,3),("E2p",(0.5,1,1),(1.5,1,1))),((4,3),("E2p",(1.5,1,1),(1.5,1,1)))),

        # 8Li transition to moment ratio
        (((3,5),("E2p",(1.0,0,1),(2.0,0,1))),((3,5),("E2p",(2.0,0,1),(2.0,0,1)))),
        ##(((3,5),("E2p",(3.0,0,1),(2.0,0,1))),((3,5),("E2p",(2.0,0,1),(2.0,0,1)))),

        # 9Be transition to moment ratio
        (((4,5),("E2p",(2.5,1,1),(1.5,1,1))),((4,5),("E2p",(1.5,1,1),(1.5,1,1)))),
        (((4,5),("E2p",(3.5,1,1),(1.5,1,1))),((4,5),("E2p",(1.5,1,1),(1.5,1,1)))),
        
    ]

    # tabulate absolute transitions
    nuclide_transition_set = {
        nuclide_transition_pair[member]
        for nuclide_transition_pair in nuclide_transition_pair_list for member in [0,1]
    }
    for interaction_coulomb in interaction_coulomb_list:
        for nuclide_transition in sorted(list(nuclide_transition_set)):
            tabulate_obs.make_rtp_table_file(
                mesh_data,interaction_coulomb,nuclide_transition,
                directory="data/8li-coulex/rtp",
                verbose=False)

    # tabulate ratios
    for interaction_coulomb in interaction_coulomb_list:
        for nuclide_transition_pair in nuclide_transition_pair_list:
            tabulate_obs.make_rtp_ratio_table_file(
                mesh_data,interaction_coulomb,nuclide_transition_pair,
                directory="data/8li-coulex/rtp-ratio",
                verbose=False
            )

def make_plots(mesh_data):
    """ Generate plots for 8Li moment and transition analysis.

    pdftk plot-*-moment-M1-02.0-0-01.pdf plot-*-rtp-M1-02.0-0-01-01.0-0-01.pdf plot-*-rtp-M1-02.0-0-01-03.0-0-01.pdf output tabulate_obs_8li-coulex_m1-COMBO_210316.pdf
    """

    # content lists
    interaction_coulomb_list = [("Daejeon16",1),("JISP16",1),("LENPICchi2bi2C",0)]
    nuclide_observable_list = {

        # energies
        ((3,5),("energy",(2.0,0,1))) : (-50.,20.),
        ("diff", ((3,5),("energy",(1.0,0,1))), ((3,5),("energy",(2.0,0,1)))) : (0.,3.),
        ("diff", ((3,5),("energy",(1.0,0,2))), ((3,5),("energy",(1.0,0,1)))) : (0.,7.),
        
        # moments (M1)
        ((3,5),("moment","M1",(2.0,0,1))) : (0.,2.),
        ## ((3,5),("M1",(1.0,0,1))),
        
        # transitions (M1 downward)
        ((3,5),("rtp","M1",(2.0,0,1),(1.0,0,1))) : (-3.,0.),
        ##((3,5),("rtp","M1",(2.0,0,1),(3.0,0,1))),

        # ratios (M1)
        ("ratio", ((3,5),("rtp","M1",(2.0,0,1),(1.0,0,1))), ((3,5),("momentsqr","M1",(2.0,0,1)))) : (0.,4.),

        # moments (E2)
        ((3,5),("moment","E2p",(2.0,0,1))) : (0,5),
        ## ((3,5),("E2p",(1.0,0,1))),
        
        # transitions (E2 upward)
        ((3,5),("rtp","E2p",(1.0,0,1),(2.0,0,1))) : (0,5),
        ##((3,5),("rtp","E2p",(3.0,0,1),(2.0,0,1))),

        # ratios (E2)
        ("ratio", ((3,5),("rtp","E2p",(1.0,0,1),(2.0,0,1))), ((3,5),("momentsqr","E2p",(2.0,0,1)))) : (0,0.3),

        ("ratio", ((3,5),("moment","E2p",(1.0,0,1))), ((3,5),("moment","E2p",(2.0,0,1)))) : (0,0.6),
        
    }

    # plotting parameters
    Nmax_max=14
    hw_range_by_interaction_coulomb = {
        ("Daejeon16",1): (5,25),
        ("JISP16",1): (10,30),
        ("LENPICchi2bi2C",0): (15,35),
        }
    
    data_directory="data"

    # tabulate observables
    for interaction_coulomb in interaction_coulomb_list:
        for nuclide_observable in nuclide_observable_list:

            (interaction,coulomb) = interaction_coulomb
            
            # generate descriptor
            descriptor=mfdnres.data.hw_scan_descriptor(interaction_coulomb,nuclide_observable)
            print(descriptor)
            
            # select interaction/coulomb
            mesh_data_selected = mfdnres.analysis.selected_mesh_data(
                mesh_data,
                {"interaction":interaction,"coulomb":coulomb}
            )

            # tabulate and plot
            observable_data = mfdnres.data.make_hw_scan_data(
                mesh_data_selected,nuclide_observable
                )
            mfdnres.data.write_hw_scan_data(
                descriptor,observable_data,
                directory=data_directory
            )
            mfdnres.data.write_hw_scan_plot(
                descriptor,
                interaction_coulomb,nuclide_observable,
                observable_data,
                hw_range=hw_range_by_interaction_coulomb[interaction_coulomb],
                observable_range=nuclide_observable_list[nuclide_observable],
                Nmax_max=Nmax_max,
                directory=data_directory
            )


################################################################
# single plot -- basic example
################################################################

def make_basic_plot(mesh_data):
    """ Provides a bare-bones example of making a single analysis plot.
    """

    # TUTORIAL: Here is where the output will go.  Be sure to make this
    # subdirectory, if it doesn't already exist.
    
    data_directory="data/single"

    # select interaction/coulomb

    # TUTORIAL: The tabulation routines will take care of selecting the nuclide
    # and iterating over (Nmax,hw).  But right now our mesh contains results
    # from many interactions.  We have to filter for the interaction of
    # interest.
    
    interaction_coulomb = ("Daejeon16",1)
    (interaction,coulomb) = interaction_coulomb
    mesh_data_selected = mfdnres.analysis.selected_mesh_data(
        mesh_data,
        {"interaction":interaction,"coulomb":coulomb}
    )

    # define observable

    # TUTORIAL: It turns out to be most convenient to group our nuclide and
    # observble together in a tuple.  We usually need to know both together.
    #
    # Simple observables have the form:
    #
    #     (<observable>, [<operator>], qn, [qn])
    #
    # Examples:
    #
    #     ("energy", (1.5,1,1))  # energy of first 3/2- state
    #
    #     ("rtp", "E2p",  (1.5,1,1),  (2.5,1,1))  # E2 (proton) reduced transition probability B(E2;5/2->3/2)
    
    nuclide_observable = ((4,5), ("energy", (1.5,1,1)))
        
    # generate descriptor

    # TUTORIAL: The "descriptor" is a standard string representation of the
    # observable, which we will use in out output filenames.
    
    descriptor=mfdnres.data.hw_scan_descriptor(interaction_coulomb,nuclide_observable)
    print(descriptor)

    # tabulate

    # TUTORIAL: Here we digest the observable data into a pandas.DataFrame
    # object.  This is basically a "spreadsheet".  The observable values are in
    # a single column.  The "index" (row label), technically a pandas
    # "multi-index", is (Nmax,hw).  The column label is "values".

    observable_data = mfdnres.data.make_hw_scan_data(
        mesh_data_selected,nuclide_observable
        )

    # TUTORIAL: To see the pandas data frame...
    if False:
        print(observable_data)
        exit()

    # prune data

    # TUTORIAL: In general, we may want to prune to a given window in hw and
    # Nmax.  E.g., get rid of a stray high-Nmax test run, or extremes in hw.
    
    hw_range = (5.,25.)
    Nmax_max = 10
    observable_data = observable_data.loc[(slice(0,Nmax_max),slice(*hw_range)),:]
        
    # write data

    # TUTORIAL: Right here and now we write out the tabular data.  This way we
    # can go back and check numerical values, and we have it for posterity, to
    # submit as suppmemental data with a publication, re-plot it years from now
    # in other software, etc.
    
    mfdnres.data.write_hw_scan_data(
        descriptor,observable_data,
        directory=data_directory
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
        directory=data_directory
    )

################################################################
# mutiple analysis plots
################################################################

def make_analysis_plots(mesh_data):
    """ Make a series of standalone analysis plots.
    """

    # TUTORIAL: When you get the error message
    #
    #    FileNotFoundError: [Errno 2] No such file or directory: 'data/analysis/hw-scan_Daejeon16-1_Z04-N05-energy-01.5-1-01_table.dat'
    #
    # don't panic.  That's just because you forgot to make this directory...
    
    data_directory="data/analysis"

    # mesh definitions
    interaction_coulomb_list = [
        ("Daejeon16",1),
        ("JISP16",1),
        ##("LENPICchi2bi2C",0)
    ]
    hw_range_by_interaction_coulomb = {
        ("Daejeon16",1): (5.,25.),
        ("JISP16",1): (10.,30.),
        ("LENPICchi2bi2C",0): (15.,35.),
    }
    marker_by_interaction_coulomb = {
        ("Daejeon16",1): "o",
        ("JISP16",1): "s",
        ("LENPICchi2bi2C",0): "D",
    }
    Nmax_max = 10

    # observable definitions

    # TUTORIAL: This was originally just a list of (nuclide,observable) pairs.
    # But we've turned it into a dictionary, which we use to manually specify a
    # plot range, as well.
    #
    # Note the use of "compound" observables, automatically constructed as
    # differences or ratios of other (simple) observables.  The descriptor
    # string, plot labels, etc., will also be automatically constructed out of
    # those for the individual simple observables, as well.
    
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
    for interaction_coulomb in interaction_coulomb_list:
        for nuclide_observable in nuclide_observable_list:

            # select interaction/coulomb
            (interaction,coulomb) = interaction_coulomb
            mesh_data_selected = mfdnres.analysis.selected_mesh_data(
                mesh_data,
                {"interaction":interaction,"coulomb":coulomb}
            )

            # generate descriptor

            descriptor=mfdnres.data.hw_scan_descriptor(interaction_coulomb,nuclide_observable)
            print(descriptor)

            # tabulate
            observable_data = mfdnres.data.make_hw_scan_data(
                mesh_data_selected,nuclide_observable
            )
            
            # TUTORIAL: To see the pandas data frames...
            if False:
                print(observable_data)

            # prune data
            hw_range = hw_range_by_interaction_coulomb[interaction_coulomb]
            observable_data = observable_data.loc[(slice(0,Nmax_max),slice(*hw_range)),:]
        
            # write data
            mfdnres.data.write_hw_scan_data(
                descriptor,observable_data,
                directory=data_directory
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
                directory=data_directory
            )
    
            
################################################################
# main
################################################################

mesh_data=read_data()

##make_basic_plot(mesh_data)
make_analysis_plots(mesh_data)
