"""read_res_transitions_spamp.py

    Provides simple example of reading and accessing spectroscopic amplitude
    results in a putative format for future spectroscopic amplitude calculations
    by mfdn-transitions.  In initial implmentation, such results are obtained by
    digesting results from rhodium.

    In practice, such results may need to be "merged" with results from mfdn.

    Required test data:
        data/mfdn-transitions-spamp/runmac0908-transitions-spamp-Z3-N3-Daejeon16-coul1-hw15.000-Nmax02-Mj0.5.res
        data/mfdn-transitions-spamp/runmac0908-transitions-spamp-Z3-N3-Daejeon16-coul1-hw15.000-Nmax02-Mj1.5.res

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 07/03/25 (mac): Created.

"""

import os

import mfdnres
import mfdnres.ncci

################################################################
# reading data
################################################################

def read_data():
    """Read results.
    """

    print("Reading input file...")
    data_dir = os.path.join("data","mfdn-transitions-spamp")
    mesh_data = mfdnres.input.slurp_res_files(
        data_dir,
        res_format="mfdn_v15",
        filename_format="mfdn_format_7_ho",
        verbose=True
    )
    print()
    
    # diagnostic output -- FOR ILLUSTRATION ONLY
    print("Raw mesh (params)")
    for results_data in mesh_data:
        print(mfdnres.analysis.dict_items(results_data.params))
    print()
    
    # merge results data
    print("Merging mesh points...")
    mesh_data = mfdnres.analysis.merged_mesh(
        mesh_data,
        ("nuclide","interaction","coulomb","hw","Nmax"),
        postprocessor=mfdnres.ncci.augment_params_with_parity,
        verbose=False
    )
    print()

    # diagnostic output -- FOR ILLUSTRATION ONLY
    print("Merged mesh (params)")
    for results_data in mesh_data:
        print(mfdnres.analysis.dict_items(results_data.params))
    print()

    return mesh_data

################################################################
# explore single mesh point
################################################################

def explore_point(results_data):
    """Examine mfdn_results_data members and results of accessors, for MFDn postprocessor results.
    """

    # pick out mesh point manually
    results_data = mesh_data[0]

    # examine data attributes
    print("Data attributes...")
    print("results_data.postprocessor_spectroscopic_amplitudes {}".format(results_data.postprocessor_spectroscopic_amplitudes))
    print()

    # access ob rmes
    print("Test accessors (spectroscopic amplitudes)...")
    amplitudes = results_data.get_spectroscopic_amplitudes((0,+1), ((0.5, 1, 1), (1.0, 0, 1)))
    for orbital, value in amplitudes.items():
        n, l, j = orbital
        print("  {:1d} {:1d} {:3.1f}    {:+e}".format(n, l, j, value))
    print()

    
################################################################
# main
################################################################

# read data
mesh_data = read_data()
explore_point(mesh_data[0])
