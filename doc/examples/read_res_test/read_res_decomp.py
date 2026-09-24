"""read_res_decomp.py

    Provides simple example of reading and accessing MFDn decomposition run
    results.

    Required test data:
        runmfdndecomp02-decomp-Z3-N3-Daejeon16-coul1-hw15.000-Nmax02-J01.0-g0-n01-S-dNmax02-dlan0100.res

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 09/24/25 (mac): Created.

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
    data_dir = os.path.join("data", "mfdn-decomp")
    mesh_data = mfdnres.input.slurp_res_files(
        data_dir,
        res_format="decomp",
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
    print("results_data.mfdn_level_lanczos_decomposition_data {}".format(results_data.mfdn_level_lanczos_decomposition_data))
    print()

    # access decomposition data
    print("Test accessors (decomposition data)...")
    ##decomposition = results_data.get_lanczos_decomposition_data("S", (1.0,0,1))
    ##print("S decomposition {}".format(decomposition))
    print()

    
################################################################
# main
################################################################

# read data
mesh_data = read_data()
explore_point(mesh_data[0])
