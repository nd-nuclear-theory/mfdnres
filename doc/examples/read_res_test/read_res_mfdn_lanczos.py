"""read_res_mfdn_lanczos.py

    Provides simple example of reading MFDn Lanczos decomposition data files.

    Required test data:

        data/mfdn/v15-lanczos/*.lanczos

        This example output is produced by mcscript-ncci/docs/examples/runmfdn14.py.

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 10/12/23 (mac): Created.

"""

import os

import numpy as np

import mfdnres
import mfdnres.ncci
import mfdnres.decomposition

################################################################
# reading data
################################################################

def read_data():
    """Read results.
    """

    print("Reading lanczos files...")
    data_dir = os.path.join("data","mfdn","v15-lanczos")
    mesh_data = mfdnres.decomposition.slurp_lanczos_files(
        data_dir,
        filename_format="mfdn_format_7_ho",
        glob_pattern="*.lanczos",
        verbose=True
    )
    print()
    
    # diagnostic output -- FOR ILLUSTRATION ONLY
    print("Raw mesh (params)")
    for results_data in mesh_data:
        print(mfdnres.analysis.dict_items(results_data.params))
    print()

    # merge results objects giving decompositions from same mesh point
    mesh_data = mfdnres.analysis.merged_mesh(
        mesh_data,
        ("nuclide","interaction","coulomb","hw","Nmax"),
        postprocessor=mfdnres.ncci.augment_params_with_parity,
        verbose=False
    )

    # diagnostic output -- FOR ILLUSTRATION ONLY
    print("Merged mesh (keys)")
    mfdnres.analysis.mesh_key_listing(
        mesh_data,
        ("nuclide","interaction","coulomb","hw","Nmax","parity"),
        verbose=True
    )
    print()
    
    return mesh_data

################################################################
# explore single mesh point
################################################################

def explore_point(results_data):
    """Examine mfdn_results_data members and results of accessors.

    """

    # set numpy formatting to abbreviate array output
    np.set_printoptions(threshold=0, edgeitems=2)
    
    # examine data attributes
    print("Data attributes...")
    print("results_data.mfdn_level_lanczos_decomposition_data {}".format(results_data.mfdn_level_lanczos_decomposition_data))
    print()

    # access data
    print("Test accessor...")
    print("Lanczos decomposition filename {}".format(results_data.get_lanczos_decomposition_filename("Nex", (1.0,0,1))))
    alpha, beta = results_data.get_lanczos_decomposition_alpha_beta("Nex", (1.0,0,1))
    print("Lanczos decomposition alpha & beta\n    {}\n    {}".format(alpha, beta))
    print()

    # do decomposition
    alpha_beta = results_data.get_lanczos_decomposition_alpha_beta("Nex", (1.0,0,1))
    eigenvalue_label_dict = mfdnres.decomposition.eigenvalue_label_dict_Nex(Nmax=2)
    decomposition = mfdnres.decomposition.generate_decomposition(alpha_beta, eigenvalue_label_dict)
    print("Nex decomposition by Lanczos")
    mfdnres.decomposition.print_decomposition(decomposition)

    
################################################################
# main
################################################################

if (__name__ == "__main__"):

    # read data
    mesh_data = read_data()
    explore_point(mesh_data[0])
