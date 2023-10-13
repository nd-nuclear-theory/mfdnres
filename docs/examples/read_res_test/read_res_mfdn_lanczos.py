"""read_res_mfdn_lanczos.py

    Provides simple example of reading MFDn Lanczos decomposition data filenames.

    Required test data:

        data/mfdn/v15-lanczos/*.lanczos

        This example output is produced by mcscript-ncci/docs/examples/runmfdn14.py.

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 10/12/23 (mac): Created.

"""

import os

import mfdnres
import mfdnres.ncci
import mfdnres.decomposition

################################################################
# reading data
################################################################

def read_data():
    """Read results.
    """

    print("Reading lanczos filenames...")
    data_dir = os.path.join("data","mfdn","v15-lanczos")
    mesh_data = mfdnres.decomposition.slurp_lanczos_filenames(
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

    # examine data attributes
    print("Data attributes...")
    print("results_data.mfdn_level_lanczos_decomposition_filenames {}".format(results_data.mfdn_level_lanczos_decomposition_filenames))
    print()

    # access filename
    print("Test accessor...")
    print("Lanczos decomposition filename {}".format(results_data.get_lanczos_decomposition_filename("Nex",(1.0,0,1))))
    print()
    
################################################################
# main
################################################################

if (__name__ == "__main__"):

    # read data
    mesh_data = read_data()
    explore_point(mesh_data[0])
