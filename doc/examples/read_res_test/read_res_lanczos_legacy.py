"""read_res_lanczos_legacy.py

    Provides simple example of reading MFDn Lanczos decomposition data files.

    This is the legacy storage scheme, supplanted by reading of a decomp res file.

    Required test data:

        data/mfdn/v15-lanczos/*.lanczos

        This example output is produced by mcscript-ncci/docs/examples/runmfdn14.py.

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 10/12/23 (mac): Created.
    - 10/02/26 (mac): Read single Lanczos file instead of slurping mesh.
        Rename from read_res_mfdn_lanczos.py to read_res_lanczos_legacy.py. 
"""

import os

import numpy as np

import mfdnres
import mfdnres.ncci
import mfdnres.decomposition

################################################################
# reading data
################################################################

def read_results():
    """Read results from single Lanczos file.
    """

    print("Reading lanczos file...")
    data_dir = os.path.join("data","mfdn","v15-lanczos")
    filename = "runmfdndecomp02-decomp-Z3-N3-Daejeon16-coul1-hw15.000-Nmax02-J01.0-g0-n01-S-dNmax02-dlan0100.lanczos" 
    
    mesh_data = mfdnres.decomposition.slurp_lanczos_files(
        data_dir,
        filename_format="mfdn_format_7_ho",
        glob_pattern=filename,
        verbose=True
    )

    results_data = mesh_data[0]
    return results_data

################################################################
# explore single mesh point
################################################################

def explore_point(results_data):
    """Examine mfdn_results_data members and results of accessors.

    """

    # examine data attributes
    print("Data attributes...")
    print("results_data.mfdn_level_lanczos_decomposition_data {}".format(results_data.mfdn_level_lanczos_decomposition_data))
    print()

    # access data
    print("Test accessor...")
    alpha_beta = results_data.get_lanczos_decomposition_alpha_beta("S", (1.0,0,1))
    print("Lanczos decomposition alpha & beta\n{}".format(alpha_beta))
    print()

    # do decomposition -- DEPRECATED syntax
    ## alpha_beta = results_data.get_lanczos_decomposition_alpha_beta("Nex", (1.0,0,1))
    ## eigenvalue_label_dict = mfdnres.decomposition.eigenvalue_label_dict_Nex(Nmax=2)
    ## decomposition = mfdnres.decomposition.generate_decomposition(alpha_beta, eigenvalue_label_dict)
    ## print("Nex decomposition by Lanczos")
    ## mfdnres.decomposition.print_decomposition(decomposition)

    
################################################################
# main
################################################################

if (__name__ == "__main__"):

    results_data = read_results()
    with np.printoptions(edgeitems=4, threshold=10):
        explore_point(results_data)
