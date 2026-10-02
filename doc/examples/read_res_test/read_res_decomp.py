"""read_res_decomp.py

    Provides simple example of reading and accessing MFDn decomposition run
    results.

    Required test data:
        runmfdndecomp02-decomp-Z3-N3-Daejeon16-coul1-hw15.000-Nmax02-J01.0-g0-n01-S-dNmax02-dlan0100.res

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 09/24/26 (mac): Created.
    - 10/02/26 (mac): Read single results file instead of slurping mesh.

"""

import os

import numpy as np

import mfdnres
import mfdnres.decomposition
import mfdnres.ncci


################################################################
# reading data
################################################################

def read_results():
    """Read results from single results file.
    """

    data_dir = os.path.join("data", "mfdn-decomp")
    filename = "runmfdndecomp02-decomp-Z3-N3-Daejeon16-coul1-hw15.000-Nmax02-J01.0-g0-n01-S-dNmax02-dlan0100.res" 

    print("Reading input file...")
    mesh_data = mfdnres.input.read_file(
        filename=os.path.join(data_dir, filename),
        ## res_format="decomp",
        ## filename_format="mfdn_format_7_ho",
        filename_format="ALL",
        verbose=True,
    )
    results_data = mesh_data[0]
    
    return results_data


################################################################
# explore single mesh point
################################################################

def explore_point(results_data):
    """Examine mfdn_results_data members and results of accessors, for MFDn postprocessor results.
    """

    # examine data attributes
    print("Data attributes...")
    print("results_data.mfdn_level_lanczos_decomposition_data {}".format(results_data.mfdn_level_lanczos_decomposition_data))
    print()

    # access decomposition data
    decomposition_type = "S"
    qn = (1.0,0,1)
    print("Test accessors (decomposition data)...")
    decomposition_data = results_data.get_lanczos_decomposition_data(decomposition_type, qn)
    print("Decomposition data object\n{}".format(decomposition_data))
    alpha_beta = results_data.get_lanczos_decomposition_alpha_beta(decomposition_type, qn)
    print("Lanczos alpha-beta coefficients\n{}".format(alpha_beta))
    print()
    num_eigenvalues = decomposition_data["statistics"]["num_eigenvalues"]
    raw_decomposition = mfdnres.decomposition.generate_raw_decomposition(alpha_beta, lanczos_iterations=num_eigenvalues)
    print("Raw decomposition\n{}".format(raw_decomposition))
    
    
################################################################
# main
################################################################

if (__name__ == "__main__"):

    results_data = read_results()
    with np.printoptions(edgeitems=4, threshold=10):
        explore_point(results_data)
