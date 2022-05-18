"""read_res_mfdn.py

    Provides simple example of reading and accessing MFDn results.

    In practice, such results may need to be "merged" with results from the MFDn
    postprocessor.

    Required test data:

        data/mfdn/v15-h2/runmfdn13-mfdn15-Z3-N3-Daejeon16-coul1-hw15.000-a_cm40-Nmax02-Mj1.0-lan200-tol1.0e-06.res

        This example output is produced by mcscript-ncci/docs/examples/runmfdn13.py.

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 09/17/20 (mac): Created.
    - 05/10/21 (mac): Update example file.
    - 05/18/22 (mac): Update example file.

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
    data_dir = os.path.join("data","mfdn","v15-h2")
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
    
    return mesh_data

################################################################
# explore single mesh point
################################################################

def explore_point(results_data):
    """Examine mfdn_results_data members and results of accessors.

    """

    # examine data attributes
    print("Data attributes...")
    print("results_data.postprocessor_ob_rmes {}".format(results_data.postprocessor_ob_rmes))
    print("results_data.postprocessor_tb_rmes {}".format(results_data.postprocessor_tb_rmes))
    print()

    # access ob moments
    print("Test accessors (one-body)...")
    print("M1 moment (native physical) {}".format(results_data.get_moment("M1-native",(1.0,0,1))))
    print("M1 moment (from dipole term rmes) {}".format(results_data.get_moment("M1",(1.0,0,1))))
    print("E2 moment (from dipole term rmes) {}".format(results_data.get_moment("E2p",(1.0,0,1))))
    print()

    # access tb expectations
    print("Test accessors (two-body)...")
    print("Rp {}".format(results_data.get_radius("rp",(1.0,0,1))))
    print()
    
################################################################
# main
################################################################

# read data
mesh_data = read_data()
explore_point(mesh_data[0])
