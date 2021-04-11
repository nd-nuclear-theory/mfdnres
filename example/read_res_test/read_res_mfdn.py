"""read_res_mfdn.py

    Provides simple example of reading and accessing MFDn run results.

    Required test data:
        type_specimens/mfdn/v15-h2/runmfdn13-mfdn15-Z3-N3-Daejeon16-coul1-hw15.000-a_cm40-Nmax02-Mj1.0-lan200-tol1.0e-06.res

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 09/17/20 (mac): Created.
    - 04/10/21 (mac): Update example file.

"""

import mfdnres
import mfdnres.ncci

################################################################
# reading data
################################################################

def read_data():
    """Read results.
    """

    print("Reading input file...")
    mesh_data = mfdnres.input.read_file(
        "type_specimens/mfdn/v15-h2/runmfdn13-mfdn15-Z3-N3-Daejeon16-coul1-hw15.000-a_cm40-Nmax02-Mj1.0-lan200-tol1.0e-06.res",
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

def explore_point(mesh_data):
    """Pick out single hw mesh point and examine spncci_results_data members and
    results of accessors.
    """

    # pick out mesh point manually
    results_data = mesh_data[0]

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
explore_point(mesh_data)
