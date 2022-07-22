"""read_res_mfdn_transitions.py

    Provides simple example of reading and accessing MFDn postprocessor results.

    In practice, such results may need to be "merged" with results from mfdn.

    Required test data:
        data/mfdn-transitions/runtransitions00-transitions-ob-Z3-N3-Daejeon16-coul1-hw15.000-Nmax02.res
        data/mfdn-transitions/runtransitions00-transitions-tb-Z3-N3-Daejeon16-coul1-hw15.000-Nmax02.res

        This example output is produced by mcscript-ncci/docs/examples/runtransitions00.py.

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 09/17/20 (mac): Created.
    - 05/18/22 (mac): Update example file.
    - 07/12/22 (mac): Provide example use of two-body RME accessor.

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
    data_dir = os.path.join("data","mfdn-transitions")
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
    print("results_data.postprocessor_ob_rmes {}".format(results_data.postprocessor_ob_rmes))
    print("results_data.postprocessor_tb_rmes {}".format(results_data.postprocessor_tb_rmes))
    print()

    # access ob rmes
    print("Test accessors (one-body)...")
    print("M1 moment (from dipole term rmes) {}".format(results_data.get_moment("M1",(1.0,0,1))))
    print("M1 rme (from dipole term rmes) {}".format(results_data.get_rme("M1",((1.0,0,1),(1.0,0,1)))))
    print("E2 moment {}".format(results_data.get_moment("E2p",(1.0,0,1))))
    print("E2 rme {}".format(results_data.get_rme("E2p",((1.0,0,1),(1.0,0,1)))))
    print()

    # access tb rmes
    print("Test accessors (two-body)...")
    print("QxQ_0 rme {}".format(results_data.get_rme("QxQ_0",((2.0,0,1),(2.0,0,1)),rank="tb")))
    print()

    print("Test get_rme verbose mode...")
    print("E2 rme {}".format(results_data.get_rme("E2p",((1.0,0,1),(1.0,0,1)),verbose=True)))
    print("E2 rme {}".format(results_data.get_rme("E2p",((1.0,0,1),(1.0,1,1)),verbose=True)))  # invalid state pair
    print()
    
################################################################
# main
################################################################

# read data
mesh_data = read_data()
explore_point(mesh_data[0])
