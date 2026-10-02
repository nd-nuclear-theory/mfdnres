"""read_res_transitions.py

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
    - 10/02/26 (mac): Read single results file instead of slurping mesh.

"""

import os

import mfdnres
import mfdnres.ncci


################################################################
# reading data
################################################################

def read_results():
    """Read results from one-body and two-body results files and merge.
    """
    
    print("Reading input file...")
    data_dir = os.path.join("data", "mfdn-transitions")
    mesh_data = mfdnres.input.slurp_res_files(
        data_dir,
        ## res_format="mfdn_v15",
        ## filename_format="mfdn_format_7_ho",
        filename_format="ALL",
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
    print("Example of failure for invalid pair of initial and final states...")
    print("E2 rme {}".format(results_data.get_rme("E2p",((1.0,0,1),(1.0,1,1)),verbose=True)))  # invalid state pair
    print()

    
################################################################
# main
################################################################

if (__name__ == "__main__"):

    results_data = read_results()
    explore_point(results_data)
