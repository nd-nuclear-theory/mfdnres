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
    - 11/26/23 (mac): Add Nex decomposition.
    - 03/13/24 (mac): Illustrate selection of mesh point by parameters.
    - 07/03/23 (mac): Add occupations.
    - 10/02/26 (mac): Read single results file instead of slurping mesh.

"""

import os

import mfdnres
import mfdnres.ncci

################################################################
# reading data
################################################################

def read_results():
    """Read results from single results file.
    """

    data_dir = os.path.join("data", "mfdn", "v15-h2")
    filename = "runmfdn13-mfdn15-Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06.res" 

    print("Reading input file...")
    mesh_data = mfdnres.input.read_file(
        filename=os.path.join(data_dir, filename),
        ## res_format="mfdn_v15",
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
    """Examine mfdn_results_data members and results of accessors.

    """

    print("Inspecting mesh point")
    print()
    
    # parameters
    print("Params")
    print(mfdnres.analysis.dict_items(results_data.params))
    print()

    # examine data attributes
    print("Data attributes...")
    print("(These are the underlying data, which we examine here for illustration,")
    print("but normally you would access data using accessors instead, as shown later below.)")
    print()
    print("results_data.energies {}".format(results_data.energies))
    print("results_data.postprocessor_ob_rmes {}".format(results_data.postprocessor_ob_rmes))
    print("results_data.postprocessor_tb_rmes {}".format(results_data.postprocessor_tb_rmes))
    print("results_data.mfdn_level_occupations {}".format(results_data.mfdn_level_occupations))
    print("results_data.mfdn_tb_expectations {}".format(results_data.mfdn_tb_expectations))
    print()

    # access levels
    print("Test accessors (basic)...")
    print("results_data.levels {}".format(results_data.levels))
    print()
    
    # access ob moments
    print("Test accessors (one-body)...")
    print("M1 moment (native physical) {}".format(results_data.get_moment("M1-native", (1.0,0,1))))
    print("M1 moment (from dipole term rmes) {}".format(results_data.get_moment("M1", (1.0,0,1))))
    print("E2 moment (from dipole term rmes) {}".format(results_data.get_moment("E2p", (1.0,0,1))))
    print()

    # access tb expectations
    print("Test accessors (two-body)...")
    print("Rp {}".format(results_data.get_radius("rp", (1.0,0,1))))
    print()

    # access Nex decomposition
    print("Test Nex decomposition...")
    decomposition = results_data.get_decomposition("Nex", (1.0,0,1))
    print("Nex decomposition {}".format(decomposition))
    print()

    # access occupations
    print("Test occupations...")
    for species_code in ["p", "n"]:
        orbitals, occupations = results_data.get_occupations(species_code, (1.0,0,1))
        for orbital, occupation in zip(orbitals, occupations):
            n, l, j = orbital
            print("  {:1s}   {:2d} {:2d} {:4.1f}    {:.6f}".format(species_code, n, l, j, occupation))
    print()


################################################################
# main
################################################################

if (__name__ == "__main__"):

    results_data = read_results()
    explore_point(results_data)
