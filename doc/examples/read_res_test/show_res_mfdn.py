"""show_res_mfdn.py

    Provide simple inspection of parsed data from MFDn results file.

    Example usage:

        python3 show_res_mfdn.py data/mfdn/v15-h2/runmfdn13-mfdn15-Z3-N3-Daejeon16-coul1-hw15.000-a_cm40-Nmax02-Mj1.0-lan200-tol1.0e-06.res

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 01/11/23 (mac): Created, based on read_res_mfdn.py.

"""

import os
import sys

import mfdnres
import mfdnres.ncci

################################################################
# explore single mesh point
################################################################

def print_point(results_data):
    """Examine mfdn_results_data members and results of accessors.

    """

    print("params")
    ## print(results_data.params)
    for key in sorted(results_data.params.keys()):
        print("  {}: {}".format(key, results_data.params[key]))
    print()

    print("levels")
    print(results_data.levels)
    print()

    print("num_eigenvalues")
    print(results_data.num_eigenvalues)
    print()

    print("energies")
    print(results_data.energies)
    print()
    
    print("mfdn_level_properties -> T")
    print(results_data.mfdn_level_properties["T"])
    print()

    print("...")
    print()

    print("mfdn_tb_expectations -> r")
    print(results_data.mfdn_tb_expectations["r"])
    print()

    print("...")
    print()

    
################################################################
# main
################################################################

def main():
    filename = sys.argv[1]
    results_data, = mfdnres.input.read_file(filename, res_format="mfdn_v15", filename_format=None)
    
    print_point(results_data)

if (__name__ == "__main__"):
    main()
