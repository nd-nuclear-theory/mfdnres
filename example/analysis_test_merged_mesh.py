"""analysis_test_merged_mesh.py

    Testbed for mesh merger.

    Required test data:
        results/mcaprio/mfdn/runmac0506

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 07/11/19 (mac): Created (tabulate_band_test), based on code in
      tabulate_band_be_bebands.py.

"""

import mfdnres
import mfdnres.ncci

################################################################
# global data file path configuration
################################################################

def read_data():
    """Read MFDn results for high-Lanczos runs, and merge mfdn and observable
    results, as well as M runs.
    """

    full_data = []

    run_list = [

        "mac0506",  # 10Be+
    ]

    data_dir_list = [
        mfdnres.res.res_file_directory("mcaprio","mfdn",run)
        for run in run_list
    ]
    res_format = "mfdn_v15"
    filename_format="mfdn_format_7_ho"
    mesh_data = mfdnres.res.slurp_res_files(
        data_dir_list,res_format,filename_format,
        glob_pattern="*-*-*.res",
        verbose=True
    )

    print("Raw mesh (params)")
    for results_data in mesh_data:
        print(mfdnres.analysis.dict_items(results_data.params))

    print("Raw mesh (keys)")
    mfdnres.analysis.mesh_key_listing(mesh_data,("nuclide","interaction","coulomb","hw","Nmax","parity","M","code_name"),verbose=True)

    # merge results data (from different M and mfdn/obscalc-ob)
    mesh_data = mfdnres.analysis.merged_mesh(
        mesh_data,
        ("nuclide","interaction","coulomb","hw","Nmax"),
        postprocessor=mfdnres.ncci.augment_params_with_parity,
        verbose=False
    )

    print("Merged mesh (keys)")
    mfdnres.analysis.mesh_key_listing(
        mesh_data,
        ("nuclide","interaction","coulomb","hw","Nmax","parity"),
        verbose=True
    )

    # select case of interest
    selector = {"interaction":"Daejeon16","hw":15.0}
    mesh_data = mfdnres.analysis.selected_mesh_data(full_data,selector,verbose=False)


    return mesh_data

################################################################
# main
################################################################

# read data
mesh_data = read_data()
