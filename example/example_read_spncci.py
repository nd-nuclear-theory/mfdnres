"""example_read_spncci.py

    Provides simple example of reading and accessing SpNCCI run results.

    Required test data:
        type_specimens/spncci/runmac0420-Z3-N3-JISP16-0-Nsigmamax02-Nmax02.res

    Mark A. Caprio
    University of Notre Dame

    Language: Python 3

    - 08/14/20 (mac): Created.

"""

import mfdnres
import mfdnres.decomposition

################################################################
# reading data
################################################################

def read_data():
    """Read spncci results.
    """

    print("Reading input file...")
    mesh_data = mfdnres.input.read_file(
        "type_specimens/spncci/runmac0420-Z3-N3-JISP16-0-Nsigmamax02-Nmax02.res",
        "spncci",
        filename_format="spncci",
        verbose=True
    )
    print()
    
    # diagnostic output -- FOR ILLUSTRATION ONLY
    print("Examining mesh_data")
    print()
    print("mesh_data {}".format(mesh_data))
    print()
    print("Mesh points (params)")
    for results_data in mesh_data:
        ## print(results_data.params)
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
    ## results_data = mesh_data[4]

    # pick out mesh point by key
    print("Selecting single results data object...")
    results_dict = mfdnres.analysis.make_results_dict(
        mesh_data,
        (("Nsigmamax",int),("Nmax",int),("hw",float))
        )
    ##print("results_dict {}".format(results_dict))
    results_data = results_dict[(2,2,20.)]
    print("results_data {}".format(results_data))
    print("results_data.params {}".format(results_data.params))
    print()

    # examine inherited data attributes
    print("Some inherited data attributes (from ResultsData)...")
    print("results_data.filename {}".format(results_data.filename))
    print("results_data.energies {}".format(results_data.energies))
    print("results_data.levels {}".format(results_data.levels))
    print()

    # examine spncci-specific data attributes
    print("Some SpNCCI-specific data attributes (from SpNCCIResultsData)...")
    print("results_data.Jgex_values {}".format(results_data.Jgex_values))
    print("results_data.baby_spncci_listing (showing first few entries)\n{}".format(results_data.baby_spncci_listing[:4]))
    print()

    # access decomposition data
    print("Decompositions -- looking at raw data (deprecated but informative!) then using accessor for proper access")
    print("results_data.decompositions -- but use accessor instead!")
    print("results_data.decompositions (type) {}".format(type(results_data.decompositions)))
    print("results_data.decompositions (keys) {}".format(results_data.decompositions.keys()))
    print("results_data.decompositions (keys for Nex) {}".format(results_data.decompositions["Nex"].keys()))
    print("results_data.decompositions (data for Nex for (1.0,0) space)\n{}".format(results_data.decompositions["Nex"][(1.0,0)]))
    print("results_data.get_decomposition (results for Nex for state (1.0,0,1))\n{}".format(results_data.get_decomposition("Nex",(1.0,0,1))))
    print("results_data.get_decomposition (results for BabySpNCCI for state (1.0,0,1))\n{}".format(results_data.get_decomposition("BabySpNCCI",(1.0,0,1))))

    # process decomposition into dictionary
    labeled_decomposition = mfdnres.decomposition.labeled_decomposition(
        results_data.get_decomposition_labels("BabySpNCCI"),
        results_data.get_decomposition("BabySpNCCI",(1.0,0,1))
        )
    print("BabySpNCCI decomposition dictionary (for state (1.0,0,1))")
    mfdnres.decomposition.print_decomposition(labeled_decomposition)
    
################################################################
# main
################################################################

# read data
mesh_data = read_data()
explore_point(mesh_data)
