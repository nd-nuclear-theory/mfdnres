""" mfdnres_test.py -- basic tests of mfdnres package

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    6/2/15 (mac): Initiated.
    Last modified 6/6/15.
    
"""

import os

import mfdnres

def test_basic_import_v14b05():
    """ Test single-file import for MFDn v14b06.
    """

    print("test_basic_import_v14b05")


    import os

    filename = os.path.join("example","type_specimens","mfdn","v14b05-vxx","MFDn.res.Z4.N5.JISP16.Nmin1.Nm13.hw20.0.La500.St06.tol1e-6") # VXX format, no trans
    res_format = "mfdn_v14b05"
    res_filename_format = None

    results_list = mfdnres.res.read_file(filename,res_format,res_filename_format,verbose=False)
    results = results_list[0]

    results = mfdnres.res.read_file(filename,res_format,verbose=False)[0]
    print(results.get_levels())
    print(results.states[(1/2,1,1)].properties)
    print(results.states[(1/2,1,1)].obo)
    ## print(results.transitions)
    print("states {}, moments {}, transitions {}".format(len(results.states),len(results.moments),len(results.transitions)))

def test_basic_import_v14b06():
    """ Test single-file import for MFDn v14b06.
    """

    print("test_basic_import_v14b06")

    ## data_dir = mfdnres.res.res_file_directory("mcaprio","mfdn","0355",run_results_are_in_subdir=False)
    ## filename = os.path.join(data_dir,"run0355-mfdn-Z2-N2-JISP16-1-hw20.000-aL20-Nmax04-MM0-lan500.res") # no M1 moments, no trans

    ## data_dir = mfdnres.res.res_file_directory("mcaprio","mfdn","0363",run_results_are_in_subdir=False)
    ## filename = os.path.join(data_dir,"run0363-mfdn-Z4-N3-N2LOopt500-0-hw20.000-aL100-Nmax10-MM1-lan1000.res") # dipole moments but not trans

    data_dir = mfdnres.res.res_file_directory("mcaprio","mfdn","0352",run_results_are_in_subdir=False)
    filename = os.path.join(data_dir,"run0352-mfdn-Z4-N5-JISP16-1-hw20.000-aL100-Nmax10-MM1-lan1000.res") # no M1 moment

    res_format = "mfdn_v14b06"

    ## filename = os.path.join("type_specimens","v14b06h2_run0381-mfdn-Z4-N6-N2LOopt500-0-hw20.000-aL100-Nmax07-Mj4.0-lan1500.res") # H2 format
    ## print(data.states[(1/2,1,1)].properties)
    ## print(data.states[(1/2,1,1)].obo)


    results = mfdnres.res.read_file(filename,res_format,verbose=False)[0]
    print(results.get_levels())
    print(results.states[(1/2,0,1)].properties)
    print(results.states[(1/2,0,1)].obo)
    ## print(results.transitions)
    print("states {}, moments {}, transitions {}".format(len(results.states),len(results.moments),len(results.transitions)))

def test_basic_import_v15():
    """ Test single-file import for MFDn v15.

pfasano/mfdn/runpjf0015
runpjf0015-mfdn15-Z3-N4-JISP16-coul1-hw20.000-a_cm40-Nmax02-Mj0.5-lan1000-tol1.0e-06-natorb-no0.res

WIP
    """

    print("test_basic_import")

    ## data_dir = mfdnres.res.res_file_directory("mcaprio","mfdn","0355",run_results_are_in_subdir=False)
    ## filename = os.path.join(data_dir,"run0355-mfdn-Z2-N2-JISP16-1-hw20.000-aL20-Nmax04-MM0-lan500.res") # no M1 moments, no trans

    ## data_dir = mfdnres.res.res_file_directory("mcaprio","mfdn","0363",run_results_are_in_subdir=False)
    ## filename = os.path.join(data_dir,"run0363-mfdn-Z4-N3-N2LOopt500-0-hw20.000-aL100-Nmax10-MM1-lan1000.res") # dipole moments but not trans

    data_dir = mfdnres.res.res_file_directory("mcaprio","mfdn","0352",run_results_are_in_subdir=False)
    filename = os.path.join(data_dir,"run0352-mfdn-Z4-N5-JISP16-1-hw20.000-aL100-Nmax10-MM1-lan1000.res") # no M1 moment

    res_format = "mfdn_v14b06"

    results = mfdnres.res.read_file(filename,res_format,verbose=False)[0]
    print(results.get_levels())
    print(results.states[(1/2,0,1)].properties)
    print(results.states[(1/2,0,1)].obo)
    ## print(results.transitions)
    print("states {}, moments {}, transitions {}".format(len(results.states),len(results.moments),len(results.transitions)))

def test_basic_import_spncci():
    """ Test single-file import for spncci.
    """

    print("test_basic_import_spncci")

    filename = os.path.join("example","type_specimens","spncci","runmac0420-Z3-N3-Nsigmamax02-Nmax02.res")
    res_format = "spncci"
    res_filename_format = None

    results_list = mfdnres.res.read_file(filename,res_format,res_filename_format,verbose=False)
    results = results_list[0]

    print(results.get_levels())
    ##print("states {}, moments {}, transitions {}".format(len(results.states),len(results.moments),len(results.transitions)))


def test_directory_slurp():
    """ Test slurping input directory and level output functions.
    """

    print("test_directory_slurp")

    data_dir = mfdnres.res.res_file_directory("mcaprio","mfdn","0352",run_results_are_in_subdir=False)
    res_format = "mfdn_v14b06"

    # slurp directory
    mesh_data = mfdnres.res.slurp_res_files(data_dir,res_format)
    KEY_DESCRIPTOR_INTERACTION_NMAX_HW =  (("interaction",str),("Nmax",int),("hw",float))
    print(mesh_data[0].params)
    results_dict = mfdnres.analysis.make_results_dict(mesh_data,KEY_DESCRIPTOR_INTERACTION_NMAX_HW)

    #("interaction","Nmax","hw")
    # dump levels
    print("mesh_data:",mesh_data)
    print("results_dict:",results_dict)
    results = results_dict[("JISP16",10,20.000)]
    mfdnres.analysis.write_level_table(
        results,
        "test-levels.dat"
    )

def test_legacy_band_output():
    """
    DEPRECATED: functions to reproduce berotor style band moment and transition files
    """

    # slurp directory
    data = {}
    mfdnres.analysis.import_res_files(
        data,
        filename = r"c:\work\research\data\mfdn\run0352",
        key_fields = ("interaction","Nmax","hw"),
        descriptor_format="format_5_ho",
        res_format="mfdn_v14b06",
        verbose=False
    )

    # write output
    members_9be_band1 = [(1.5,0,1),(2.5,0,1),(3.5,0,1),(4.5,0,1)]
    mfdnres.analysis.write_level_table(
        "mfdnres_tests_level_table_band1.dat",
        results,
        levels=members_9be_band1
    )
    mfdnres.analysis.write_band_trans_table(
        "mfdnres_tests_trans_band1-m1-dj1.dat",
        results,
        trans="M1",dJ=1,Mj=1/2,
        levels=members_9be_band1,
        signs = None
    )
    mfdnres.analysis.write_band_trans_table(
        "mfdnres_tests_trans_band1-e2-dj1.dat",
        results,
        trans="E2",dJ=1,Mj=1/2,
        levels=members_9be_band1,
        signs = None
    )
    mfdnres.analysis.write_band_trans_table(
        "mfdnres_tests_trans_band1-e2-dj2.dat",
        results,
        trans="E2",dJ=2,Mj=1/2,
        levels=members_9be_band1,
        signs = None
    )

def test_band_file():
    """
    test of input of band configuration file
    """

    # reading band configuration file
    filename = r"band-Z4-N3-JISP16-1-hw20.000-Nmax10-band1.cfg"
    mfdnres.analysis.BandDefinition(filename)

def test_band_output():
    """
    validation of new consolidated band output
    on old JISP+NC run0339
    """

    # slurp directory
    data = {}
    mfdnres.analysis.import_res_files(
        data,
        filename = r"c:\work\research\data\mfdn\run0339",
        key_fields = ("interaction","Nmax","hw"),
        descriptor_format="format_5_ho",
        res_format="mfdn_v14b06",
        verbose=False
    )

    # reading band configuration file
    band = mfdnres.analysis.BandDefinition("band-Z4-N3-JISP16-0-hw22.500-Nmax10-band1.cfg")

    # writing band transition output
    results = data[("JISP16","10","22.500")]
    mfdnres.analysis.write_band_table(results,"test-band.dat",band)

    # writing band fit output
    results = data[("JISP16","10","22.500")]
    mfdnres.analysis.write_band_fit_parameters(results,"test-band-fit.dat",band)

def test_am_output():
    """
    test of angular momentum output
    run0359 for 3He Nmax04
    """

    # slurp directory
    data = {}
    mfdnres.analysis.import_res_files(
        data,
        filename = r"c:\work\research\data\mfdn\run0359",
        key_fields = ("interaction","Nmax","hw"),
        descriptor_format="format_5_ho",
        res_format="mfdn_v14b06",
        verbose=False
    )

    # dump levels
    results = data[("JISP16","04","20.000")]
    mfdnres.analysis.write_level_am_table(
        results,
        "test-levels-am.dat"
    )


if (__name__ == "__main__"):

    # basic import tests
    test_basic_import_v14b05()
    test_basic_import_v14b06()
    test_basic_import_spncci()

    # directory slurp test
    ##test_directory_slurp()

    # analysis tests
    ##test_band_file()
    ##test_band_output()
    ##test_am_output()

