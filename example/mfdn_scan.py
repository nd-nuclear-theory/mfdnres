""" mfdn_scan.py

    Example convergence analysis of MFDn runs.

    Mark A. Caprio
    University of Notre Dame

    10/12/17 (mac): Created, based on analysis/spncci:convergence_analysis.py.
    10/23/17 (mac): Add radius scan example.
    04/02/19 (mac): Add angular momentum scan example.

"""

import mfdnres

# standard sorting and tabulation key
KEY_DESCRIPTOR_NMAX_HW = (("Nmax",int),("hw",float))
KEY_DESCRIPTOR_NORB_NMAX_HW = (("natural_orbital_iteration",int),("Nmax",int),("hw",float))

################################################################
# output functions to generate table and write to file
################################################################

def make_energy_table_file(mesh_data,nuclide,interaction_coulomb,qn):
    """ Tabulate energies for hw scan.

    Tabulation format:
        Nmax hw value

    Arguments:
        mesh_data (list of ResultsData): data set to include
        nuclide (tuple): (Z,N)
        interaction_coulomb (tuple): (interaction,use_coulomb)
        qn (tuple): (J,g,n) state label

    """

    (interaction,coulomb) = interaction_coulomb
    (J,g,n) = qn
    output_file_name="data/data-Z{nuclide[0]:02d}-N{nuclide[1]:02d}-{interaction:s}-{coulomb:1d}-scan-e-{J:04.1f}-{g:1d}-{n:02d}.dat".format(
        nuclide=nuclide,
        interaction=interaction,coulomb=coulomb,
        J=J,g=g,n=n
    )

    table = mfdnres.analysis.make_energy_table(mesh_data,KEY_DESCRIPTOR_NMAX_HW,qn)
    mfdnres.tools.write_table(
        output_file_name,"{:2d} {:7.3f} {:7.3f}",
        table
    )

def make_radius_table_file(mesh_data,nuclide,interaction_coulomb,radius_type,qn):
    """ Tabulate radius for hw scan.

    Tabulation format:
        Nmax hw value

    Arguments:
        mesh_data (list of ResultsData): data set to include
        nuclide (tuple): (Z,N)
        interaction_coulomb (tuple): (interaction,use_coulomb)
        radius_type (str): radius type code (rp/rn/r)
        qn (tuple): (J,g,n) state label
    
    """

    (interaction,coulomb) = interaction_coulomb
    (J,g,n) = qn
    output_file_name="data/data-Z{nuclide[0]:02d}-N{nuclide[1]:02d}-{interaction:s}-{coulomb:1d}-scan-{radius_type}-{J:04.1f}-{g:1d}-{n:02d}.dat".format(
        nuclide=nuclide,
        interaction=interaction,coulomb=coulomb,
        radius_type=radius_type,
        J=J,g=g,n=n
    )

    table = mfdnres.analysis.make_radius_table(mesh_data,KEY_DESCRIPTOR_NMAX_HW,"r",qn)
    mfdnres.tools.write_table(
        output_file_name,"{:2d} {:7.3f} {:7.3f}",
        table
    )

def make_am_table_file(mesh_data,nuclide,interaction_coulomb,qn):
    """ Tabulate effective angular momenta for hw scan.

    Tabulation format:
        Nmax hw value

    Arguments:
        mesh_data (list of ResultsData): data set to include
        nuclide (tuple): (Z,N)
        interaction_coulomb (tuple): (interaction,use_coulomb)
        qn (tuple): (J,g,n) state label
    
    """

    (interaction,coulomb) = interaction_coulomb
    (J,g,n) = qn
    output_file_name="data/data-Z{nuclide[0]:02d}-N{nuclide[1]:02d}-{interaction:s}-{coulomb:1d}-scan-am-{J:04.1f}-{g:1d}-{n:02d}.dat".format(
        nuclide=nuclide,
        interaction=interaction,coulomb=coulomb,
        J=J,g=g,n=n
    )

    table = mfdnres.analysis.make_am_table(mesh_data,KEY_DESCRIPTOR_NMAX_HW,qn)
    mfdnres.tools.write_table(
        output_file_name,"{:2d} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f}",
        table
    )

################################################################
# example control code
################################################################

def make_tabulations_pjf0007():
    """ Examples from mfdn_v14 run.  Tabulations are for 3He J=1/2 ground state.
    """

    # import data
    data_dir = mfdnres.res.res_file_directory("pfasano","mfdn","pjf0007")
    res_format = "mfdn_v14b06"
    filename_format="mfdn_format_7_ho"
    print("directory {}".format(data_dir))
    # Note: Only read "*-no0.res" files (oscillator runs).  The
    # "*-no1.res" files will crash the parser due to missing "One-body
    # radii & E2 moments".
    mesh_data = mfdnres.res.slurp_res_files(data_dir,res_format,filename_format,glob_pattern="*-no0.res",verbose=False)

    print("directory {}".format(data_dir))
    print("full mesh {}".format(len(mesh_data)))
    ## print([mesh_point.params for mesh_point in mesh_data])

    # select and sort
    nuclide = (2,1)
    interaction_coulomb = ("JISP16",1)
    selector = {"nuclide":nuclide,"interaction":interaction_coulomb[0],"coulomb":interaction_coulomb[1],"natural_orbital_iteration":0}
    mesh_data = mfdnres.analysis.selected_mesh_data(mesh_data,selector,verbose=False)
    print("selected mesh {}".format(len(mesh_data)))
    mesh_data = mfdnres.analysis.sorted_mesh_data(mesh_data,KEY_DESCRIPTOR_NMAX_HW)
    ##print([mesh_point.params for mesh_point in mesh_data])
    
    # make energy table
    qn = (0.5,0,1)
    make_energy_table_file(mesh_data,nuclide,interaction_coulomb,qn)

    # make radius table
    qn = (0.5,0,1)
    make_radius_table_file(mesh_data,nuclide,interaction_coulomb,"rp",qn)


def make_tabulations_pjf0015():
    """ Examples from mfdn_v15 run.  Tabulations are for 7Li J=3/2 ground state.
    """

    # import data
    data_dir = mfdnres.res.res_file_directory("pfasano","mfdn","pjf0015")
    res_format = "mfdn_v15"
    filename_format="mfdn_format_7_ho"
    mesh_data = mfdnres.res.slurp_res_files(data_dir,res_format,filename_format,verbose=False)
    print("directory {}".format(data_dir))
    print("full mesh {}".format(len(mesh_data)))
    ## print([mesh_point.params for mesh_point in mesh_data])

    # select and sort
    ##mesh_data = selected_mesh_data(mesh_data,(("interaction","JISP16"),("natural_orbital_iteration",0)))
    nuclide = (3,4)
    interaction_coulomb = ("JISP16",1)
    selector = {"nuclide":nuclide,"interaction":interaction_coulomb[0],"coulomb":interaction_coulomb[1],"natural_orbital_iteration":0}
    mesh_data = mfdnres.analysis.selected_mesh_data(mesh_data,selector,verbose=False)
    print("selected mesh {}".format(len(mesh_data)))
    mesh_data = mfdnres.analysis.sorted_mesh_data(mesh_data,KEY_DESCRIPTOR_NMAX_HW)
    ##print([mesh_point.params for mesh_point in mesh_data])
    
    # make energy table
    qn = (1.5,1,1)
    make_energy_table_file(mesh_data,nuclide,interaction_coulomb,qn)

    # make radius table
    qn = (1.5,1,1)
    make_radius_table_file(mesh_data,nuclide,interaction_coulomb,"rp",qn)

    # make angular momentum table
    qn = (1.5,1,1)
    make_am_table_file(mesh_data,nuclide,interaction_coulomb,qn)

################################################################
# main
################################################################

make_tabulations_pjf0007()
make_tabulations_pjf0015()

