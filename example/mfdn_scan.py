""" mfdn_scan.py

   Example convergence analysis of MFDn runs.

   10/12/17 (mac): Created, based on analysis/spncci::convergence_analysis.py.

"""

import mfdnres

# standard sorting and tabulation key
KEY_DESCRIPTOR_NMAX_HW = (("Nmax",int),("hw",float))
KEY_DESCRIPTOR_NORB_NMAX_HW = (("natural_orbital_iteration",int),("Nmax",int),("hw",float))


def make_energies(mesh_data,nuclide,interaction_coulomb,qn):
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

    energy_table = mfdnres.analysis.make_energy_table(mesh_data,KEY_DESCRIPTOR_NMAX_HW,qn)
    mfdnres.tools.write_table(
        output_file_name,"{:2d} {:7.3f} {:7.3f}",
        energy_table
    )

def make_tabulations_pjf0015():
    """
    """

    # import data
    data_dir = mfdnres.res.res_file_directory("pfasano","mfdn","pjf0015")
    res_format = "mfdn_v15"
    filename_format="mfdn_format_7_ho"
    mesh_data = mfdnres.res.slurp_res_files(data_dir,res_format,filename_format,verbose=False)
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
    
    # make energy tables
    qn = (1.5,1,1)
    ##print(mesh_data[0].params)
    ##print(mesh_data[0].energies)
    make_energies(mesh_data,nuclide,interaction_coulomb,qn)

    ## # make radius tables
    ## radius_table = mfdnres.analysis.make_radius_table(mesh_data,KEY_DESCRIPTOR_NNHW,"r",(1.0,0,1))
    ## mfdnres.tools.write_table(
    ##     "data/data-Z3-N3-JISP16-0-scan-r-1-0-1.dat","{:2d} {:2d} {:7.3f} {:7.3f}",
    ##     radius_table
    ## )
    ## 
    ## # make B(E2) tables
    ## #
    ## # Beware this is the "matter" B(E2).  In isoscalar case, divide by 4 for proton B(E2).
    ## rtp_table = mfdnres.analysis.make_rtp_table(mesh_data,KEY_DESCRIPTOR_NNHW,"Qintr",(1.0,0,1),(3.0,0,1))
    ## mfdnres.tools.write_table(
    ##     "data/data-Z3-N3-JISP16-0-scan-be2-1-0-1-3-0-1.dat","{:2d} {:2d} {:7.3f} {:7.3f}",
    ##     rtp_table

make_tabulations_pjf0015()
## def make_tabulations_3_3():
##     """
##     """
## 
##     # import data
##     slurp_directory = mfdnres.res.res_file_directory("mcaprio","spncci","mac0425")
##     mesh_data = mfdnres.res.slurp_res_files(slurp_directory,"spncci",verbose=False)
##     mesh_data = mfdnres.analysis.sorted_mesh_data(mesh_data,KEY_DESCRIPTOR_NNHW)
## 
##     # make energy tables
##     energy_table = mfdnres.analysis.make_energy_table(mesh_data,KEY_DESCRIPTOR_NNHW,(1.0,0,1))
##     mfdnres.tools.write_table(
##         "data/data-Z3-N3-JISP16-0-scan-e-1-0-1.dat","{:2d} {:2d} {:7.3f} {:7.3f}",
##         energy_table
##     )
## 
##     # make radius tables
##     radius_table = mfdnres.analysis.make_radius_table(mesh_data,KEY_DESCRIPTOR_NNHW,"r",(1.0,0,1))
##     mfdnres.tools.write_table(
##         "data/data-Z3-N3-JISP16-0-scan-r-1-0-1.dat","{:2d} {:2d} {:7.3f} {:7.3f}",
##         radius_table
##     )
## 
##     # make B(E2) tables
##     #
##     # Beware this is the "matter" B(E2).  In isoscalar case, divide by 4 for proton B(E2).
##     rtp_table = mfdnres.analysis.make_rtp_table(mesh_data,KEY_DESCRIPTOR_NNHW,"Qintr",(1.0,0,1),(3.0,0,1))
##     mfdnres.tools.write_table(
##         "data/data-Z3-N3-JISP16-0-scan-be2-1-0-1-3-0-1.dat","{:2d} {:2d} {:7.3f} {:7.3f}",
##         rtp_table
##     )

## def make_radii(r,Z,N,n,J_values,slurp_directories):
##     """
##     """
## 
##     # import data
##     #slurp_directory = mfdnres.res.res_file_directory("mcaprio","spncci","mac0423")
##     mesh_data = mfdnres.res.slurp_res_files(slurp_directories,"spncci",verbose=False)
## 
##     #J_values=[0.0,2.0,4.0]#,6.0,8.0]
##     for J in J_values:
##         # define output filename
##         output_file_name="data/data-Z{:02d}-N{:02d}-JISP16-0-scan-r-{:04.1f}-0-{:02d}.dat".format(Z,N,J,n)
##         # make radius tables
##         radius_table = mfdnres.analysis.make_radius_table(mesh_data,KEY_DESCRIPTOR_NNHW,r,(J,0,n))
##         mfdnres.tools.write_table(
##             output_file_name,"{:2d} {:2d} {:7.3f} {:7.3f}",
##             radius_table
##         )
## 
## def make_be2s(Z,N,nf,ni,J_values,slurp_directories):
##     """
##     """
## 
##     # import data
##     #slurp_directory = mfdnres.res.res_file_directory("mcaprio","spncci","mac0423")
##     mesh_data = mfdnres.res.slurp_res_files(slurp_directories,"spncci",verbose=False)
## 
##     for Ji in J_values:
##         for Jf in J_values:
##             if Jf > Ji :
##                 continue
##             if (Jf+Ji)<2.0 or abs(Jf-Ji)>2.0:
##                 continue
##             output_file_name="data/data-Z{:02d}-N{:02d}-JISP16-0-scan-be2-{:04.1f}-0-{:02d}-{:04.1f}-0-{:02d}.dat".format(Z,N,Jf,nf,Ji,ni)
##             # make B(E2) tables
##             #
##             # Beware this is the "matter" B(E2).  In isoscalar case, divide by 4 for proton B(E2).
##             rtp_table = mfdnres.analysis.make_rtp_table(mesh_data,KEY_DESCRIPTOR_NNHW,"Qintr",(Jf,0,nf),(Ji,0,ni))
##             mfdnres.tools.write_table(
##                 output_file_name,"{:2d} {:2d} {:7.3f} {:7.3f}",
##                 rtp_table
##             )
