import mfdnres
## import mfdnres.spncci_analysis

## data=mfdnres.spncci_analysis.spncci_slurp(directory)
## mfdnres.spncci_analysis.e_vs_hw(data, output_file_name, output_graph_name)

# Required environment configuration:
#
# setenv GROUP_HOME /afs/crc.nd.edu/group/nuclthy

## directory = '/home/mcaprio/results/mcaprio/spncci/runmac0415/results'
slurp_directory = mfdnres.res.res_file_directory("mcaprio","spncci","mac0415")
mesh_data = mfdnres.res.slurp_res_files(slurp_directory,"spncci",verbose=True)

results_dict = mfdnres.analysis.make_results_dict(mesh_data,("Nsigmamax","Nmax","hw"),verbose=True)

