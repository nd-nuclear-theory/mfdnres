import mfdnres.spncci_analysis

directory = '/home/mcaprio/results/mcaprio/spncci/runmac0415/results/*.res'
output_file_name = 'spncci_analysis_test_text.txt'
output_graph_name = 'spncci_analysis_test_graph.png'

## data=mfdnres.spncci_analysis.spncci_slurp(directory)
## mfdnres.spncci_analysis.e_vs_hw(data, output_file_name, output_graph_name)

mesh_data = mfdnres.res.slurp_res_files(directory,"spncci")

results_dict = mfdnres.analysis.make_results_dict(mesh_data,("Nsigmamax","Nmax","hw"),verbose=True)

