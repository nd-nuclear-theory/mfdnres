import mfdnres.spncci_analysis

directory = '/afs/crc.nd.edu/user/j/jbutler7/runmac0415/results/*.res'
output_file_name = 'spncci_analysis_test_text.txt'
output_graph_name = 'spncci_analysis_test_graph.png'

data=mfdnres.spncci_analysis.spncci_slurp(directory)
mfdnres.spncci_analysis.e_vs_hw(data, output_file_name, output_graph_name)
