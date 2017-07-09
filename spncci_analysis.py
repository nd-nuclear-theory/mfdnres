"""
    SpNCCI Analysis
    Julie Butler
    July 5, 2017
    Python3

    Set of analysis function for SpNCCI results files.

    Methods:
    spncci_slurp: Parses a directory of SpNCCI results files.  Takes
        as an argument the directory where the results files are located.
        Returns a nested dictionary.  The outer keys of the dictionar
        are tuples of (Nmax, Nsigmamax).  The inner keys are h-bar omega
        values.  Each h-bar omega value maps to an instance of
        SpNCCIMeshPointData.  
    e_vs_hw: creates a table and graph of energy versus h-bar omega.  Takes
        as arguments a dictionary made from spncci_slurp, a filename for the
        outputted table, and a filename for the outputted graph.  Returns nothing.
        Produces a table where the lines have the format 'Nmax  Nsigmamax  hw  E',
        a saves the table to the filename specified in the arguments.  This table
        has a header line beginning with a '#', denoting what the columns in the 
        table are.  This method also produces a graph, with energy on the y-axis,
        and h-bar omega on the x-axis.  Different series on the graph are denoted
        in the legend by Nmax and Nsigmamax.  The graph is both displayed and 
        automatically saved to the filename specified in the arguments.
"""
# Allows for all files in a directory to be imported
import glob
# Imports the parser
import mfdnres.res
# Allows for graphing
##import matplotlib.pyplot as plt

def spncci_slurp (directory):
    """
        Arguments:
            directory (string): Location of the SpNCCI results files. 
                Should be of the form '/location/of/files/*.res'
        Returned:
            data (nested dictionary): Maps from (Nmax, Nsigmamax) to hw to 
                SpNCCIMeshPointData instance.  The outer keys of the dictionary
                are tuples of the form (Nmax, Nsigmamax).  The inner keys are 
                h-bar omega values.  Each h-bar omega value maps to its correspionding
                instance of SpNCCIMeshPointData.

        Takes all the SpNCCI results files from a specified directory and parses them into
        a dictionary for later analysis.  This dictionary is returned at the end of the 
        method.
    """
    # imports all the files from the specified directory
    files = glob.glob(directory)
    # Holds the parsed data
    data = {}
    for fle in files:
        # instances is a list of all the SpNCCIMeshPoint instances created from a single results
        # file
        instances = mfdnres.res.read_file(fle, res_format='spncci', verbose=True)
        # Assuming all instances of SpNCCIMeshPoint data from the same results file will have
        # identical values fo Nmax and Nsigmamax
        Nmax = instances[0].params['Nmax']
        Nsigmamax = instances[0].params['Nsigmamax']
        temp = {}
        for x in instances:
            temp[x.hw] = x
        data[(Nmax, Nsigmamax)] = temp
    # data now has the following format (Nmax, Nsigmamax): hw: instance of SpNCCIMeshPointData
    return data

def e_vs_hw(mesh_data,output_file_name):
    """
        Arguments:
            data (nested dictionary): The nested dictionary created by spncci_slurp
            output_file_name (string): The filename where the table is to be written
            output_graph_name (string): The filename where the graph is to be saved
        Returned:
            None.

        Takes the parsed data from spncci_slurp and created a table and a graph of
        energy versus h-bar omega.  Each line of the table has the form 'Nmax  Nsigmamax
        hw  energy'.  There is a header line at the beginning of the table file, identified
        by a preceding '#', that denotes the columns.  The table is automatically saved to the 
        supplied filename.  The graph created has energy on the y-axis and h-bar omega on the 
        x-axis.  Each series of the graph is identified in the legend by Nmax and Nsigmamax.  The
        graph is both displayed and automatically saved to the supplied filename.
    """
    # The contents of output_lines are written to the table file
    output_lines = ['#Nmax  Nsigmamax  hw  energy']
    for truncation, hw_dictionary in mesh_data.items():
        (Nmax, Nsigmamax) = truncation
        for hw, mesh_point_data in hw_dictionary.items():
            # The key is hw
            energies = mesh_point_data.energies
            # Finds the quantum number tuple associated with the lowest energy
            lowest_energy_qn = min(energies, key=energies.get)
            # E is set to the lowest energy
            E = energies[lowest_energy_qn]
            # Tabulation to be written to the table file
            string = "{:2d} {:2d} {:7.3} {:e}".format(Nmax,Nsigmamax,hw,E)
            output_lines.append(string)
        # How the legend of each series will appear on the graph
        ##legend_label = 'Nmax: ' + str(Nmax) + '; Nsigmamax: ' + str(Nsigmamax)

    # Creates the table
    res.tools.write_lines(output_file_name,output_lines)

