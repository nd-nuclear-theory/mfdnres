"""spncii Parser Revision 1
    Julie Butler
    June 21, 2017
    Python3

    The controlling file for the spncci results file
    parser.  The parser used here is the method
    make_dict_spncci, stored in the file make_dict.py
    (it is imported below).  Currently this file is set up
    to parse the file whos name is provided below and make
    a graph and table of energy versus hw.  This can be
    expanded later to include more analysis.

    Once the format of the results file is on a more
    finalized form, this should be rearranged into a
    class structure for ease of access to variables.

TODO:

   Needed fields:

     need to count (J,g) ->  num e-states
       can extract by tallying from [Energies]
       but maybe SpNCCI should just output table of num actual eigenstates
         by (J,g) as min::(dimension,num_eigenstates)



     params: key-value pairs from -- Space, Interaction, Calculation
     (i.e., hw) (if not counting run)

     basis: SpJ listing, BabySpNCCI listing

     For actual calculation runs:

         energies

         state.amplitudes [from Decompositions: Nex]
           => state.decomposition_Nex

         state.decompositin_baby_spncci  [from Decompositions: BabySpNCCI]

         [rms radii] -- from diagonals of the r^2 observable
            can store in tbo until split out to more appropriate name?
            no, maybe treat as a transition observable from the get-go

         [quadrupole transitions] -- actually, don't use the MFDn
         transitions structure here, either...

         observables matrices -- store as matrices between
           (observable_filename,(J,g)_final,(J,g)_initial)  [from Observables::filenames]

         Modify accessors accordingly for spncci!

"""


# Imports the parser from the make_dict file
from mfdnres.make_dict import make_dict_spncii

# Move this too
"""
# Imported to allow for graphing
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
"""

def spnciiParser(file_name):
    """
        Arguments:
            file_name: a string.  The filename of the results
            file that needs to be analyzed.
        Returned:
            results_dict: a dictionary.  What is returned
            from make_dic_spncci.

        This method currently takes in a file names and
        converts the file contents into a list of strings,
        where each string is a line from the file.  It then
        takes the list and passes it to the parser
        make_dict_spncci.  The dictionary returned from the
        parser is then printed and returned by the method.
        This method will be expanded and possibly incoporated
        into a class structure when the format of the results
        files is finalized.
    """
    # Stores the results from inputting the file
    results = []

    # Inputs the file passed as an argument and stores the
    # parsed file contents in results
    with open(file_name, 'rt') as fin:
        for row in fin:
            results.append(row)

    # Passes results to the parser make_dict_spncci and
    # stores the returned dictionary and list
    results_dict, order = make_dict_spncii(results)

    for key,value in results_dict.items():
        print('###KEY###')
        print(key, '\n')
        print('#####VALUE#####')
        if isinstance(value,dict):
            for k,v in value.items():
                print('    ###KEY###  ', k)
                print('    #####VALUE####')
                if isinstance(key,float) or key == 'BabySpNCCI (listing)':
                    for x in v:
                        print('    ',x)
                else:
                    print('    ', v)
        else:
            print(value)
        print('\n\n') 

    # Returns the dictionary (need to remoce this in final version)
    return results_dict


## res = spnciiParser('type_specimens/runmac0405-Z3-N3-Nsigmamax02-Nmax02.res')
## 
## # Move this to its own file with a class of spncci caculation methods
## """
## # This section of code makes the graph of energy versus h bar omega.
## # It also makes a table of the form nsigmamax,nmax,hw,energy for each
## # hw value and writes those to a file
## # The list of hw values used in the file are contained in the hw
## # variable of the Mesh section of the results file.  This line
## # retrieves those values as a list
## hws = res['Mesh']['hw']
## hws = [float(i) for i in hws]
## 
## # Retrieves the energies based on the hw value (which is the key)
## # and stores them in the list energies.
## energies = []
## for x in hws:
##     energies.append(res[x][0][1])
## 
## # Retrieves the value of the variables nmax and nsigmamax, both contained
## # in the section 'Space'
## nmax =  res['Space']['Nmax']
## nsigmamax = res['Space']['Nsigmamax']
## 
## # The names for the outputted graph and table
## graphName = 'spncciTest.png'
## output_file_name = 'spncciTest.txt'
## 
## # The lable for the graph
## legend_label = 'Nsigmamax = ' + str(nsigmamax) + '; Nmax = ' + str(nmax)
## 
## # This section creates the graph and table and saves them as the file names
## # specified above, if there are the same number of elements in hws and energies
## if len(hws) == len(energies):
##     # Formats the data to be written to the file and then writes the data to the
##     # file with a lable, denoted by the hashtag.  This label can be removed if
##     # needed.
##     to_file = []
##     for i in range(len(hws)):
##         temp = str(nsigmamax)  + ',' + str(nmax) + ',' + str(hws[i]) + ',' + str(energies[i]) + '\n'
##         to_file.append(temp)
##     with open(output_file_name, 'w') as fout:
##         fout.write('#nsigmamax,nmax,hw,Energy\n')
##         for x in to_file:
##             fout.write(x)
## 
##     # This section created the graph and saves it to the file name specified
##     # above.  To display the graph instead of saving it, change the last line
##     # to 'plt.plot()' (no quotes).
##     plt.figure()
##     plt.plot (hws, energies, 'r-', label=legend_label)
##     plt.xlabel('h-bar omega (MeV)')
##     plt.ylabel('Ground State Energy (Mev)')
##     plt.legend(bbox_to_anchor=(1, 1))
##     plt.savefig(graph_name, bbox_inches='tight')
## """
