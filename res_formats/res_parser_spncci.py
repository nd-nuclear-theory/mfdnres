"""
    spncii Parser Revision 1
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
"""


import mfdnres.res

# Imports the parser from the make_dict file
from make_dict import make_dict_spncii


def res_parser_spncci(file_name):
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

    # Prints the returned dictionary.  For debugging purposes
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


# Register the parser
mfdnres.res_format('spncci', res_parser_spncci)
res = spnciiParser('type_specimens/runmac0405-Z3-N3-Nsigmamax02-Nmax02.res')
