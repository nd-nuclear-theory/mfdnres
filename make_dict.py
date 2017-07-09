"""
    make_dict.py
    Julie Butler
    June 22, 2017
    Python3

    This file parses results files from either MFDn15 or spncci.
    The two methods contained within are:
        make_dict_spncci
        make_dict_mfdn15
    Both methods take as an argument a list of strings, where each
    string corresponds to a line of the results file.  Both methods
    return a dictionary and a list.  The dictionary contains the
    parsed data separated by section label (more details in the
    docstrings of the methods).  This list contains tuples which
    specify the order section label appeared in the file and the
    type of data each contained.

    This file also contains the class Data, which uses Enum to
    specify the data types of the various sections within the
    results files.

    Classes:
        Data: Defines the different types of data found in results
            files.
            Does so through the use of the import enum.
    Methods:
        make_dict_spncci: Parser for SpNCCI results files.  Takes
            as an argument the list of strings resulting from
            importing the results file.  Returns a dictionary of the
            parsed data and a list denoting the order of headings in the
            results file.
        spncci_results_section: Helper method for make_dict_spncci.
            Takes as an argument a list of list, which is the data 
            contained within a single results heading of a SpNCCI
            results file (without the headings 'Calculations' or 
            'Results', and without the hw = line). Returns a 
            dictionary with the formatted data for the results section.
        make_dict_mfdn15: Parser for MFDn Version 15 results files.
            Takes as an argument the list of strings resulting from
            importing the results file.  Returns a dictionary of the
            parsed data and a list denoting the order of headings in the
            results file.
            
"""

# Used in Data
import enum


@enum.unique
class Data(enum.Enum):
    """
        This class defines the different types of data that are
        expected in mfdn15 and spcci results files.

        super_header: a super header is a header directly above
        another header:
            [Super header]
            [Header]
        It contains no data of its own but contains one or more
        data sections.

        variable: Varaible data is data contained within a section
        that has the following form:
            VariableName = Value
        Variable data must include an equals sign to be identified
        correctly.

        row: Row data is data contained within a section that has
        the following form:
            x    x    x    x   x
        where x is a number.  The number of spaces between entries
        in a row does not matter.

        spncci_calculation: This data type is contained under a
        'Calculations' heading in a spncci results file.  It
        contains an line of the form 'hw = x' where x is a number,
        the section header 'Energies', and tabulated energy data.
    """
    super_header = 0
    variable = 1
    row = 2
    spncci_calculation = 3


# Unicode representation of the number sign.  Used later in the
# program to remove rows that begin with the number sign.  Can
# not use the character because it causes a comment.  Used to
# remove the commented lines
# from the code.
num_sign = u"\u0023"

# Really?  I don't see the problem.
this_is_also_a_num_sign = "#"

def make_dict_spncci(results):
    """
        Arguments:
            results (list of strings:) The results of importing a results file.
                Each string in results was a line in the results file.
        Returned:
            results_dic (dictionary):  This dictionary holds the parsed data.
                There are four types of entries possible in this dictionary,
                corresponding to the four data types defined in the class
                Data.  If a super header is found, the name of the super
                header is entered into the dictionary as the key and the
                value is the string 'None'.  Variable data is entered as a
                nested dictionary, where the main key is the section label.
                The value is a dictionary of the variable data where the keys
                are the variable names and the values are the values of the
                variables.  Row data is entered with the section header as the
                key and a nested list as the value.  Each of the inner list
                contains one line of data from the results file.  The string
                as separated into elements of a list using white space.  The
                final type of entry into results_dict is from the spncci_calculations
                data type, designed to deal with the multiple 'RESULTS' sections
                in the results file.  It deals it the Results sections by creating
                nested dictionary entries in results_dict.  The main key is the 
                value of hw for that Results section (as a floating point number).
                The value is then a dictionary formatted with the keys being 
                section headings and the values being the data in that section.
            order (list of tuples): Denotes the order headings were located in the file.
                This list contains typles of the type (headerName, dataType)
                where headerName is the section label and dataType is from the
                class Data. This list is used to keep the order the data sections
                appeared in the results file since they will get scrambled in the
                dictionary

        This method parses a spncci results file.  It takes in a list of strings,
        where each string corresponds to a line of the results file.  It returns
        the parsed data in a dictionary.  It also returned a list of tuples for
        information of the order of the data sections and the type of data
        each section contains.
    """
    # These two items will hold the parsed data.
    results_dict = {}
    order = []
    # This removes leading white space from all the elements of results, making it
    # easier to find the [Header name] notation of section headers
    for i in range(0, len(results)):
        results[i] = results[i].lstrip()
    # This section removes any blank or commented lines from results, which are not
    # needed in the final result
    i = 0
    while i < len(results):
        if results[i] == '' or results[i][0] == num_sign:
            results.pop(i)
        else:
            i = i + 1

    # This section of code creates a list of all of the headings under the first
    # occurence of the super header '[RESULTS]'.  The headings under all other
    # occurences of '[RESULTS]' should be the same.
    headers_under_results = []
    found = False
    index = 0
    while not found:
        if len(results[index]) > 7 and results[index][0:9] == '[RESULTS]':
            found = True
        else:
            index = index + 1
    first_occurence_results = index
    index = first_occurence_results + 1
    while index < len(results) and (not results[index][0:9] == '[RESULTS]'):
        if results[index][0] == '[':
            headers_under_results.append(results[index])
        index = index + 1

    # This is the main loop of the method.  It goes through the line section by section.
    while len(results) > 0:
        # Used to store each section of data as it is being parsed
        temp_array = []
        # These booleans could be set to true later if the section being analyzed is
        # determined to be the data tyoe indicated by the boolean name.
        is_variable = False
        is_spncci_results = False
        # This section deals with the unlikely possibility that there is some data
        # at the beginning of the file that is not contained in a section header.
        # It deals with this by removing entries from results until it finds a
        # section header, denoted by a '['.
        if not results[0][0] == '[':
            print ('Beginning data not contained in a section.  Removing it')
            while not results[0][0] == '[':
                results.pop(0)
        else:
            # The first enntry in results is a section header
            temp_array.append(results[0])
            results.pop(0)
            # This loops through results, moving data to temp_array
            # and removing it from results until either the end
            # of the file is reached of a section header is found.
            while len(results) > 0 and (not results[0][0] == '['):
                temp_array.append(results[0])
                results.pop(0)

            # This is one of the differences between the mfdn15
            # parser and the spncci parser.  This section checks to see
            # if the next heading is one of the known headings under '[RESULTS]',
            # and it it is, it adds all the data under that section to temp_array.
            #  This is done to keep the data from each mesh point together.
            if results[0] in headers_under_results and temp_array[0][0:9] == '[RESULTS]':
                temp_array.append(results[0])
                results.pop(0)
                while len(results) > 0 and (not results[0][0:9] == '[RESULTS]'):
                    temp_array.append(results[0])
                    results.pop(0)

            # This loop converts the individual strings into list, where
            # the elements are separated by white space. It excludes the
            # first element of temp_array, which is the section label.  This
            # should remain a string because it will be used as a key
            # in the dictionary.
            for i in range(1, len(temp_array)):
                temp_array[i] = [str(j) for j in temp_array[i].split()]

            # This line makes key the first element of temp_array, which is the
            # section header, but it removes the beginning and ending brackets
            key = temp_array[0][1:-2]

            # This section deals with the case that there was no data following
            # a section header, probably because it was a super header.  If there
            # is not data remaining in temp_array after the section header.  The
            # value in the dictionary is set to the string 'None'.
            if len(temp_array[1:]) > 0:
                value = temp_array[1:]
            else:
                value = 'None'

            # This section determines what type of data is contained in each
            # section be analyzing value
            # The if handles the cause where value is 'None' which means it is
            # most likely a super header.
            if value == 'None':
                order.append((key, Data.super_header))
            # This elif deals with the case there an equal sign is contained
            # in value, so it is variable data.  But it also adds a check to
            # make sure that the section being anlyzed is not 'Calculations'
            # which is handled differently
            elif '=' in value[0] and (not key == 'RESULTS'):
                order.append((key, Data.variable))
                is_variable = True
            # This elif filters out the 'Calculations' sections.
            elif key == 'RESULTS':
                order.append((key, Data.spncci_calculation))
                is_spncci_results = True
            # The else case assumes that if the data does not fall into
            # one of the above categories, then it is row data.
            else:
                order.append((key, Data.row))

            # If the data in the section is determined to be variable data,
            # then values is converted from a lists of lists to a dictionary
            # with the variable name as the key and the value of the variable
            # as the value
            if is_variable:
                temp = {}
                for x in value:
                    if len(x) == 3:
                        temp[x[0]] = x[2]
                    elif len(x) > 3:
                        temp[x[0]] = x[2:]
                    else:
                        temp[x[0]] = 'None'
                value = temp

            # If the data in the section is determined to be spincci calculation
            # then the list of list is converted into a dictionary where the key
            # is the value of hw (as a float) and the value is a dictionary where
            # the keys are the section headings and the values are the data within 
            # those sections.
            if is_spncci_results:
                key = float(value[1][2])
                value.pop(0)    # Removes '[Calculation]' 
                value.pop(0)    # Removes 'hw = x' (data stored as the key)
                value = spncci_results_section (value)

            # Adds the formatted key and value to results_dict
            results_dict[key] = value

    return results_dict, order

def spncci_results_section (results):
    """
        Arguments:
            results (list of list)
        Returned:
            results_dict (dictionary)

        Reformats the value of results_dict for a
        SpNCCI Results sections.
    """
    for i in range(0,len(results)):
        if results[i][0][0] == '[':
            results[i] = ' '.join(results[i])
    results_dict = {}
    while len(results) > 0:
        key = results[0][1:-1]
        results.pop(0)
        value = []
        while len(results) > 0 and (not results[0][0] == '['):
            value.append(results[0])
            results.pop(0)
        results_dict[key] = value
    return results_dict    


def make_dict_mfdn15(results):
    """
        Arguments:
            results (list of strings):  The results of importing a results file.
                Each string in results was a line in the results file.
        Returned:
            results_dict (dictionary): This dictionary holds the parsed data.
                There are three types of entries possible in this dictionary
                corresponding to the four data types defined in the class
                Data.  If a super header is found, the name of the super
                header is entered into the dictionary as the key and the
                value is the string 'None'.  Variable data is entered as a
                nested dictionary, where the main key is the section label.
                The value is a dictionary of the variable data where the keys
                are the variable names and the values are the values of the
                variables.  Row data is entered with the section header as the
                key and a nested list as the value.  Each of the inner list
                contains one line of data from the results file.  The string
                as separated into elements of a list using white space.
            order (list of tuples): Denotes the order headings were located in
                the results file. This list contains typles of the type
                (headerName, dataType) where headerName is the section
                label and dataType is from the class Data. This list is used
                to keep the order the data sections appeared in the results
                file since they will get scrambled in the dictionary

        Takes in a list of strings (the lines from the results file) and formats them
        into a dictionary.  The keys are the section headers and the values are the
        data in that section.  The format of the values in the dictionary is determined
        by the type of data in the section being analyzed.  This method also returns a
        list of tuples of the form (header,type of data) where header is the section label
        and type of data is from the class Data.
    """
    # These two items will hold the parsed data.
    results_dict = {}
    order = []

    # This removes leading white space from all the elements of results, making it
    # easier to find the [Header name] notation of section headers
    for i in range(0, len(results)):
        results[i] = results[i].lstrip()

    # This section removes any blank or commented lines from results, which are not
    # needed in the final result
    i = 0
    while i < len(results):
        if results[i] == '' or results[i][0] == num_sign:
            results.pop(i)
        else:
            i = i + 1

    # This is the main loop of the method.  It goes through the line section by section.
    while len(results) > 0:
        # Used to store each section of data as it is being parsed
        temp_array = []
        # This boolean could be set to true later if the section being analyzed is
        # determined to be variable data.
        is_variable = False
        # This section deals with the unlikely possibility that there is some data
        # at the beginning of the file that is not contained in a section header.
        # It deals with this by removing entries from results until it finds a
        # section header, denoted by a '['.
        if not results[0][0] == '[':
            print ('Beginning data not contained in a section.  Removing it')
            while not results[0][0] == '[':
                results.pop(0)

        else:
            # The first entry in results is a section header
            temp_array.append(results[0])
            results.pop(0)
            # This loops through results, moving data to temp_array
            # and removing it from results until either the end
            # of the file is reached of a section header is found.
            while len(results) > 0 and (not results[0][0] == '['):
                temp_array.append(results[0])
                results.pop(0)

            # This loop converts the individual strings into list, where
            # the elements are separated by white space. It excludes the
            # first element of temp_array, which is the section label.  This
            # should remain a string because it will be used as a key
            # in the dictionary.
            for i in range(1, len(temp_array)):
                temp_array[i] = [str(j) for j in temp_array[i].split()]

            # This line makes key the first element of temp_array, which is the
            # section header, but it removes the beginning and ending brackets
            key = temp_array[0][1:-2]

            # This section deals with the case that there was no data following
            # a section header, probably because it was a super header.  If there
            # is not data remaining in temp_array after the section header.  The
            # value in the dictionary is set to the string 'None'.
            if len(temp_array[1:]) > 0:
                value = temp_array[1:]
            else:
                value = 'None'

            # This section determines what type of data is contained in each
            # section be analyzing value
            # The if handles the cause where value is 'None' which means it is
            # most likely a super header.
            if value == 'None':
                order.append((key, Data.super_header))
            # This elif deals with the case there an equal sign is contained
            # in value, so it is variable data.
            elif '=' in value[0]:
                order.append((key, Data.variable))
                is_variable = True
            # The else case assumes that if the data does not fall into
            # one of the above categories, then it is row data.
            else:
                order.append((key, Data.row))

            # If the data in the section is determined to be variable data,
            # then values is converted from a lists of lists to a dictionary
            # with the variable name as the key and the value of the variable
            # as the value
            if is_variable:
                temp = {}
                for x in value:
                    if len(x) == 3:
                        temp[x[0]] = x[2]
                    else:
                        temp[x[0]] = 'None'
                value = temp

            # Formatted key and value are placed in results_dic
            results_dict[key] = value
    return results_dict, order
