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
"""


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


def make_dict_spncii(results):
    """
        Arguments:
            results: a list of strings.  The results of importing a results file.
                Each string in results was a line in the results file.
        Returned:
            results_dic: a dictionary.  This dictionary holds the parsed data.
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
                data type, designed to deal with the multiple 'Calculations' and
                'Energies' section headings.  It does so by making the key the
                 value of hw for that particular section.  The value is a list
                 of tuples, with each tuple having the form ((J, gex, i,), E).
                 There should be one type per line of data in the energies
                 sections.
            order: a list of tuples.  This list contains typles of the type
                (headerName, dataType) where headerName is the section
                label and dataType is from the class Data. This list is used
                to keep the order the data sections appeared in the results
                file since they will get scrambled in the dictionary

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
    # This is the main loop of the method.  It goes through the line section by section.
    while len(results) > 0:
        # Used to store each section of data as it is being parsed
        temp_array = []
        # These booleans could be set to true later if the section being analyzed is
        # determined to be the data tyoe indicated by the boolean name.
        isVariable = False
        isSpncciCalculation = False
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

            # This if state if one of the differences between the mfdn15
            # parser and the spncci parser.  This section check to see
            # if the next heading is 'Energies', and it it is, it adds
            # all the data under that section to temp_array.  This is done
            # to keep the energy data together with its corresponding
            # hw value, stored under the preceeding 'Calculations' label.
            if results[0][1:9] == 'Energies':
                temp_array.append(results[0])
                results.pop(0)
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
            # in value, so it is variable data.  But it also adds a check to
            # make sure that the section being anlyzed is not 'Calculations'
            # which is handled differently
            elif '=' in value[0] and (not key == 'Calculation'):
                order.append((key, Data.variable))
                isVariable = True
            # This elif filters out the 'Calculations' sections.
            elif value[0][0] == 'hw':
                order.append((key, Data.spncci_calculation))
                isSpncciCalculation = True
            # The else case assumes that if the data does not fall into
            # one of the above categories, then it is row data.
            else:
                order.append((key, Data.row))

            # If the data in the section is determined to be variable data,
            # then values is converted from a lists of lists to a dictionary
            # with the variable name as the key and the value of the variable
            # as the value
            if isVariable:
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
            # is the value of hw (as a float) and the value is a list of tuples.
            # Each tuple  correspnds to a line of data and has the form
            # ((J, gex, i), E)
            if isSpncciCalculation:
                key = float(value[0][2])
                value.pop(0)
                value.pop(0)
                temp = []
                for x in value:
                    temp1 = (float(x[0]), float(x[1]), float(x[2]))
                    temp2 = float(x[3])
                    temp.append((temp1, temp2))
                value = temp

            # Adds the formatted key and value to results_dict
            results_dict[key] = value

    return results_dict, order


def make_dict_mfdn15(results):
    """
        Arguments:
            results: a list of strings.   The results of importing a results file.
                Each string in results was a line in the results file.
        Returned:
            results_dict: a dictionary. This dictionary holds the parsed data.
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
            order: a list of tuples. This list contains typles of the type
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
        isVariable = False
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
                isVariable = True
            # The else case assumes that if the data does not fall into
            # one of the above categories, then it is row data.
            else:
                order.append((key, Data.row))

            # If the data in the section is determined to be variable data,
            # then values is converted from a lists of lists to a dictionary
            # with the variable name as the key and the value of the variable
            # as the value
            if isVariable:
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
