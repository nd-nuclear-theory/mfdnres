""" tools.py -- internal tools for use by data file parsers

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame
    5/31/15 (mac): Initiated (as mfdn_res.py).
    6/5/15 (mac): Restructure as subpackage.
    7/26/15 (mac): Allow mismatch in line parser.
    7/8/17 (mac): Add write_lines and write_table.
    7/9/17 (mac): Add parsing tools for structured results files.
    
"""

import re

################################################################
# line parser for free-form res files
################################################################

def parse_line(line,pattern,strict=True):
    """Matches given input line to regex.

    Line is stripped of leading/trailing whitespace.

    Args:
        line (string) : input line to parse
        pattern (string) : string containing regex
        strict (bool) : raise exception if mismatch (default: True)

    Returns:
        The match object obtained from match.

    Raises:
        ValueError if input line does not match given regex.

    """

    # read and strip string
    in_str = line.strip()

    # attempt match
    match = re.match(pattern,in_str)
    if (match is None):
        print("Expected:",pattern)
        print("Read:",in_str)
        raise(ValueError("unexpected string"))
    
    return match

################################################################
# parsing for structured files
################################################################

#    After tokenization on whitespace, the following special line
#    types are recognized:
#    
#       - []: Empty line -- ignored.
#       - ["#",...]: Comment line -- ignored.
#       - ["[section]"]: Section header, where <section> denotes an
#         arbitrary string.

def is_active_line(tokens):
    """ Identify nonempty, noncomment line.

    Helper function for file parsing.

    Arguments:
        (list of str): tokenized line

    Returns:
        (bool)
    """
    return bool(tokens) and (tokens[0]!="#")

section_header_regex = re.compile(r"\[(.*)\]")

def extract_section_name(tokens):
    """ Identify header line and (if header line) extract section name.

    Helper function for file parsing.
    
    Sections names can contain whitespace, but whitespace will be regularized
    to single spacing, as an artifact of tokenization followed by rejoining.

    Arguments:
        (list of str): tokenized line

    Returns:
        (str): Header line name (if header) else None.
    """

    # Beware that, after tokenization, the section name may have been split
    # over multiple tokens:
    #
    #      ('[BabySpNCCI', '(listing)]')

    # test for header (with short circuit trap for empty line)
    is_header_line = bool(tokens) and (tokens[0][0]=="[") and (tokens[-1][-1]=="]")

    if (is_header_line):
        spliced_line=" ".join(tokens)
        return section_header_regex.match(spliced_line).group(1)  # or just chop off the bracket on each end...
    else:
        return None
    

def split_and_prune_lines(lines):
    """Split input lines into tokenized lines, suppressing comment
    (beginning with hash token) or empty lines.

    Tokenization is by whitespace, i.e., with split.

    The tokenized lines are returned as tuples rather than lists to
    support downstream processing.  For instance, structured array
    creation with np.array requires the entries to be tuples.

    Arguments:
       (iterable of str): input lines

    Returns:
       (iterator of tuple of str): split and filtered lines

    """
    
    tokenized_lines = map(lambda s : tuple(s.split()),lines)
    return filter(is_active_line,tokenized_lines)

def extracted_sections(tokenized_lines):
    """Provide iterator yielding succesive sections from given tokenized lines.

    This is a generator function.

        >>> tokenized_lines = map(lambda s : s.split(),["[A]","a b","c","[D]","[E]"])
        >>> list(section_generator(tokenized_lines))
        [('A', [['a', 'b'], ['c']]), ('D', []), ('E', [])]

    """
    
    # convert lines to iterator
    #
    # This ensures that the lines are represented as an iterator, not
    # just an interable, so that we can call next on them.
    tokenized_line_iterator = iter(tokenized_lines)

    # extract first header line (to "prime" the loop)
    header_tokens = next(tokenized_line_iterator,None)
    section_name = extract_section_name(header_tokens)
    if (not section_name):
        raise ValueError("expected section header line but found {}".format(header_tokens))

    # loop over sections
    #
    # Termination is when section_name==None.
    while (section_name):
        
        # accumulate non-header lines
        #
        # Note: We could almost use itertools.takewhile with
        # "is_not_section_header_line" as predicate function, but this
        # would discard the next header line, or by adding an extra
        # filtering layer which inserts an end-of-section flag into
        # the iteration over lines.
        section_lines = []
        done_with_section = False
        while (not done_with_section):
            # get line
            line_tokens = next(tokenized_line_iterator,None)

            # check for next section header or end of input
            next_section_name = extract_section_name(line_tokens)
            done_with_section = bool(next_section_name) or (not line_tokens)

            # if line is regular line: append it to this section
            if (not done_with_section):
                section_lines.append(line_tokens)

        # yield section
        yield (section_name,section_lines)

        # advance to next section
        section_name = next_section_name

################################################################
# key-value conversion
################################################################

def singleton_of(conversion):
    """Generate conversion function to convert a single-entry list of
    strings to a single value of given type.

    >>> a = ["1"]
    >>> singleton_of(int)(a)
        1

    Arguments:
        conversion (function): type converstion function for single entry
    
    Returns:
        (function): function extract such entry

    """

    def f(data):
        if (len(data)!=1):
            raise ValueError("expecting list of length 1 but found {}".format(data))
        return conversion(data[0])

    return f

def list_of(conversion):
    """Generate conversion function to convert list of strings to list
    of given type.

    >>> a = ["1","2"]
    >>> list_of(int)(a)
        [1,2]

    Arguments:
        conversion (function): type converstion function for single entry
    
    Returns:
        (function): function to convert list of such entries

    """

    def f(data):
        return list(map(conversion,data))
    return f

def extract_key_value_pairs(tokenized_lines,conversions):
    """ Parse tokenized lines as key-value pairs.
    
    A valid key-value line is of the form:

       [<key>,"=",<v1>,...]

    Values are only retained if a conversion is specified
    for that key string.

    >>> test_lines = ["a = 1","b = 1 2 3","c = 42"]
    >>> tokenized_lines = split_and_prune_lines(test_lines)
    >>> conversions = {"a" : singleton_of(int), "b" : list_of(int)}
    >>> extract_key_value_pairs(tokenized_lines,conversions)
    {'b': [1, 2, 3], 'a': 1}

    Arguments:
        tokenized_lines (iterator): tokenized input lines
        conversions (dict): conversion functions for recognized key strings

    Returns:
        (dict): key value pairs obtained through given conversions
    """

    results = dict()
    for tokenized_line in tokenized_lines:

        # validate line format
        valid_line = (len(tokenized_line)>=3) and (tokenized_line[1]=="=")
        if (not valid_line):
            raise ValueError("expected key-value line but found {}".format(tokenized_line))

        # extract line parts
        key = tokenized_line[0]
        value_strings = tokenized_line[2:]
        
        # convert and store value
        if (key in conversions):
            results[key] = conversions[key](value_strings)

    return results

################################################################
# table output
################################################################

def write_lines(filename,lines):

    """ Write lines of text to file.

    Arguments:
        filename (str): name for output file
        lines (list of str): output lines (sans newlines)
    """

    output_string = "\n".join(lines) + "\n"
    with open(filename, 'wt') as out_file:
        out_file.write(output_string)

def write_table(out_filename,format_descriptor,data,verbose=False):
    """Write output tabulation to file.
    
    Data may be in any form accessible by double indexing as
    data[row][col], e.g., a simple list of lists or a numpy array.

    Arguments:
       out_filename (str): output filename
       format_descriptor (str): format descriptor for single line
       data (array): data values
       verbose (bool, optional): verbose output flag
    """

    data_lines = [
        format_descriptor.format(*entry)
        for entry in data
    ]
    if (verbose):
        print(data_lines)
    write_lines(out_filename,data_lines)

################################################################
# range construction utilities
################################################################

def value_range(x1,x2,dx,epsilon=0.00001):
    """ value_range(x1,x2,dx) produces float or integer results [x1, ..., x2] with step dx

    This is meant to emulate Mathematica's Range[x1,x2,x], with a "<=" upper bound, in
    contrast to Python's range, which is limited to integers and has an "<" upper bound.

    Borrowed from mcscript package.
    
    Limitations: Currently assumes dx is positive.  Upper cutoff is
    slightly fuzzy due to use of epsilon in floating point comparison.

    epsilon: tolerance for cutoff, as fraction of dx
    """

    value_list = []
    x=x1
    while (x < (x2 + dx * epsilon)):
        value_list.append(x)
        x += dx
    return value_list


################################################################
# test
################################################################

if (__name__=="__main__"):

    # test structured file parsing

    test_lines = ["[A]","# comment","   ","a b","c","[D]","[E]","[A Z]"]
    print("Raw lines:",test_lines)

    active_lines_iterator = split_and_prune_lines(test_lines)
    print("Tokenized and filtered lines:",list(active_lines_iterator))

    extracted_sections_iterator = extracted_sections(split_and_prune_lines(test_lines))
    print("Sections:",list(extracted_sections_iterator))

    print()

    # test key-value conversions

    conversions = {"a" : singleton_of(int), "b" : list_of(int)}

    test_lines = ["a = 1","b = 1 2 3","c = 42"]
    print("Raw lines:",test_lines)
    tokenized_lines = split_and_prune_lines(test_lines)
    results = extract_key_value_pairs(tokenized_lines,conversions)
    print("Key-value pairs:",results)
    
    # test key-value conversions again with "bad" data
    if (False):
        test_lines = ["a = 1 2","b = 1 2 3","c = 42"]
        tokenized_lines = split_and_prune_lines(test_lines)
        results = extract_key_value_pairs(tokenized_lines,conversions)

    if (False):
        test_lines = ["not a valid line","b = 1 2 3","c = 42"]
        tokenized_lines = split_and_prune_lines(test_lines)
        results = extract_key_value_pairs(tokenized_lines,conversions)
