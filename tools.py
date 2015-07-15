""" tools.py -- internal tools for use by data file parsers

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame
    5/31/15 (mac): Initiated (as mfdn_res.py).
    6/5/15 (mac): Restructure as subpackage.
    Last modified 6/5/15.
    
"""

import re

################################################################
# line parser for res file input
################################################################

def parse_line(line,pattern):
    """Matches given input line to regex.

    Line is stripped of leading/trailing whitespace.

    Args:
        line: input line to parse
        pattern: string containing regex

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
