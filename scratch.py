def round_half_integer(x):
    """ Rounds number to nearest floating point half integer.

    Args:
        x (number): quantity to be rounded

    Returns:
        float: number rounded to (generally exact) integer or half integer
    """
    
    return round(2*float(x))/2


    
    
def extract_tabular_data(tokenized_lines,dtype):
    """Parse lines as numpy structured array.

    DEPRECATED in favor of using np.array directly.

    The data types specified in dtype are used for string conversion,
    so they must be actual Python types which support string
    conversion, rather than the more general numpy data type
    specifications.  Thus, e.g., the numpy tutorial example for structured
    arrays would not be useable:

        dtype=[('foo', 'i4'),('bar', 'f4'), ('baz', 'S10')])  # can't use here

    But the following example would work:

        dtype=[('foo', int),('bar', float), ('baz',str)])

    Arguments:
        tokenized_lines (iterator): tokenized input lines
        dtype (list of tuple): list of (tag,conversion) pairs

    Returns:
        (array): structured array representation of input table
    """

    return np.array(tokenized_lines,dtype=dtype)
