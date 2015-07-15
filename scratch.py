def round_half_integer(x):
    """ Rounds number to nearest floating point half integer.

    Args:
        x (number): quantity to be rounded

    Returns:
        float: number rounded to (generally exact) integer or half integer
    """
    
    return round(2*float(x))/2

