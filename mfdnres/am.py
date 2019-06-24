""" am.py -- angular momentum utilities for mfdnres

    Loosely inspired by the C++ am/am.h.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    07/15/17 (mac): Created.
    06/24/19 (mac): Add parity_grade.
    
"""

import math

def hat(j):
    """ Angular momentum hat symbol.

    j-hat = (2*j+1)**(1/2)

    Arguments:
        j (int or float): angular momentum

    Returns:
        (float): j-hat
    """

    return math.sqrt(2*j+1)

def is_integer(j):
    """ Tests if half-integer j is an integer.
    """
    return not(int(2*j)%2)
  
def parity_sign(j):
    """ Parity sign corresponding to integer angular momentum.

    Arguments:
        j (int): *integer* angular momentum

    Returns:
        (int): (-)**j
    """

    # From am.h:
    #
    # DEBUGGING: Note that C++ % operator is not guaranteed positive
    # definite unless both arguments are nonnegative.  Although taking
    # abs() does not in general preserve modular equivalency class, it
    # does suffice for the present purpose (checking evenness or
    # oddness).  Without abs(), e.g., ParitySign(-1) can result in
    # failure, from a remainder result other than 0 or 1.

    if (not(is_integer(j))):
        raise(ValueError("noninteger angular momentum {}".format(j)))
    remainder = abs(int(j)) % 2
    sign = 1 - 2*remainder
    return sign

def parity_grade(parity):
    """ Recover grade from parity.

    Arguments:
        parity (int): parity +1 or -1

    Returns:
        grade (int): parity grade 0 or 1
    """

    if (not ((parity==+1) or (parity==-1))):
        raise(ValueError("invalid parity {}".format(parity)))
    return (1-parity)//2

def allowed_triangle(j1,j2,j3):
    """ Test if three angular momenta are coupled legally, i.e., they form a closed triangle.
  
    Also checks if their combined "integrity" is valid.
  
    Returns True if both triangle condition and parity condition are met, False if either is not met.

    Arguments:
        j1, j2, j3 (int or float): angular momenta
    
    Return:
        (bool): allowed coupling
    """
    
    triangular = ((abs(j1-j2) <= j3) and (j3 <= (j1+j2)))
    proper_integrity = is_integer(j1+j2+j3)
    return (triangular and proper_integrity)

def product_angular_momenta(j1,j2):
    """Create a range representing the angular momenta that j1 and j2 can
    be coupled to under the triangle inequality.

    Arguments:
        j1, j2 (int or float): angular momenta

    Returns:
        (range): triangle-allowed range
    """

    ## return range(abs(j1-j2),j1+j2+1,1)  # range only works on integers
    values = [abs(j1-j2) + dj for dj in range(int((j1+j2)-abs(j1-j2))+1)]
    return values

if (__name__=="__main__"):

    # test is_integer
    for j in [0,0.5,1,1.5]:
        print("j {} is_integer {}".format(j,is_integer(j)))

    # test parity_sign
    for j in range(0,5):
        print("j {} parity_sign {}".format(j,parity_sign(j)))
    if (False):
        for j in [0.5]:
            print("j {} parity_sign {}".format(j,parity_sign(j)))
    
    # test parity_grade
    for P in [+1,-1]:
        print("P {:+2d} parity_grade {}".format(P,parity_grade(P)))
    if (False):
        for P in [0]:
            print("P {:+2d} parity_grade {}".format(P,parity_grade(P)))
