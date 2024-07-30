""" Arithmetic with uncertainties.

    Numbers with uncertainties are treated as tuples (x,dx).  A more natural
    approach will be to define them as object with overloaded arithmetic
    operators.

    Mark A. Caprio
    University of Notre Dame

    - 10/07/21 (mac): Created, inspired by old HP48S RPL routines SIGMAMATH.
    - 07/30/24 (mac): Force uncertainties to be positive.

"""

import numpy as np

def add(x_sigma,y_sigma):
    """ Add numbers with uncertainties.

    Arguments:

       x_sigma (tuple): (x,dx)

       y_sigma (tuple): (y,dy)

    Returns:

        (tuple): (z,dz)
    """

    x, dx = x_sigma
    y, dy = y_sigma

    z = x+y

    dz = np.sqrt(dx**2+dy**2)

    return (z,dz)

def mul(x_sigma,y_sigma):
    """ Multiply numbers with uncertainties.

    Arguments:

       x_sigma (tuple): (x,dx)

       y_sigma (tuple): (y,dy)

    Returns:

        (tuple): (z,dz)
    """

    x, dx = x_sigma
    y, dy = y_sigma

    z = x*y

    dz = abs(z*np.sqrt((dx/x)**2+(dy/y)**2))

    return (z,dz)

def pow(x_sigma,p):
    """ Exponentiate number with uncertainty.

    Arguments:

       x_sigma (tuple): (x,dx)

       p (float): power

    Returns:

        (tuple): (z,dz)
    """

    x, dx = x_sigma

    z = x**p

    dz = abs(z*p*(dx/x))

    return (z,dz)

def scale(x_sigma,c):
    """ Scale number (i.e., multiply by scalar) with uncertainty.

    Derived operation.

    Arguments:

       x_sigma (tuple): (x,dx)

       c (float): scale

    Returns:

        (tuple): (z,dz)
    """

    return mul(x_sigma,(c,0))

def sub(x_sigma,y_sigma):
    """ Subtract numbers with uncertainties.

    Derived operation.

    Arguments:

       x_sigma (tuple): (x,dx)

       y_sigma (tuple): (y,dy)

    Returns:

        (tuple): (z,dz)
    """

    return add(x_sigma,scale(y_sigma,-1))

def div(x_sigma,y_sigma):
    """ Divide numbers with uncertainties.

    Derived operation.

    Arguments:

       x_sigma (tuple): (x,dx)

       y_sigma (tuple): (y,dy)

    Returns:

        (tuple): (z,dz)
    """

    return mul(x_sigma,pow(y_sigma,-1))
    


################################################################
# main
################################################################
  
def main():

    print(add((1.,0.5),(1,0.5)))
    print(scale((1,0.5),-1))
    print(pow((100.,0.5),1/2))
    print(sub((1.,0.5),(1,0.5)))

if __name__ == "__main__":
    main()
