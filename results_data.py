"""results_data.py

    Define base (interface) class for results data storage.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    10/10/17 (mac): Extracted from res.py.
"""

import numpy as np

#################################################
# ResultsData
#################################################

class ResultsData (object):
    """Object to hold results for a single run "mesh point".

    This is an interface class.  That is, ResultsData should never be
    instantiated.  Children should be defined to hold results from
    each specific code (such as mfdn or spncci), and only these
    children should be instantiated.

    Attributes:
        params (dict): Mesh point run parameters, as dictionary of keyword->value.

            Ex: For spncci, 

        energies (dict): Energy eigenvalues, as dictionary qn->energy.

            Here qn is a tuple of quantum numbers, typically (J,g,i) or (J,gex,i).

        num_eigenvalues (dict):  Number of eigenvalues, as dictionary subspace->number.

            Here subspace is a tuple of quantum numbers, typically (J,g) or (J,gex).

        filename (str): Source results filename for this mesh point.

    Accessors:
            get_levels:  Accessor for all quantum numbers.
                Takes no arguments are returns a list of all quantum numbers produced
                by the run, sorted based on the energy associated with each set of quantum numbers.
            get_energy:  Accessor for energy by quantum number tuple.
                Takes as an argument a tuple of quantum numbers.  The the set of quantum numbers
                is valid, it returns the energy associated with those quantum numbers.  If the quantum
                numbers are not valid, it returns None and prints a message to the console.        
        Methods:

    """
    ########################################
    # Initializer
    ########################################
    def __init__ (self):
        """
            Arguments:
                None.
            Returned:
                None.

            Initializes self.params, self.energies, self.states, and self.properties.  
        """
        self.params = {}
        self.energies = {}
        self.num_eigenvalues = {}
        self.filename = ""

    ########################################
    # Accessors
    ########################################        
    def get_levels(self):
        """
            Arguments:
                None.
            Returned:
                qn_list (list): A list of quantum numbers sorted by their
                    associated energy

            Returns a list of quantum numbers ((J, g, n) for MFDn or (J, gex, i)
            for SpNCCI), sorted by the ground state energy associated with set of
             quantum numbers.
        """
        # Makes a list of unsorted quantum number tuples
        raw_qn_list = list(self.energies.keys())
        # Sorts the quantum numbers based on their associated energy
        qn_list = sorted(raw_qn_list,key=(lambda qn : self.energies[qn]))
        return qn_list

    def get_energy(self,qn,default=np.nan):
        """ Retrieve the energy of level with given quantum numbers.

            Arguments:
                qn (tuple): A tuple of quantum numbers.
            Returned:
                value(varies): The energy for a valid set of quantum numbers.
                    If the set of quantum numbers is valid, value is set to the
                    ground state energy associated with the quantum numbers.  If the set of
                    quantum numbers is not valid, value is set to None. 

            Returns the ground state associated with the quantum numbers associated with 'qn', 
            if they are valid.  If they are not valid, it returns the None, and prints
            a message to the console.
        """
        # Check to be sure the quantum numbers supplied are in self.energies
        try:
            value = self.energies[qn]
        except:
            return default
        return value

#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
