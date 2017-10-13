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

            Ex: For spncci (runmac0420)...

            {'hw': 10.0, 'Nsigmamax': 2, 'nuclide': (3, 3),
            'nuclide-N': 3, 'A': 6, 'use_coulomb': False, 'nuclide-Z':
            3, 'Nmax': 2, 'observable_names': ['hamiltonian',
            'r2intr', 'Qintr'], 'interaction': 'JISP16', 'nuclide.N':
            3, 'nuclide.Z': 3, 'N1v': 1}

        energies (dict): Energy eigenvalues, as dictionary qn->energy.

            Here qn is a tuple of quantum numbers, typically (J,g,i) or (J,gex,i).

        num_eigenvalues (dict):  Number of eigenvalues, as dictionary subspace->number.

            Here subspace is a tuple of quantum numbers, typically (J,g) or (J,gex).

        filename (str): Source results filename for this mesh point.

    Accessors:

        get_levels: List of qn tuples for all levels.  Values are sorted by increasing energy eigenvalue.
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
        """ Retrieve list of quantum number tuples (J,g,n) for levels, sorted by increasing energy eigenvalue.

        Note some parsers (legacy) may store (J,gex,n), in which case that is the interpretation of the tuple.

        Returns:
            (list of tuples): list of quantum numbers

        """
        # Makes a list of unsorted quantum number tuples
        raw_qn_list = list(self.energies.keys())
        # Sorts the quantum numbers based on their associated energy
        qn_list = sorted(raw_qn_list,key=(lambda qn : self.energies[qn]))
        return qn_list

    def get_energy(self,qn,default=np.nan):
        """Retrieve the energy of the level with given quantum numbers.

        Returns a default "flag" value if the quantum numbers are not
        found among the calculate levels.

        Arguments:
            qn (tuple): tuple of quantum numbers (J,g,n)
            default (numeric,optional): default flag value to use

        Returns:
            value (float): energy eigenvalue, or default if qn not found

        """

        try:
            ##print("get_energy",self.energies,qn,qn in self.energies)
            value = self.energies[qn]
        except:
            value = default
        return value

#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
