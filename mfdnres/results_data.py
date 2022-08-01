"""results_data.py

    Define base (interface) class for results data storage.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    10/10/17 (mac): Extracted from res.py.
    09/06/18 (pjf): Replace get_levels() with levels property.
    05/29/19 (mac): Add update method.
    01/04/20 (mac): Convert num_eigenvalues from static data to property.
"""

from __future__ import annotations

import collections
import typing
import numpy as np

# import for type hinting
from .tools import (
    SubspaceType,
    LevelQNType,
    LevelQNPairType,
    OperatorQNType,
)

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

        filename (str): Source results filename for this mesh point.

    Properties:

        levels (dict): list of quantum number tuples (J,g,n)

        num_eigenvalues (dict):  Number of eigenvalues, as dictionary subspace->number.

            Here subspace is a tuple of quantum numbers, typically (J,g) or (J,gex).

    Accessors:

        [See definitions.]

    Manipulators:

        [See definitions.]


    """

    params:dict
    energies:dict[LevelQNType,float]
    filename:str

    ########################################
    # Initializer
    ########################################

    def __init__(self):
        """ Null-initialize standard attributes.
        """
        self.params = {}
        self.energies = {}
        self.filename = ""

    ########################################
    # Accessors
    ########################################

    @property
    def levels(self):
        """List of quantum number tuples (J,g,n) for levels, sorted by increasing energy eigenvalue.

        Note some parsers (legacy) may store (J,gex,n), in which case that is
        the interpretation of the tuple.

        Returns:
            (list of tuples): list of quantum numbers

        """
        # Make a list of unsorted quantum number tuples
        raw_qn_list = list(self.energies.keys())
        # Sort the quantum numbers based on their associated energy
        qn_list = sorted(raw_qn_list, key=(lambda qn: self.energies[qn]))
        return qn_list

    @property
    def num_eigenvalues(self):
        """Dictionary of subspace dimensions (J,g)->dim for levels.

        Note some parsers (legacy) may store (J,gex), in which case that is
        the interpretation of the tuple.

        Returns:
            (dict): subspace dimension dictionary

        """
        # tally up levels by subspace
        num_eigenvalues_by_subspace:dict[SubspaceType,int]={}
        for qn in self.energies.keys():
            (J,g,_)=qn
            num_eigenvalues_by_subspace[(J,g)]=num_eigenvalues_by_subspace.setdefault((J,g),0)+1
        return num_eigenvalues_by_subspace

    def get_levels(self):
        """DEPRECATED -- Retrieve list of quantum number tuples (J,g,n) for levels, sorted by increasing energy eigenvalue.

        Note some parsers (legacy) may store (J,gex,n), in which case that is the interpretation of the tuple.

        Returns:
            (list of tuples): list of quantum numbers

        """
        raise DeprecationWarning("accessor get_levels() is deprecated; use property instead")
        return self.levels

    def get_energy(self,qn:LevelQNType,default=np.nan):
        """Retrieve the energy of the level with given quantum numbers.

        Returns a default "flag" value if the quantum numbers are not
        found among the calculate levels.

        Arguments:
            qn (tuple): tuple of quantum numbers (J,g,n)
            default (numeric,optional): default flag value to use

        Returns:
            value (float): energy eigenvalue, or default if qn not found

        """
        return self.energies.get(qn, default)

    ########################################
    # Manipulators
    ########################################

    def update(self,other:ResultsData):
        """Merge in data from other ResultsData object.

        Merges energy listing and updates eigenvalue count.  Metadata (params
        and filename) are retained from self.

        Arguments:
            other (ResultsData): other results set to merge in

        """

        self.energies.update(other.energies)

#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
