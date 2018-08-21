""" mfdn_results_data.py

    Result storage and access for mfdn runs.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    10/6/17 (mac): Extract MFDnResultsData from res.py.
    10/23/17 (mac): Add get_radius accessor.
"""

import math

import numpy as np

from . import (
    am,
    results_data,
    )


#################################################
# MFDnResultsData
#################################################

class MFDnResultsData(results_data.ResultsData):
    """Container for MFDn results.

    TODO clean up docstring

    TODO: finish neatening transition phases
        - check GT conjugation phase
        - define setter method for transitions, to only store in
           canonically decreasing (f<i) order
        - reduce states data member to sorted list

    Inherited attributes:
        params (dict)
        energies (dict)
        num_eigenvalues (dict)
        filename (str)

    Attributes:
        decompositions (dict): wave function amplitude decompositions

            Mapping: decomposition_name -> qn -> np.array (vector)

            decomposition name:

              "Nex": decomposition in excitation quanta

        native_static_properties (dict): native-calculated static properties

            Mapping: property_name -> qn -> value

            property_name:

                "T": native-calculated isospin (only valid in pn-symmetric basis)

                "M1<type>": M1 moments (<type> in [mu,lp,ln,sp,sn])

                "E2<type>": E2 moments (<type> in [p,n])

                "am<type>": native-calculated angular momentum (only
                valid in pn-symmetric basis???) (<type> in [L,S,Sp,Sn,J])

        two_body_static_observables (dict): two-body static observables

            Mapping: observable_name -> qn -> value

            observable_name:

                "r<type>": radii (<type> in [p,n,""])

                others taken from TBME filenames, e.g.,

                    "TBMEfile(2) = tbme-Trel.bin" -> observable name "Trel"


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

            Initializes self.params, self.energies from BaseResultsData.  Also initializes
            self.states, self.properties, self.transitions, self.tbo, self.moments,
            and self.orbital_occupations.
        """
        super().__init__()
        self.decompositions = {}
        self.native_static_properties = {}
        self.two_body_static_observables = {}


    ########################################
    # Accessors
    ########################################

    def get_radius(self,radius_type,qn,default=np.nan):
        """Retrieve rms radius value.

        We use the value from two-body "Relative radii" calculation
        (not the one-body radius reported by MFDn in oscillator runs,
        if it still even does that in v15).

        Arguments:
           radius_type (str): radius type rp/rn/r
           qn (tuple): quantum numbers for state
           default (float,optional): default value to return for missing radius

        Returns
           (float): rms radius

        """

        # extract labels
        (J,gex,n) = qn

        rms_radius = self.two_body_static_observables.get(radius_type, {}).get(qn, default)

        return rms_radius

#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
