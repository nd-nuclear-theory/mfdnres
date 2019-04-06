""" mfdn_results_data.py

    Result storage and access for mfdn runs.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    10/06/17 (mac): Extract MFDnResultsData from res.py.
    10/23/17 (mac): Add get_radius accessor.
    09/06/18 (pjf):
        + Add native_transition_properties attribute.
        + Implement get_rme() and get_rtp().
    12/14/18 (mac): Add get_moment accessor.
    04/02/19 (mac): Add get_am accessor.
"""

import math

import numpy as np

from . import (
    am,
    results_data,
    tools,
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

        native_transition_properties (dict): MFDn-native transitions

            observable_name -> (qnf,qni) -> value


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
        self.native_transition_properties = {}


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

        rms_radius = self.two_body_static_observables.get(radius_type,{}).get(qn,default)

        return rms_radius

    def get_moment(self,moment_type,qn,default=np.nan):
        """Retrieve moment value.

        This accessor actually could retrieve *any* native static property, but
        we keep a name which reflects semantics rather than the internal
        represenation.

        Arguments:
           moment_type (str): moment type ("E2p","E2n","M1mu","M1lp","M1ln","M1sp","M1sn")
           qn (tuple): quantum numbers for state
           default (float,optional): default value to return for missing radius

        Returns
           (float): observable value

        """

        # extract labels
        (J,gex,n) = qn

        value = self.native_static_properties.get(moment_type,{}).get(qn,default)

        return value

    def get_am(self,am_type,qn,default=np.nan):
        """Extract effective angular momentum value.

        Can be converted to effective angular momentum value with
        mfdnres.analysis.effective_am.

        Arguments:
           am_type (str): am type ("L","Sp","Sn","S")
           qn (tuple): quantum numbers for state
           default (float,optional): default value to return for missing radius

        Returns
           (float): observable value

        """

        # extract labels
        (J,gex,n) = qn

        # validate am type name (to protect user from nonsensical Lp and Ln)
        if (am_type not in {"L","Sp","Sn","S"}):
            raise ValueError("Invalid angular momentum type name {}".format(am_type))

        # extract effective am value
        value_sqr = self.native_static_properties.get(am_type+"_sqr",{}).get(qn,default)
        if (value_sqr<0):
            raise ValueError("Invalid squared angular momentum value: {:e}".format(value_sqr))
        value = tools.effective_am(value_sqr)

        return value


    def get_rme(self,observable,qn_pair,default=np.nan,verbose=False):
        """Retrieve reduced matrix element (RME).

        Returns RME in Edmonds convention, as common for spectroscopic data
        analysis in the shell model community.

        """

        # canonicalize labels
        (qn_pair_canonical,flipped,canonicalization_factor) = tools.canonicalize_Jgn_pair(
            qn_pair,tools.RMEConvention.kEdmonds
        )

        # extract labels
        (qn_bra,qn_ket) = qn_pair_canonical
        (J_bra,g_bra,n_bra) = qn_bra
        (J_ket,g_ket,n_ket) = qn_ket

        # retrieve underlying rme
        try:
            rme = (canonicalization_factor
                * self.native_transition_properties[observable][qn_pair_canonical])
        except KeyError:
            return default

        return rme

    def get_rtp(self,observable,qn_pair,default=np.nan):
        """ Retrieve reduced transition probability (RTP).
        """

        # extract labels
        (qn_bra,qn_ket) = qn_pair
        (J_bra,g_bra,n_bra) = qn_bra
        (J_ket,g_ket,n_ket) = qn_ket

        # retrieve underlying rme
        try:
            rme = self.get_rme(observable,qn_pair)
        except KeyError:
            return default

        # derive final value from rme
        rtp = 1/(2*J_ket+1)*rme**2

        return rtp

#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
