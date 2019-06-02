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
    04/05/19 (pjf):
        + Add one_body_transition_properties attribute.
        + Give get_rme() access to one-body observables.
    05/29/19 (mac): Add update method.
    06/02/19 (mac): Add deduced isoscalar and isovector E2 observable
      support.
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
    """Container for results from single MFDn run.

    Recall that there are inherited attributes:

        params (dict)
        energies (dict)
        num_eigenvalues (dict)
        filename (str)

    Data attributes:

        decompositions (dict): wave function probability decompositions

            Mapping: decomposition_name -> qn -> values

                decomposition name:

                  "Nex": decomposition in excitation quanta

                qn: (J,g,n)

                values (np.array vector): probabilities
 

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
    
                        "TBMEfile(2) = tbme-Trel.bin"
    
                    yields observable name "Trel"

        native_transition_properties (dict): MFDn-native transitions

            observable_name -> (qnf,qni) -> value

        one_body_static_properties (dict): ...

        one_body_transition_properties (dict): ...

    Accessors:
       [See definitions below.]

    """

    ########################################
    # Initializer
    ########################################

    def __init__ (self):
        """ Null-initialize standard attributes.
        """
        super().__init__()
        self.decompositions = {}
        self.native_static_properties = {}
        self.two_body_static_observables = {}
        self.native_transition_properties = {}
        self.one_body_static_properties = {}
        self.one_body_transition_properties = {}

    ########################################
    # Accessors
    ########################################

    def get_decomposition(self,decomposition_type,qn,verbose=False):
        """ Retrieve decomposition ("Nex") as np.array.
        """

        # validate decomposition type argument
        if (decomposition_type!="Nex"):
            raise(ValueError("invalid decomposition type {}".format(decomposition_type)))

        # retrieve decomposition
        try:
            decomposition = self.decompositions[decomposition_type][qn]
        except:
            return None

        return decomposition

    def get_radius(self,radius_type,qn,default=np.nan):
        """Retrieve rms radius value.

        We use the value from two-body "Relative radii" calculation
        (not the one-body radius reported by MFDn in oscillator runs,
        if it still even does that in v15).

        Arguments:
           radius_type (str): radius type rp/rn/r
           qn (tuple): quantum numbers for state
           default (float,optional): default value to return for missing radius

        Returns:
           (float): rms radius

        """

        # extract labels
        (J,gex,n) = qn

        rms_radius = self.two_body_static_observables.get(radius_type,{}).get(qn,default)

        return rms_radius

    def get_moment(self,observable,qn,default=np.nan,verbose=False):
        """Retrieve moment value.

        This accessor actually could retrieve *any* native static property, but
        we keep a name which reflects semantics rather than the internal
        represenation.

        Arguments:
           observable (str): moment type ("E2p","E2n","M1mu","M1lp","M1ln","M1sp",
               "M1sn"), as well as deduced cases ("E20", "E22")
           qn (tuple): quantum numbers for state
           default (float,optional): default value to return for missing radius

        Returns
           (float): observable value

        TODO(pjf): extend to give access to one_body_static_properties

        """

        # trap deduced isoscalar/isovector observables
        if (observable in {"E20","E21"}):
            E2p = self.get_moment("E2p",qn,default,verbose)
            E2n = self.get_moment("E2n",qn,default,verbose)
            if (observable =="E20"):
                value = 1/2*(E2p+E2n)
            else:
                value = 1/2*(E2p-E2n)
            return value

        # extract labels
        (J,gex,n) = qn

        value = self.native_static_properties.get(observable,{}).get(qn,default)

        return value

    def get_am(self,am_type,qn,default=np.nan):
        """Extract effective angular momentum value.

        Can be converted to effective angular momentum value with
        mfdnres.analysis.effective_am.

        Arguments:
           am_type (str): am type ("L","Sp","Sn","S")
           qn (tuple): quantum numbers for state
           default (float,optional): default value to return for missing radius

        Returns:
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

        Adds support for "E20" and E21" as isoscalar E2 and isovector E2.

        Arguments:
           observable (str): operator type ("E2p","E2n",...), as well as
               deduced cases ("E20", "E22")
           qn_pair (tuple): quantum numbers for states (qn_bra,qn_ket)
           default (float,optional): default value to return for missing radius

        Returns
           (float): observable value

        """

        # trap deduced isoscalar/isovector observables
        if (observable in {"E20","E21"}):
            E2p = self.get_rme("E2p",qn_pair,default,verbose)
            E2n = self.get_rme("E2n",qn_pair,default,verbose)
            if (observable =="E20"):
                value = 1/2*(E2p+E2n)
            else:
                value = 1/2*(E2p-E2n)
            return value

        # canonicalize labels
        (qn_pair_canonical,flipped,canonicalization_factor) = tools.canonicalize_Jgn_pair(
            qn_pair,tools.RMEConvention.kEdmonds
        )
        if (verbose):
            print("get_rme: observable {}, qn_pair {} => qn_pair_canonical {}".format(observable,qn_pair,qn_pair_canonical))

        # extract labels
        (qn_bra,qn_ket) = qn_pair_canonical
        (J_bra,g_bra,n_bra) = qn_bra
        (J_ket,g_ket,n_ket) = qn_ket

        # retrieve underlying rme
        if self.params.get("hw", 0) != 0:
            try:
                rme = (canonicalization_factor
                    * self.native_transition_properties[observable][qn_pair_canonical])
            except KeyError:
                try:
                    rme = (canonicalization_factor
                        * self.one_body_transition_properties[observable][qn_pair_canonical])
                except KeyError:
                    return default
        else:
            try:
                rme = (canonicalization_factor
                    * self.one_body_transition_properties[observable][qn_pair_canonical])
            except KeyError:
                return default

        if (verbose):
            print("    rme {:e}".format(rme))


        return rme

    def get_rtp(self,observable,qn_pair,default=np.nan,verbose=False):
        """ Retrieve reduced transition probability (RTP).
        """

        # extract labels
        (qn_bra,qn_ket) = qn_pair
        (J_bra,g_bra,n_bra) = qn_bra
        (J_ket,g_ket,n_ket) = qn_ket

        # retrieve underlying rme
        try:
            rme = self.get_rme(observable,qn_pair,default,verbose)
        except KeyError:
            return default

        # derive final value from rme
        rtp = 1/(2*J_ket+1)*rme**2

        return rtp

    ########################################
    # Manipulators
    ########################################
        
    def update(self,other):
        """Merge in data from other MFDnResultsData object.

        Behavior for inherited attributes (e.g., energies) is as defined by ResultsData.update().

        Arguments:
            other (MFDnResultsData): other results set to merge in

        """
        super().update(other)

        # merge observable dictionaries


        def update_observable_dictionary(self_dict,other_dict):
            """ Merge two observable dictionaries of type name->qn->value.

            As currently implemented, data for a named observable is only merged if
            this named observable arises in "self", not if it is only found in
            "other".

            Arguments:
                self_dict (dict): observable dictionary to be updated
                other_dict (dict): observable dictionary providing update data
            """
            
            # add dictionaries for any missing named observables
            for name in other_dict.keys() - self_dict.keys():
                self_dict[name] = {}
            
            # update dictionaries for all named
            for name in self_dict.keys():
                self_dict[name].update(other_dict[name])

        update_observable_dictionary(self.decompositions,other.decompositions)
        update_observable_dictionary(self.native_static_properties,other.native_static_properties)
        update_observable_dictionary(self.two_body_static_observables,other.two_body_static_observables)
        update_observable_dictionary(self.native_transition_properties,other.native_transition_properties)
        update_observable_dictionary(self.one_body_static_properties,other.one_body_static_properties)
        update_observable_dictionary(self.one_body_transition_properties,other.one_body_transition_properties)

#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
