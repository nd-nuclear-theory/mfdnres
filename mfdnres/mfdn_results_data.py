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
    06/24/19 (mac): Make update_observable_dictionary robust against
        missing observables.
    06/25/19 (mac): Add get_isospin accessor.
    07/26/19 (mac): Support "M1" as deduced observable in get_moment
        and get_rme.
    09/24/19 (mac): Move update_observable_dictionary out to results_data.py.
    11/17/19 (mac): Extend get_moment to use self-transition rme first if available.
    01/04/20 (mac): Add "M1-native" special case to get_moment.
    04/27/20 (mac): Correct normalization of deduced isocalar/isovector quadrupole RMEs
        to match definitions in intrinsic.

"""

import math

import numpy as np

from . import (
    am,
    results_data,
    tools,
    )

#################################################
# helper function for merging observable dictionaries
#################################################

def update_observable_dictionary(self_dict,other_dict):
    """Merge two observable dictionaries of type name->qn->value.

    Allows for possibility that an observable name may exist only in
    self_dict or only in other_dict.

    Arguments:
        self_dict (dict): observable dictionary to be updated
        other_dict (dict): observable dictionary providing update data

    """
    
    observables = set(self_dict.keys()).union(other_dict.keys())
    for name in observables:
        # ensure observable exists in self_dict
        if (name not in self_dict.keys()):
            self_dict[name] = {}
        # update observable from other_dict
        if name in other_dict.keys():
            self_dict[name].update(other_dict[name])

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

    ================================================================

    05/03/20 (mac): TODO proposed cleanup of dictionary names
          
    Here is a proposed renaming based on:

       <source>_<particle-rank>_<observable-or-property-type>

    with the exceptional case of

       <source=mfdn>_<particle-rank=level>_properties

    for extracted level quantum numbers.

    - - - -

    native_static_properties
      => mfdn_level_properties (T, am)
         AND
         mfdn_one_body_moments (E2, M1)

    native_transition_properties
      => mfdn_one_body_transitions

    two_body_static_observables
      => mfdn_two_body_expectation_values

    one_body_static_properties
      => UNUSED, since postprocessors will always output RMEs labeled by initial and final state

    one_body_transition_properties
      => postprocessor_one_body_rmes

    PROPOSED
      => postprocessor_two_body_rmes

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

    def get_isospin(self,qn,default=np.nan,verbose=False):
        """ Retrieve effective isospin T.

        Arguments:
            qn (tuple): quantum numbers for state
            default (float,optional): default value to return for missing observable

        Returns:
            (float): observable value
        """

        # retrieve decomposition
        try:
            value = self.native_static_properties["T"][qn]
        except:
            return default

        return value

    def get_decomposition(self,decomposition_type,qn,verbose=False):
        """ Retrieve decomposition ("Nex") as np.array.

        Arguments:
            decomposition_type (str): decomposition type ("Nex")
            qn (tuple): quantum numbers for state

        Returns:
            (np.array): decomposition probabilities
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
            default (float,optional): default value to return for missing observable

        Returns:
            (float): rms radius

        """

        # extract labels
        (J,gex,n) = qn

        rms_radius = self.two_body_static_observables.get(radius_type,{}).get(qn,default)

        return rms_radius

    def get_moment(self,observable,qn,default=np.nan,verbose=False):
        """Retrieve moment value.

        The preferred method of deducing the moment is from the RME from the
        state to itself, if available via get_rme.  (This RME may be either a
        native # calculated rme or from the postprocessor.)  If the RME is not
        available, the accessor looks for a native calculated moment.  But
        beware that mfdnv15b02 outputs junk native moments if OBDME calculation
        has been turned out, and native calculated moments (as well as native
        calculated transitions) are output with fixed low precision.  Hence our
        preference for using the RME to calculate the moment.

        Adds support for deduced cases: "E20" and E21" for isoscalar E2 and
        isovector E2, "M1" deduced from dipole terms.  (Any native-calculated
        "M1" from MFDn is ignored, as this is redundant to the dipole terms and
        is not provided by obscalc-ob.)

        Arguments:
            observable (str): moment type ("E2p", "E2n", "Dlp", "Dln", "Dsp",
                "Dsn", ...), as well as deduced cases ("E20", "E21", "M1",
                "M1-native")
            qn (tuple): quantum numbers for state
            default (float,optional): default value to return for missing observable

        Returns
            (float): observable value

        From berotor2 footnote 3:

            mu(J)=math.sqrt(4*math.pi/3)/am.hat(J)*(JJ10|JJ)*<J||M1||J>

            (JJ10|JJ) = math.sqrt((J)/((J + 1)))

        From berotor2 footnote 1:

            eQ(J)=math.sqrt(16*math.pi/5)/am.hat(J)*(JJ20|JJ)*<J||Q2||J>
            (JJ20|JJ) = math.sqrt((J*(2*J - 1))/((J + 1)*(2*J + 3)))

        """

        # trap deduced observables (isoscalar/isovector E2 or physical M1)
        if (observable in {"E20","E21"}):
            E2p = self.get_moment("E2p",qn,default,verbose)
            E2n = self.get_moment("E2n",qn,default,verbose)
            if (observable =="E20"):
                value = E2p+E2n
            else:
                value = E2p-E2n
            return value
        elif (observable == "M1"):
            Dsp = self.get_moment("Dsp",qn,default,verbose)
            Dsn = self.get_moment("Dsn",qn,default,verbose)
            Dlp = self.get_moment("Dlp",qn,default,verbose)
            gp = 5.586
            gn = -3.826
            value = Dlp+gp*Dsp+gn*Dsn
            return value
        elif (observable == "M1-native"):
            # force retrieval of MFDn.res precomputed physical value
            #
            # If no RMEs are available, and only the MFDn.res static M1 moment
            # observables are available, this may yield higher precision than
            # the M1 value this function would usually compute as a linear
            # combination of the D terms taken from MFDn.res.  These terms
            # typically have smaller magnitudes but are saved at the same fixed
            # point precision at the physical M1.
            value = self.native_static_properties.get("M1",{}).get(qn,default)
            return value

        # extract labels
        (J,gex,n) = qn

        # get moment from rme, if available
        if (observable in ["Dlp", "Dln","Dsp", "Dsn"]):
            rme_prefactor = math.sqrt(4*math.pi/3)/am.hat(J)*math.sqrt((J)/((J + 1)))
        elif (observable in ["E2p", "E2n"]):
            rme_prefactor = math.sqrt(16*math.pi/5)/am.hat(J)*math.sqrt((J*(2*J - 1))/((J + 1)*(2*J + 3)))
        else:
            raise(ValueError("unrecognized moment type {}".format(observable)))
        value = rme_prefactor*self.get_rme(observable,(qn,qn))
        
        # else revert to native static moment (may be garbage if obmes were turned off)
        if (np.isnan(value)):
            value = self.native_static_properties.get(observable,{}).get(qn,default)

        return value

    def get_am(self,am_type,qn,default=np.nan):
        """Extract effective angular momentum value.

        Can be converted to effective angular momentum value with
        mfdnres.analysis.effective_am.

        Arguments:
            am_type (str): am type ("L","Sp","Sn","S")
            qn (tuple): quantum numbers for state
            default (float,optional): default value to return for missing observable

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

        If hw parameter is nonzero, first recourse is to native calculated rme,
        then fallback is postprocessor rme.  But for hw=0, assumption is that
        non-oscillator basis was used, and native rmes are not meaningful, so
        accessor only looks for postprocessor rme.

        Adds support for deduced cases: "E20" and E21" for isoscalar E2 and
        isovector E2, "M1" deduced from dipole terms.  (Any native-calculated
        "M1" from MFDn is ignored, as this is redundant to the dipole terms and
        not provided by obscalc-ob.)

        Arguments:
            observable (str): operator type ("E2p", "E2n", "Dlp", "Dln",
                "Dsp", "Dsn", ...), as well as deduced cases ("E20", "E21","M1")
            qn_pair (tuple): quantum numbers for states (qn_bra,qn_ket)
            default (float,optional): default value to return for missing observable

        Returns
            (float): observable value

        """

        # trap deduced observables (isoscalar/isovector E2 or physical M1)
        if (observable in {"E20","E21"}):
            E2p = self.get_rme("E2p",qn_pair,default,verbose)
            E2n = self.get_rme("E2n",qn_pair,default,verbose)
            if (observable =="E20"):
                value = E2p+E2n
            else:
                value = E2p-E2n
            return value
        elif (observable == "M1"):
            Dsp = self.get_rme("Dsp",qn_pair,default,verbose)
            Dsn = self.get_rme("Dsn",qn_pair,default,verbose)
            Dlp = self.get_rme("Dlp",qn_pair,default,verbose)
            gp = 5.586
            gn = -3.826
            value = Dlp+gp*Dsp+gn*Dsn
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

    def get_rme_matrix(self,observable,subspace_pair,dim_pair,default=np.nan,verbose=False):
        """Construct matrix of reduced matrix elements (RMEs).

        Resulting matrix may be useful either for diagnostic purposes (to review
        available transitions) or in block calculations with RMEs.

        See get_rme for conventions regarding RME itself.
    
        Example:
            >>> mesh_point.get_rme_matrix("E2p",((2,0),(2,0)),(4,4),verbose=True)


        Arguments:
            observable (str): operator type ("E2p", ...) accepted by get_rme
            subspace_pair (tuple): quantum numbers for subspaces (Jg_bra,Jg_ket)
            dim_pair (tuple): dimensions (dim_bra,dim_ket)
            default (float,optional): default value to return for missing observable

        Returns
            (np.array of float): observable values

        """

        # unpack arguments
        (Jg_bra,Jg_ket) = subspace_pair
        (J_bra,g_bra) = Jg_bra
        (J_ket,g_ket) = Jg_ket
        (bra_dim,ket_dim) = dim_pair

        # assemble matrix
        rme_matrix = np.empty((bra_dim,ket_dim))
        for bra_index in range(bra_dim):
            for ket_index in range(ket_dim):
                qn_bra = (J_bra,g_bra,bra_index+1)  # convert to spectroscopic 1-based numbering
                qn_ket = (J_ket,g_ket,ket_index+1)  # convert to spectroscopic 1-based numbering
                rme_matrix[bra_index,ket_index] = self.get_rme(
                    observable,(qn_bra,qn_ket),
                    default=default,verbose=verbose
                )

        if (verbose):
            print("    rme matrix")
            print(rme_matrix)
            
        return rme_matrix
                
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
    # Updating method
    ########################################
        
    def update(self,other):
        """Merge in data from other MFDnResultsData object.

        Behavior for inherited attributes (e.g., energies) is as defined by ResultsData.update().

        Arguments:
            other (MFDnResultsData): other results set to merge in

        """
        super().update(other)

        # merge observable dictionaries
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
