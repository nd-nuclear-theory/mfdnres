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
    06/02/19 (mac): Add deduced isoscalar and isovector E2 observable support.
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
    04/27/20 (zz): Add deduced isoscalar and isovector M1 observable support.
    09/17/20 (mac): Overhaul observable data attribute naming scheme.
    05/25/21 (mac): Add "E2" as alias for "E2p", for consistency in plot observable labeling.
    09/25/21 (mac): Add single-species radii relative to own center of mass in get_radius().
    10/20/21 (mac): Trap spurious one-body E0 in get_rme()
    02/24/22 (zz): Add E1 special case in get_rme().
    07/12/22 (mac): Provide support for two-body observables in get_rme().
    07/13/22 (mac): Add get_expectation_value().
    07/17/22 (mac): Deduce diagonal E0 RMEs from radius (or else suppress) in get_rme().
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

OBSERVABLE_ALIASES = {
    # physical electric operators
    "E0": "E0p",
    "E1": "E1p",
    "E2": "E2p",
    # sensible names for M1 components (a.k.a. "dipole terms")
    "M1lp": "Dlp",
    "M1ln": "Dln",
    "M1sp": "Dsp",
    "M1sn": "Dsn",
}

class MFDnResultsData(results_data.ResultsData):
    """Container for results from single MFDn run.

    Recall that there are inherited attributes:

        params (dict)
        energies (dict)
        num_eigenvalues (dict)
        filename (str)

    Naming scheme for obserable attributes:

        <source=mfdn|postprocessor>_<particle-rank=level|ob|tb>_<observable-or-property-type>

            particle-rank:
               level = level quantum numbers (no associated particle rank)
               ob = one-body
               tb = two-body

            Note: The same observable as an "expectation value" or an "RME" will
            differ by an angular momentum factor, as RMEs are stored under the
            Edmonds convention for the Wigner-Eckart theorem.


    Data attributes:

        mfdn_level_decompositions (dict): wave function probability decompositions

            Mapping: decomposition_name -> qn -> values

                decomposition name:

                  "Nex": decomposition in excitation quanta

                qn: (J,g,n)

                values (np.array vector): probabilities

        mfdn_level_residuals (dict): mfdn Lanczos residuals

            Mapping: qn -> value

        mfdn_level_properties (dict): mfdn extracted mean quantum numbers (J, T)

            Mapping: property_name -> qn -> value

                property_name:

                    "T": native-calculated isospin (only valid in pn-symmetric basis)

                    "am<type>": native-calculated angular momentum (only
                    valid in pn-symmetric basis???) (<type> in [L,S,Sp,Sn,J])

        mfdn_ob_moments (dict): mfdn one-body moments (only well-defined for oscillator basis)

            Mapping: observable_name -> qn -> value

                observable_name:

                    "M1", "Dlp", "Dln", "Dsp", "Dsn": M1 moments

                       PROPOSED ALTERNATIVE: "M1<type>": M1 moments (<type> in [mu,lp,ln,sp,sn])

                    "E2<type>": E2 moments (<type> in [p,n])

        mfdn_tb_expectations (dict): two-body static observables

            Mapping: observable_name -> qn -> value

                observable_name:

                    "r<type>": radii (<type> in [p,n,""])

                    others taken from TBME filenames, e.g.,

                        "TBMEfile(2) = tbme-Trel.bin"

                    yields observable name "Trel"

        mfdn_ob_rmes (dict): MFDn-native transitions

            Mapping: observable_name -> (qnf,qni) -> value

        postprocessor_ob_rmes (dict): One-body RMEs calculated by postprocessor and obscalc-ob.

            Mapping: observable_name -> (qnf,qni) -> value

        postprocessor_tb_rmes (dict): Two-body RMEs calculated by postprocessor.

            Mapping: observable_name -> (qnf,qni) -> value

    Accessors:
       [See definitions below.]
    """

    # Data attribute renaming 09/17/20 (mac)
    #
    # decompositions
    #     => mfdn_level_decompositions
    #
    # native_static_properties  TODO
    #     => mfdn_level_properties (T, am)
    #        AND
    #        mfdn_ob_moments (E2, M1)
    #
    # native_transition_properties
    #     => mfdn_ob_rmes
    #
    # two_body_static_observables
    #     => mfdn_tb_expectations
    #
    # one_body_static_properties
    #     => REMOVED
    #
    # one_body_transition_properties
    #     => postprocessor_ob_rmes
    #
    # NEW
    #     => postprocessor_tb_rmes

    ########################################
    # Initializer
    ########################################

    def __init__ (self):

        """ Null-initialize standard attributes.
        """
        super().__init__()
        self.mfdn_level_decompositions = {}
        self.mfdn_level_residuals = {}
        self.mfdn_level_properties = {}
        self.mfdn_ob_moments = {}
        self.mfdn_ob_rmes = {}
        self.mfdn_tb_expectations = {}
        self.postprocessor_ob_rmes = {}
        self.postprocessor_tb_rmes = {}

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
            value = self.mfdn_level_properties["T"][qn]
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
            decomposition = self.mfdn_level_decompositions[decomposition_type][qn]
        except:
            return None

        if verbose:
            print("{} {} {}".format(decomposition_type,qn,decomposition))

        return decomposition

    def get_radius(self,radius_type,qn,default=np.nan):
        """Retrieve rms radius value.

        We use the value from two-body "Relative radii" calculation
        (not the one-body radius reported by MFDn in oscillator runs,
        if it still even does that in v15).

        Arguments:

            radius_type (str): radius type ("rp", "rn", "r"), as well as deduced
                cases ("rp-ss", "rn-ss")

            qn (tuple): quantum numbers for state

            default (float,optional): default value to return for missing observable

        Returns:
            (float): rms radius

        """

        # extract labels
        (J,gex,n) = qn

        # trap deduced observables (single species radii relative to own center
        # of mass)
        #
        # Relation to MFDn output observables "r_pp" and "r_nn" deduced from
        # cshalo [PRC 90, 034305 (2014)] (A5).
        if (radius_type == "rp-ss"):
            nuclide = self.params["nuclide"]
            Np, Nn = nuclide
            A = sum(nuclide)
            rpp = self.mfdn_tb_expectations.get("rpp",{}).get(qn,default)
            value = A/Np*rpp
            return value
        elif (radius_type == "rn-ss"):
            nuclide = self.params["nuclide"]
            Np, Nn = nuclide
            A = sum(nuclide)
            rnn = self.mfdn_tb_expectations.get("rnn",{}).get(qn,default)
            value = A/Nn*rnn
            return value

        rms_radius = self.mfdn_tb_expectations.get(radius_type,{}).get(qn,default)

        return rms_radius

    def get_moment(
        self, observable, qn,
        allow_mfdn_native=True, allow_moment_from_rme=True, allow_e0_from_radius=True,
        default=np.nan, verbose=False
    ):
        """Retrieve moment value.

        The preferred method of deducing the moment is from the RME from the
        state to itself, if available via get_rme.  (This RME may be either a
        native calculated rme as provided by older versions of mfdn or from the
        postprocessor.)  If the RME is not available, the accessor looks for a
        native calculated moment.  But beware that mfdn15 outputs junk
        native moments if OBDME calculation has been turned off, and native
        calculated moments (as well as native calculated transitions) are output
        with fixed low precision.  Hence our preference for using the RME to
        calculate the moment.

        Adds support for deduced cases: "E20" and E21" for isoscalar E2 and
        isovector E2, "M1" deduced from dipole terms.  (Any native-calculated
        "M1" from MFDn is ignored, as this is redundant to the dipole terms and
        is not provided by obscalc-ob.)

        Arguments:

            observable (str): moment type ("E2p", "E2n", "Dlp", "Dln", "Dsp",
                "Dsn", ...), as well as deduced cases ("E20", "E21", "E2" as alias for "E2p",
                "M1", "M1-native")

            qn (tuple): quantum numbers for state

            allow_mfdn_native (bool, optional): whether or not to permit
                fallback to mfdn native value for moment,
                if not available as rme

            allow_moment_from_rme (bool, optional): whether or not to
                preferentially deduce moments from calcuated RME, if available,
                rather than using native calculated moment (which may have been
                truncated to lower precision on output)

            allow_e0_from_radius (bool, optional): whether or not to enable
                calculation of diagonal E0 rme from radius
 
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

        # TODO 07/27/22 (mac): bundle pass-through arguments into kwargs dict
        
        if verbose:
            print("get_moment {} {}".format(observable,qn))

        # apply aliases
        if observable in OBSERVABLE_ALIASES:
            observable = OBSERVABLE_ALIASES[observable]
            if verbose:
                print("    observable aliased to {}".format(observable))

        # trap deduced observables
        if (observable in {"E20","E21"}):
            # isoscalar/isovector E2
            E2p = self.get_moment(
                "E2p", qn,
                allow_mfdn_native=allow_mfdn_native, allow_moment_from_rme=allow_moment_from_rme, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            E2n = self.get_moment(
                "E2n", qn,
                allow_mfdn_native=allow_mfdn_native, allow_moment_from_rme=allow_moment_from_rme, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            if (observable =="E20"):
                value = E2p+E2n
            else:
                value = E2p-E2n
            return value
        elif (observable == "M1"):
            # physical M1
            Dsp = self.get_moment(
                "Dsp", qn,
                allow_mfdn_native=allow_mfdn_native, allow_moment_from_rme=allow_moment_from_rme, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            Dsn = self.get_moment(
                "Dsn", qn,
                allow_mfdn_native=allow_mfdn_native, allow_moment_from_rme=allow_moment_from_rme, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            Dlp = self.get_moment(
                "Dlp", qn,
                allow_mfdn_native=allow_mfdn_native, allow_moment_from_rme=allow_moment_from_rme, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            gp = 5.586
            gn = -3.826
            value = Dlp+gp*Dsp+gn*Dsn
            return value
        elif (observable == "M1-native"):
            # force retrieval of MFDn native precomputed "physical M1" linear combination
            #
            # If no RMEs are available, and only the MFDn.res static M1 moment
            # observables are available, this may yield higher precision than
            # the M1 value this function would usually compute as a linear
            # combination of the D terms taken from MFDn.res.  These terms
            # typically have smaller magnitudes but are saved at the same fixed
            # point precision as the physical M1.
            value = self.mfdn_ob_moments.get("M1",{}).get(qn,default)
            return value
        elif (observable in {"Dl0","Dl1"}):
            # isoscalar/isovector orbital M1
            Dlp = self.get_moment(
                "Dlp", qn,
                allow_mfdn_native=allow_mfdn_native, allow_moment_from_rme=allow_moment_from_rme, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            Dln = self.get_moment(
                "Dln", qn,
                allow_mfdn_native=allow_mfdn_native, allow_moment_from_rme=allow_moment_from_rme, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            if (observable =="Dl0"):
                value = Dlp+Dln
            else:
                value = Dlp-Dln
            return value
        elif (observable in {"Ds0","Ds1"}):
            # isoscalar/isovector spin M1
            Dsp = self.get_moment(
                "Dsp", qn,
                allow_mfdn_native=allow_mfdn_native, allow_moment_from_rme=allow_moment_from_rme, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            Dsn = self.get_moment(
                "Dsn", qn,
                allow_mfdn_native=allow_mfdn_native, allow_moment_from_rme=allow_moment_from_rme, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            if (observable =="Ds0"):
                value = Dsp+Dsn
            else:
                value = Dsp-Dsn
            return value

        # extract labels
        (J,gex,n) = qn

        # get moment from rme, if available
        value = np.nan
        if allow_moment_from_rme:
            if (observable in ["Dlp", "Dln","Dsp", "Dsn"]):
                rme_prefactor = math.sqrt(4*math.pi/3)/am.hat(J)*math.sqrt((J)/((J + 1)))
            elif (observable in ["E2p", "E2n"]):
                rme_prefactor = math.sqrt(16*math.pi/5)/am.hat(J)*math.sqrt((J*(2*J - 1))/((J + 1)*(2*J + 3)))
            else:
                raise(ValueError("unrecognized moment type {}".format(observable)))
            # disable allow_rme_from_moment to avoid circular call
            rme = self.get_rme(observable, (qn,qn), allow_rme_from_moment=False)
            if verbose:
                print("    rme {}".format(rme))
            value = rme_prefactor*rme

        # else revert to native static moment
        #
        # Caveat: Native static moment may be reported but with garbage value
        # (contents of uninitialized memory) if obme calculation was disabled in
        # MFDn.
        if allow_mfdn_native and np.isnan(value):
            try:
                value = self.mfdn_ob_moments[observable][qn]
                if verbose:
                    print("    mfdn moment native lookup succeeded ({})".format(value))
            except KeyError as err:
                if verbose:
                    print("    mfdn moment native lookup failed ({})".format(err))

        if np.isnan(value):
            value = default
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
        value_sqr = self.mfdn_level_properties.get(am_type+"_sqr",{}).get(qn,default)
        if (value_sqr<0):
            raise ValueError("Invalid squared angular momentum value: {:e}".format(value_sqr))
        value = tools.effective_am(value_sqr)

        return value

    def get_rme(
            self, observable, qn_pair, rank="ob",
            allow_mfdn_native=True, allow_rme_from_moment=True, allow_e0_from_radius=True,
            default=np.nan, verbose=False
    ):
        """Retrieve reduced matrix element (RME).

        Returns RME in Edmonds convention, as common for spectroscopic data
        analysis in the shell model community.

        The first choice of source for the rme is a postprocessor rme.

        For hw=0, the assumption is that non-oscillator basis was used, and
        native calculated RMEs or moments from MFDn are not meaningful.

        However, for one-body operators, in an oscillator run ("hw!=0"), the fallback
        options for deducing the RME are:
        
           - MFDn native transition RME (mfdn v15b00) -- if allow_mfdn_native
        
           - MFDn native moment (for M1/E2 diagonal matrix element) -- if allow_rme_from_moment
        
        See discussion under get_moment for the relation to moments.  Note that
        the present get_rme invokes get_moment with allow_moment_from_rme
        disabled, and get_moment invokes the present get_rme with
        allow_rme_from_moment disabled, to prevent circular referencing.

        Adds support for deduced cases: "E20" and E21" for isoscalar E2 and
        isovector E2, "M1" deduced from dipole terms.  (Any native-calculated
        "M1" from MFDn is ignored, as this is redundant to the dipole terms and
        not provided by obscalc-ob.)

        Note that the same string identifier (e.g., "E2p") might in principle be
        used for both a one-body operator and a two-body operator in the
        postprocessor output.  To avoid any such ambiguity, the particle rank is
        explicitly specified through the rank argument, and not left for the
        accessor to figure out.

        Arguments:

            observable (str): operator type ("E2p", "E2n", "Dlp", "Dln", "Dsp",
                "Dsn", ...), as well as deduced cases ("E20", "E21", "E2" as
                alias for "E2p", "M1")

            qn_pair (tuple): quantum numbers for states (qn_bra,qn_ket)

            rank (str, optional): particle rank ("ob" or "tb") of observable

            allow_mfdn_native (bool, optional): whether or not to permit
                fallback to mfdn native value for an ob rme (from mfdn v15b00),
                if not available from postprocessor

            allow_rme_from_moment (bool, optional): whether or not to permit
                fallback to mfdn native calculated moment as means of deducing
                M1 or E2 ob rme, if not available from postprocessor

            allow_e0_from_radius (bool, optional): whether or not to enable
                calculation of diagonal E0 rme from radius

            default (float, optional): default value to return for missing observable

        Returns:

            (float): observable value (or default if missing)

        """

        # TODO 07/12/22 (mac): It would be useful to be able to deduce generic
        # TBO rmes from mfdn native two-body observable expectation values,
        # similar to the fall-throughs done for the ob observables.
        # 
        # TODO 07/17/22 (mac): For diagonal E0 rme in CMF calculation, could apply
        # analytic cm correction, rather than simply suppressing.  However, this
        # is at least notionally redundant to using the calculated rms radius.
        # (There might be differences in available precision.)
        #
        # TODO 07/27/22 (mac): bundle pass-through arguments into kwargs dict

        if verbose:
            print("get_rme {} {} {}".format(observable, qn_pair, rank))


        # apply aliases
        if observable in OBSERVABLE_ALIASES:
            observable = OBSERVABLE_ALIASES[observable]
            if verbose:
                print("    observable aliased to {}".format(observable))
            
        # trap deduced observables
        if (observable in {"E20","E21"}):
            # isoscalar/isovector E2
            E2p = self.get_rme(
                "E2p", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            E2n = self.get_rme(
                "E2n", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            if (observable =="E20"):
                value = E2p+E2n
            else:
                value = E2p-E2n
            return value
        elif (observable in {"E00","E01"}):
            # isoscalar/isovector E0
            E0p = self.get_rme(
                "E0p",
                qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            E0n = self.get_rme(
                "E0n", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            if (observable =="E00"):
                value = E0p+E0n
            else:
                value = E0p-E0n
            return value
        elif (observable in {"E0p","E0n"} and qn_pair[0]==qn_pair[1]):
            # diagonal E0 -- deduce from radius or suppress due to spurious contribution
            #
            # Note: Diagonal E0 matrix element (i.e., "squared radius")
            # calculated with one-body lab-frame operator has spurious
            # contribution.  This could be corrected in the case where cm motion
            # is known 0s [see intrinsic (A7)].  However, we either deduce the
            # correct value from the two-body radius (if permitted by the option
            # allow_e0_from_radius), or else just suppress retrieval, to
            # prevent its accidental use.
            if allow_e0_from_radius:
                Z, N = self.params["nuclide"]
                qn = qn_pair[0]
                J, g, n = qn
                if observable=="E0p":
                    expectation_value = Z*self.get_radius("rp", qn, default=default)**2
                elif observable=="E0n":
                    expectation_value = N*self.get_radius("rn", qn, default=default)**2
                value = am.hat(J)*expectation_value
            else:
                # suppress diagonal E0, since has CM contamination
                value = default
            return value
        elif (observable == "M1"):
            # physical M1
            Dsp = self.get_rme(
                "Dsp", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            Dsn = self.get_rme(
                "Dsn", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            Dlp = self.get_rme(
                "Dlp", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            gp = 5.586
            gn = -3.826
            value = Dlp+gp*Dsp+gn*Dsn
            return value
        elif (observable in {"Dl0","Dl1"}):
            # isoscalar/isovector orbital M1
            Dlp = self.get_rme(
                "Dlp", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            Dln = self.get_rme(
                "Dln", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            if (observable =="Dl0"):
                value = Dlp+Dln
            else:
                value = Dlp-Dln
            return value
        elif (observable in {"Ds0","Ds1"}):
            # isoscalar/isovector spin M1
            Dsp = self.get_rme(
                "Dsp", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            Dsn = self.get_rme(
                "Dsn", qn_pair, rank=rank,
                allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                default=default, verbose=verbose,
            )
            if (observable =="Ds0"):
                value = Dsp+Dsn
            else:
                value = Dsp-Dsn
            return value


        # canonicalize labels
        (qn_pair_canonical,flipped,canonicalization_factor) = tools.canonicalize_Jgn_pair(
            qn_pair,tools.RMEConvention.kEdmonds
        )
        if (verbose):
            print("    canonicalizing: observable {}, qn_pair {} => qn_pair_canonical {}".format(observable,qn_pair,qn_pair_canonical))

        # extract labels
        (qn_bra,qn_ket) = qn_pair_canonical
        (J_bra,g_bra,n_bra) = qn_bra
        (J_ket,g_ket,n_ket) = qn_ket

        # retrieve underlying rme

        # Historical fallback sequence (prior to 07/27/22):
        #
        ## try:
        ##     if rank == "ob":
        ##         if allow_mfdn_native and (self.params.get("hw", 0) != 0):
        ##             try:
        ##                 rme = self.mfdn_ob_rmes[observable][qn_pair_canonical]
        ##             except KeyError:
        ##                 rme = self.postprocessor_ob_rmes[observable][qn_pair_canonical]
        ##         else:
        ##             rme = self.postprocessor_ob_rmes[observable][qn_pair_canonical]
        ##     elif rank == "tb":
        ##         rme = self.postprocessor_tb_rmes[observable][qn_pair_canonical]
        ## except KeyError as e:
        ##     if verbose:
        ##         print("    lookup failed ({})".format(e))
        ##     return default

        # get rme from postprocessor, if available
        rme = np.nan
        try:
            if rank == "ob":
                rme = self.postprocessor_ob_rmes[observable][qn_pair_canonical]
            elif rank == "tb":
                rme = self.postprocessor_tb_rmes[observable][qn_pair_canonical]
            if verbose:
                print("    postprocessor lookup succeeded ({})".format(rme))
        except KeyError as err:
            if verbose:
                print("    postprocessor lookup failed ({})".format(err))

        # else fall back on MFDn native transition
        if allow_mfdn_native and np.isnan(rme) and (self.params.get("hw", 0) != 0):
            try:
                rme = self.mfdn_ob_rmes[observable][qn_pair_canonical]
                if verbose:
                    print("    mfdn rme native lookup succeeded ({})".format(rme))
            except KeyError as err:
                if verbose:
                    print("    native rme lookup failed ({})".format(err))

        # else fall back on MFDn native moment
        if (
                allow_rme_from_moment and np.isnan(rme) and (self.params.get("hw", 0) != 0)
                and qn_pair[0]==qn_pair[1] and observable in ["Dlp", "Dln","Dsp", "Dsn", "E2p", "E2n"]
        ):
            qn, _ = qn_pair
            J, _, _ = qn
            if (observable in ["Dlp", "Dln","Dsp", "Dsn"]):
                rme_prefactor = math.sqrt(4*math.pi/3)/am.hat(J)*math.sqrt((J)/((J + 1)))
            elif (observable in ["E2p", "E2n"]):
                rme_prefactor = math.sqrt(16*math.pi/5)/am.hat(J)*math.sqrt((J*(2*J - 1))/((J + 1)*(2*J + 3)))
            # disable allow_moment_from_rme to avoid circular call
            moment = self.get_moment(observable, qn, allow_moment_from_rme=False, verbose=verbose)
            rme = moment / rme_prefactor
            if verbose:
                print("    mfdn moment retrieval completed (moment {} -> rme {})".format(moment, rme))
                    
        # apply phase factor to decanonicalize labels
        rme *= canonicalization_factor
        
        if verbose:
            print("    rme {:e}".format(rme))

        if np.isnan(rme):
            rme = default
        return rme

    def get_rme_matrix(
            self, observable, subspace_pair, dim_pair, rank="ob",
            allow_mfdn_native=True, allow_rme_from_moment=True, allow_e0_from_radius=True,
            default=np.nan, verbose=False
    ):
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

            See get_rme for optional pass-through arguments.

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
                qn_pair = (qn_bra, qn_ket)
                rme_matrix[bra_index,ket_index] = self.get_rme(
                    observable, qn_pair, rank=rank,
                    allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
                    default=default, verbose=verbose,
                )

        if (verbose):
            print("    rme matrix")
            print(rme_matrix)

        return rme_matrix

    def get_rtp(
            self, observable, qn_pair, rank="ob",
            allow_mfdn_native=True, allow_rme_from_moment=True, allow_e0_from_radius=True,
            default=np.nan, verbose=False
    ):
        """ Retrieve reduced transition probability (RTP).


        Arguments:

            observable (str): operator type ("E2p", ...) accepted by get_rme

            qn_pair (tuple): quantum numbers for states (qn_bra,qn_ket)

            See get_rme for optional pass-through arguments.

        Returns:

            (float): observable value (or default if missing)

        """

        # extract labels
        (qn_bra,qn_ket) = qn_pair
        (J_bra,g_bra,n_bra) = qn_bra
        (J_ket,g_ket,n_ket) = qn_ket

        # retrieve underlying rme
        rme = self.get_rme(
            observable, qn_pair, rank=rank,
            allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
            default=default, verbose=verbose,
        )

        # derive final value from rme
        rtp = 1/(2*J_ket+1)*rme**2

        return rtp

    def get_expectation_value(
            self, observable, qn, rank="ob",
            allow_mfdn_native=True, allow_rme_from_moment=True, allow_e0_from_radius=True,
            default=np.nan, verbose=False
    ):
        """ Retrieve expectation value (deduced from RME).

        Arguments:

            observable (str): operator type ("E2p", ...) accepted by get_rme

            qn (tuple): quantum numbers for state

            See get_rme for optional pass-through arguments.

        Returns:

            (float): observable value (or default if missing)

        """

        # extract labels
        (J,g,n) = qn

        # retrieve underlying rme
        rme = self.get_rme(
            observable, (qn,qn), rank=rank,
            allow_mfdn_native=allow_mfdn_native, allow_rme_from_moment=allow_rme_from_moment, allow_e0_from_radius=allow_e0_from_radius,
            default=default, verbose=verbose,
        )

        # derive final value from rme
        expectation_value = 1/am.hat(J)*rme

        return expectation_value
    
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
        update_observable_dictionary(self.mfdn_level_decompositions,other.mfdn_level_decompositions)
        update_observable_dictionary(self.mfdn_level_properties,other.mfdn_level_properties)
        update_observable_dictionary(self.mfdn_ob_moments,other.mfdn_ob_moments)
        update_observable_dictionary(self.mfdn_ob_rmes,other.mfdn_ob_rmes)
        update_observable_dictionary(self.mfdn_tb_expectations,other.mfdn_tb_expectations)
        update_observable_dictionary(self.postprocessor_ob_rmes,other.postprocessor_ob_rmes)
        update_observable_dictionary(self.postprocessor_tb_rmes,other.postprocessor_tb_rmes)

#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
