""" spncci_results_data.py

    Result storage and access for spncci runs.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    07/09/17 (mac): Extract SpNCCIMeshPointData from res.py.
    07/15/17 (mac): Implement approximate shape invariants.
    07/22/17 (mac): Move out approximate shape invariants.
    10/10/17 (mac): Rename SpNCCIMeshPointData to SpNCCIResultsData.
    09/21/19 (mac): Fill in docstring for SpNCCIResultsData.

"""

import math

import numpy as np

from . import (
    am,
    results_data,
    tools,
    )

#################################################
# SpNCCIResultsData
#################################################

class SpNCCIResultsData(results_data.ResultsData):
    """Container for results data for single spncci mesh point.

    Recall that there are inherited attributes:

        params (dict)
        energies (dict)
        num_eigenvalues (dict)
        filename (str)

    Data attributes:

        WARNING: The "g" quantum number is currently interpreted as the gex
        excitation quantum number and hard coded as 0 in the parser (9/21/19).

        spj_listing (list of tuple): table of J subspace dimensions for the run,
        as list of (subspace_index,J,dim)

        Jgex_values (list of tuple): listing of J subpsaces for the run, as
        (J,g)

        baby_spncci_listing (array): table of "BabySpNCCI" subspace dimensions
        for the run, as an array with rows having the following numpy dtype:

            [
                ("subspace_index",int),("irrep_family_index",int),
                ("Nsigmaex",int),("sigma.N",float),("sigma.lambda",int),("sigma.mu",int),
                ("Sp",float),("Sn",float),("S",float),
                ("Nex",int),("omega.N",float),("omega.lambda",int),("omega.mu",int),
                ("gamma_max",int),("upsilon_max",int),("dim",int)
            ]

        decompositions (dictionary): wave function probability decompositions

            Mapping: decomposition_type -> (J,g) -> matrix 
    
                decomposition_type (str): identifier for decomposition type
                ("Nex", "BabySpNCCI")
    
                (J,g) (tuple): (J,g) for space

                matrix (matrix): matrix representation of decomposition
                probabilities, i.e., where the probability histogram for each
                eigenstate is stored as the column vector matrix[:,n0], where n0
                is the zero-based state index

        observables (dictionary): observable RMEs on eigenbasis

            Mapping: (observable,(Jg_bra,Jg_ket)) -> matrix

                observable (str): observable identifier

                (Jg_bra,Jg_ket) (tuple): (J,g) for bra and ket spaces, respectively

                     Note that (Jg_bra,Jg_ket) is always stored canonically, and
                     the observable matrix elements for the reverse pairing is
                     obtained through the conjugation relation.

                matrix (matrix): matrix representation of observable on energy
                eigenbasis, i.e., indexed by zero-based state indices
                (n0_bra,n0_ket)

            The RMEs in an observable matrix are stored in group theoretical
            (Rose) convention, but the accessor for individual RMEs converts to
            Edmonds convention to meet the expectation of traditional nuclear
            theory users.

    Accessors:

        get_baby_spncci_subspace_label
        get_rme_matrix
        get_radius
        get_rme
        get_rtp
        get_decomposition

    """

    ########################################
    # Initializer
    ########################################
    def __init__ (self):
        """Initialize attributes as empty containers or None.

        Note: Attributes from parent type (params, energies) are implicitly
        initialized by calling the parent class's __init__.
        """
        super().__init__()
        self.Jgex_values = []
        self.spj_listing = None
        self.baby_spncci_listing = None
        self.decompositions = {}
        self.observables = {}

    ########################################
    # Accessors
    ########################################

    def get_baby_spncci_subspace_label(self,baby_spncci_subspace_index,label):
        """
        Arguments:
            subspace_index (int): subspace index
            label (str): one of the dtype labels for the basis listing
                structured array (e.g., "Nex", "omega.mu", ...)
        """
        return self.baby_spncci_listing[baby_spncci_subspace_index][label]

    def get_rme_matrix(self,observable,Jg_pair,verbose=False):
        """Retrieve RME matrix for observable.

        Assumes stored matrices are between (J,g) subspaces in
        canonical order.  Takes care of canonicalization on retrieval.

        Assumes matrix elements are in group-theory (Rose) convention.

        Assumes matrices on diagonal sector are completely filled in,
        rather than stored just as upper triangles.

        ...

        """

        # determine canonicalization
        (Jg_pair_canonical,flipped,canonicalization_factor) = tools.canonicalize_Jg_pair(
            Jg_pair,tools.RMEConvention.kRose
        )
        if (verbose):
            print("Jg_pair_canonical {} flipped {} canonicalization_factor {}".format(Jg_pair_canonical,flipped,canonicalization_factor))

        # retrieve underlying matrix
        key = (observable,Jg_pair_canonical)
        try:
            matrix = self.observables[observable][Jg_pair_canonical]
        except:
            return None

        # derive canonicalized matrix
        if (flipped):
            matrix = canonicalization_factor*matrix.transpose()

        return matrix

    def get_radius(self,radius_type,qn,default=np.nan):
        """Retrieve rms radius value.

        Note: Raw group-theory convention RME from spncci is intrinsic
        squared radius, i.e., summed over particles.

        Arguments:
           radius_type (str): radius type rp/rn/r
           qn (tuple): quantum numbers for state
           default (float,optional): default value to return for missing radius

        Returns:
           (float): rms radius

        """

        # extract labels
        (J,gex,n) = qn
        n0 = n-1

        # retrieve underlying rme
        if (radius_type=="r"):
            key = ("r2intr",(J,gex),(J,gex))
            ## if (key not in self.observables):
            ##     return np.nan
            try:
                Jg_pair = ((J,gex),(J,gex))
                sum_sqr_radius = self.get_rme_matrix("r2intr",Jg_pair)[n0,n0]
            except:
                return default
        elif (radius_type in {"rp","rn"}):
            sum_sqr_radius = default
        else:
            raise ValueError("radius type code")

        # derive final value from rme
        A = self.params["A"]
        rms_radius = math.sqrt(sum_sqr_radius)

        return rms_radius

    def get_rme(self,observable,qn_pair,default=np.nan,verbose=False):
        """Retrieve reduced matrix element (RME).

        Returns RME in Edmonds convention, as common for spectroscopic data
        analysis in the shell model community.

        <Jf||op||Ji>_Racah = sqrt(2*Jf+1) * <Jf||op||Ji>_gt

        Relies on get_rme_matrix for canonicalization of bra-ket order.

        """

        # extract labels
        (qn_bra,qn_ket) = qn_pair
        (J_bra,gex_bra,n_bra) = qn_bra
        (J_ket,gex_ket,n_ket) = qn_ket
        n0_bra = n_bra-1
        n0_ket = n_ket-1

        # retrieve underlying rme
        try:
            Jg_pair = ((J_bra,gex_bra),(J_ket,gex_ket))
            if (verbose):
                print("  Looking up rme matrix {} {} ->  {}[{}]".format(observable,qn_pair,Jg_pair,(n0_bra,n0_ket)))
            matrix = self.get_rme_matrix(observable,Jg_pair,verbose=verbose)
            rme_gt = matrix[n0_bra,n0_ket]
        except:
            return default

        # derive final value from rme
        rme_racah = math.sqrt(2*J_bra+1)*rme_gt

        return rme_racah

    def get_rtp(self,observable,qn_pair,default=np.nan):
        """ Retrieve reduced transition probability (RTP).
        """

        # extract labels
        (qn_bra,qn_ket) = qn_pair
        (J_bra,gex_bra,n_bra) = qn_bra
        (J_ket,gex_ket,n_ket) = qn_ket

        # retrieve underlying rme
        try:
            rme = self.get_rme(observable,qn_pair)
        except:
            return default

        # derive final value from rme
        rtp = 1/(2*J_ket+1)*rme**2

        return rtp

    def get_decomposition(self,decomposition_type,qn):
        """ Retrieve decomposition ("Nex", "BabySpNCCI", ...) as np.array.
        """

        # extract labels
        (J,gex,n) = qn
        n0 = n-1

        # retrieve decomposition
        try:
            decomposition = self.decompositions[decomposition_type][(J,gex)][:,n0]
        except:
            return None

        return decomposition

#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
