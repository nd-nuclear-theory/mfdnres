""" spncci.py

    

    TODO: finish neatening transition phases
         -- check GT conjugation phase
         -- define setter method for transitions, to only store in
            canonically decreasing (f<i) order

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    7/9/17 (mac): Extract SpNCCIMeshPointData from res.py.

"""

import math

import numpy as np

import mfdnres.res

#################################################
# SpNCCIMeshPointData (Child of BaseResultsData)
#################################################
class SpNCCIMeshPointData(mfdnres.res.BaseResultsData):
    """
        Child of BaseResultsData
        Attributes:
            self.params (dictionary):  Container for properties of run.
                Inherited from BaseRunData.
                Params holds various properties of the
                run, but the keys depend on rather it the run is MFDn of SpNCCI.
                There are only four entries in params for MFDnRunData: hw, Nmin,
                Nmax, and the tuple (Z, N).  The entries in params for SpNCCIRunData
                are all the data stored under the headings 'Space', 'Interaction', 
                and 'Mesh', which are currently nuclide, A, Nsigma0, Nsigmamax,
                N1v, Nmax, interaction, use_coulomb, and hw.
            self.energies (dictionary):  Maps from quantum number tuple to energy.
                Inherited from BaseRunData.
                The keys are the identifiers for a particualar
                state and the values are the ground state energy for that state.  For
                MFDnRunData, the keys are of the form (J, g, n) (or MFDnStateData.qn).  
                For SpNCCIRunData, they keys are of the form (hw, (J, gex, i)) (or
                SpNCCIStateData.qn)).  
            self.spj_listing (list of tuples): List of tuples of the form (J, dim).
                Stores the information under the SpJ (listing) 
                data section.  Each tuple has the format (J, dim), where J is a float and dim is
                an int.
            self.baby_spncci_listing (list of list):
            self.dimensions_by_omega (dictionary):
            self.decompositions (dictionary):
            self.observables (dictionary):
        Accessors:
            get_levels: Accessor for all quantum numbers.
                Inherited from BaseRunData.
                Takes no arguments are returns a list of all
                quantum numbers produced by the run, sorted based on the energy associated with
                each set of quantum numbers.
            get_energy: Accessor for energy by quantum number tuple.
                Inherited from BaseRunData.
                Takes as an argument a tuple of quantum numbers.
                The the set of quantum numbers is valid, it returns the energy associated with those
                quantum numbers.  If the quantum numbers are not valid, it returns None and prints a
                message to the console.      
        Methods:
       
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
        self.num_eigenvalues = {}
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

    def get_radius(self,radius_type,qn,default=np.nan):
        """
        Note: Raw gt-convention RME is intrinsic squared radius, i.e., summed over particles.

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
                sum_sqr_radius = self.observables[("r2intr",(J,gex),(J,gex))][n0,n0]
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
        

    def get_rme(self,observable,qn_bra,qn_ket,default=np.nan):
        """

        
        <Jf||op||Ji>_Edmonds = sqrt(2*Jf+1) * <Jf||op||Ji>_gt

        TODO:
          - implement bra-ket conjugation flip
          - fail gracefully with default

            Retrieves reduced matrix element (RME) of transition operator,
            regardless of which direction it was calculated in the data set.

            Obtains RME <Jf||op||Ji>, using relation
                <Jf||op||Ji> = (-)^(Jf-Ji) <Ji||op||Jf>
            which applies to both the M1 and E2 operators under Condon-Shortley phase
            conventions (Suhonen Ch. 6) for these operators.  Derived from
            W-E theorem, symmetry of CG coefficient, and conjugation
            properties of these operators.

            Note that for MFDn the reference state is the "final" state.

            Limitations: Conjugation relations need to be checked for
            other operators (e.g., GT).

        """

        # extract labels
        (J_bra,gex_bra,n_bra) = qn_bra
        (J_ket,gex_ket,n_ket) = qn_ket
        n0_bra = n_bra-1
        n0_ket = n_ket-1

        # retrieve underlying rme
        try:
            rme_gt = self.observables[(observable,(J_bra,gex_bra),(J_ket,gex_ket))][n0_bra,n0_ket]
        except:
            return default

        # derive final value from rme
        rme_edmonds = math.sqrt(2*J_bra+1)*rme_gt

        return rme_edmonds

    def get_rtp(self,observable,qn_bra,qn_ket,default=np.nan):
        """
        """ 

        # extract labels
        (J_bra,gex_bra,n_bra) = qn_bra
        (J_ket,gex_ket,n_ket) = qn_ket

        # retrieve underlying rme
        try:
            rme = self.get_rme(observable,qn_bra,qn_ket)
        except:
            return default

        # derive final value from rme
        rtp = 1/(2*J_ket+1)*rme**2

        return rtp

    def get_decomposition(self,decomposition_type,qn):
        """
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
    


    ########################################
    # Methods
    ########################################


#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
