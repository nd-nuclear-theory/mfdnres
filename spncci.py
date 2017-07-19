""" spncci.py

    Result storage and access for spncci runs.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    7/9/17 (mac): Extract SpNCCIMeshPointData from res.py.
    7/15/17 (mac): Implement approximate shape invariants.

"""

import math

import numpy as np

import mfdnres.am
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

    def get_rme_matrix(self,observable,Jg_pair,verbose=False):
        """Retrieve RME matrix for observable.

        Assumes stored matrices are between (J,g) subspaces in
        canonical order.  Takes care of canonicalization on retrieval.

        Assumes matrix elementrs are in group-theory convention.

        Assumes matrices on diagonal sector are completely filled in,
        rather than stored just as upper triangles.

        ...

        """

        # determine canonicalization
        (Jg_pair_canonical,flipped,canonicalization_factor) = mfdnres.tools.canonicalize_Jg_pair(
            Jg_pair,mfdnres.tools.RMEConvention.kGroupTheory
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
        """

        
        <Jf||op||Ji>_Racah = sqrt(2*Jf+1) * <Jf||op||Ji>_gt

        TODO:
          - implement bra-ket conjugation flip
          - fail gracefully with default


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
        """
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

################################################################
# sorting
################################################################

# standard spncci sorting and tabulation key
KEY_DESCRIPTOR_NNHW = (("Nsigmamax",int),("Nmax",int),("hw",float))

################################################################
# Q-invariant analysis
################################################################

def shape_invariant_vector_q2(mesh_point,Jgex,selection_dim,intermediate_dim):
    """Generate matrix elements of Kumar-Cline q_2 operator.

    The shape invariants q_n are implicitly defined on page 698 of
    D. Cline, ARNPS 36, 683 (1986):

        q_2 = sqrt(5)*(QxQ)_0 ~ q^2

    If the requested dimensions are greater than the number of
    calculated eigenstates for any applicable J space, they will be
    adjusted downward accordingly.

    Arguments:
        mesh_point (SpNCCIMeshPointData): mesh point containing E2 RMEs
        Jgex (tuple): (J,gex) labels for subspace on which to calculate Q-invariant
        selection_dim (int): requested dimension of final matrix
        intermediate_dim (int): requested dimension for intermediate resolutions of the identity

    Returns:
        (np.array): vector of Q-invariant values on given subspace
    """
    
    # set up target matrix
    (J,gex) = Jgex
    dim_J = min(mesh_point.num_eigenvalues[(J,gex)],selection_dim)
    target_matrix = np.zeros((dim_J,dim_J))

    # accumulate target matrix
    operator = "Qintr"
    prefactor = math.sqrt(5)
    for J_bar in mfdnres.am.product_angular_momenta(J,2):
        coefficient = (
            (mfdnres.am.parity_sign(J_bar+J)*mfdnres.am.hat(J_bar))
            / (mfdnres.am.hat(J)*mfdnres.am.hat(2))
        )
        dim_J_bar = min(mesh_point.num_eigenvalues[(J_bar,gex)],intermediate_dim)
        matrix1 = mesh_point.get_rme_matrix(operator,((J,gex),(J_bar,gex)))[:dim_J,:dim_J_bar]
        matrix2 = mesh_point.get_rme_matrix(operator,((J_bar,gex),(J,gex)))[:dim_J_bar,:dim_J]
        target_matrix += prefactor * coefficient* np.dot(matrix1,matrix2)

    return np.diag(target_matrix)

# hard-coded 6j table for q3 calculation
#
# constructed with spncci/make_6j_table_q3
g_6j_table_q3 = {
    ( 0.0 ,  2.0 ,  2.0) :  +0.20000000,
    ( 0.5 ,  1.5 ,  1.5) :  +0.18708287,
    ( 0.5 ,  1.5 ,  2.5) :  +0.10000000,
    ( 0.5 ,  2.5 ,  2.5) :  -0.16329932,
    ( 1.0 ,  1.0 ,  1.0) :  +0.15275252,
    ( 1.0 ,  1.0 ,  2.0) :  +0.15275252,
    ( 1.0 ,  1.0 ,  3.0) :  +0.04364358,
    ( 1.0 ,  2.0 ,  2.0) :  -0.10000000,
    ( 1.0 ,  2.0 ,  3.0) :  -0.10690450,
    ( 1.0 ,  3.0 ,  3.0) :  +0.13997084,
    ( 1.5 ,  1.5 ,  1.5) :  +0.00000000,
    ( 1.5 ,  1.5 ,  2.5) :  -0.13363062,
    ( 1.5 ,  1.5 ,  3.5) :  -0.05345225,
    ( 1.5 ,  2.5 ,  2.5) :  +0.05832118,
    ( 1.5 ,  2.5 ,  3.5) :  +0.10497813,
    ( 1.5 ,  3.5 ,  3.5) :  -0.12371791,
    ( 2.0 ,  2.0 ,  2.0) :  -0.04285714,
    ( 2.0 ,  2.0 ,  3.0) :  +0.11428571,
    ( 2.0 ,  2.0 ,  4.0) :  +0.05714286,
    ( 2.0 ,  3.0 ,  3.0) :  -0.03499271,
    ( 2.0 ,  3.0 ,  4.0) :  -0.10101525,
    ( 2.0 ,  4.0 ,  4.0) :  +0.11167657,
    ( 2.5 ,  2.5 ,  2.5) :  +0.05832118,
    ( 2.5 ,  2.5 ,  3.5) :  -0.09914601,
    ( 2.5 ,  2.5 ,  4.5) :  -0.05832118,
    ( 2.5 ,  3.5 ,  3.5) :  +0.02061965,
    ( 2.5 ,  3.5 ,  4.5) :  +0.09671474,
    ( 2.5 ,  4.5 ,  4.5) :  -0.10235326,
    ( 3.0 ,  3.0 ,  3.0) :  -0.06415330,
    ( 3.0 ,  3.0 ,  4.0) :  +0.08748178,
    ( 3.0 ,  3.0 ,  5.0) :  +0.05832118,
    ( 3.0 ,  4.0 ,  4.0) :  -0.01116766,
    ( 3.0 ,  4.0 ,  5.0) :  -0.09258201,
    ( 3.0 ,  5.0 ,  5.0) :  +0.09489114,
    ( 3.5 ,  3.5 ,  3.5) :  +0.06598289,
    ( 3.5 ,  3.5 ,  4.5) :  -0.07835468,
    ( 3.5 ,  3.5 ,  5.5) :  -0.05773503,
    ( 3.5 ,  4.5 ,  4.5) :  +0.00465242,
    ( 3.5 ,  4.5 ,  5.5) :  +0.08876254,
    ( 3.5 ,  5.5 ,  5.5) :  -0.08876254,
    ( 4.0 ,  4.0 ,  4.0) :  -0.06599070,
    ( 4.0 ,  4.0 ,  5.0) :  +0.07106691,
    ( 4.0 ,  4.0 ,  6.0) :  +0.05685352,
    ( 4.0 ,  5.0 ,  5.0) :  +0.00000000,
    ( 4.0 ,  5.0 ,  6.0) :  -0.08528029,
    ( 4.0 ,  6.0 ,  6.0) :  +0.08362420,
    ( 4.5 ,  4.5 ,  4.5) :  +0.06513389,
    ( 4.5 ,  4.5 ,  5.5) :  -0.06513389,
    ( 4.5 ,  4.5 ,  6.5) :  -0.05582905,
    ( 4.5 ,  5.5 ,  5.5) :  -0.00341394,
    ( 4.5 ,  5.5 ,  6.5) :  +0.08211734,
    ( 4.5 ,  6.5 ,  6.5) :  -0.07924289,
    ( 5.0 ,  5.0 ,  5.0) :  -0.06386904,
    ( 5.0 ,  5.0 ,  6.0) :  +0.06021938,
    ( 5.0 ,  5.0 ,  7.0) :  +0.05474489,
    ( 5.0 ,  6.0 ,  6.0) :  +0.00597316,
    ( 5.0 ,  6.0 ,  7.0) :  -0.07924289,
    ( 5.0 ,  7.0 ,  7.0) :  +0.07545432,
    ( 5.5 ,  5.5 ,  5.5) :  +0.06242640,
    ( 5.5 ,  5.5 ,  6.5) :  -0.05608622,
    ( 5.5 ,  5.5 ,  7.5) :  -0.05364769,
    ( 5.5 ,  6.5 ,  6.5) :  -0.00792429,
    ( 5.5 ,  6.5 ,  7.5) :  +0.07662422,
    ( 5.5 ,  7.5 ,  7.5) :  -0.07213932,
    ( 6.0 ,  6.0 ,  6.0) :  -0.06092620,
    ( 6.0 ,  6.0 ,  7.0) :  +0.05256378,
    ( 6.0 ,  6.0 ,  8.0) :  +0.05256378,
    ( 6.0 ,  7.0 ,  7.0) :  +0.00943179,
    ( 6.0 ,  7.0 ,  8.0) :  -0.07423075,
    ( 6.0 ,  8.0 ,  8.0) :  +0.06920922,
    ( 6.5 ,  6.5 ,  6.5) :  +0.05943216,
    ( 6.5 ,  6.5 ,  7.5) :  -0.04952680,
    ( 6.5 ,  6.5 ,  8.5) :  -0.05150788,
    ( 6.5 ,  7.5 ,  7.5) :  -0.01060872,
    ( 6.5 ,  7.5 ,  8.5) :  +0.07203524,
    ( 6.5 ,  8.5 ,  8.5) :  -0.06659660,
    ( 7.0 ,  7.0 ,  7.0) :  -0.05797777,
    ( 7.0 ,  7.0 ,  8.0) :  +0.04688154,
    ( 7.0 ,  7.0 ,  9.0) :  +0.05048782,
    ( 7.0 ,  8.0 ,  8.0) :  +0.01153487,
    ( 7.0 ,  8.0 ,  9.0) :  -0.07001400,
    ( 7.0 ,  9.0 ,  9.0) :  +0.06424926,
    ( 7.5 ,  7.5 ,  7.5) :  +0.05657986,
    ( 7.5 ,  7.5 ,  8.5) :  -0.04455664,
    ( 7.5 ,  7.5 ,  9.5) :  -0.04950738,
    ( 7.5 ,  8.5 ,  8.5) :  -0.01226780,
    ( 7.5 ,  8.5 ,  9.5) :  +0.06814663,
    ( 7.5 ,  9.5 ,  9.5) :  -0.06212607,
    ( 8.0 ,  8.0 ,  8.0) :  -0.05524596,
    ( 8.0 ,  8.0 ,  9.0) :  +0.04249689,
    ( 8.0 ,  8.0 , 10.0) :  +0.04856787,
    ( 8.0 ,  9.0 ,  9.0) :  +0.01284985,
    ( 8.0 ,  9.0 , 10.0) :  -0.06641557,
    ( 8.0 , 10.0 , 10.0) :  +0.06019422,
    ( 8.5 ,  8.5 ,  8.5) :  +0.05397830,
    ( 8.5 ,  8.5 ,  9.5) :  -0.04065898,
    ( 8.5 ,  8.5 , 10.5) :  -0.04766915,
    ( 8.5 ,  9.5 ,  9.5) :  -0.01331273,
    ( 8.5 ,  9.5 , 10.5) :  +0.06480575,
    ( 8.5 , 10.5 , 10.5) :  -0.05842713,
    ( 9.0 ,  9.0 ,  9.0) :  -0.05277618,
    ( 9.0 ,  9.0 , 10.0) :  +0.03900848,
    ( 9.0 ,  9.0 , 11.0) :  +0.04681017,
    ( 9.0 , 10.0 , 10.0) :  +0.01368050,
    ( 9.0 , 10.0 , 11.0) :  -0.06330420,
    ( 9.0 , 11.0 , 11.0) :  +0.05680306,
    ( 9.5 ,  9.5 ,  9.5) :  +0.05163726,
    ( 9.5 ,  9.5 , 10.5) :  -0.03751769,
    ( 9.5 ,  9.5 , 11.5) :  -0.04598943,
    ( 9.5 , 10.5 , 10.5) :  -0.01397170,
    ( 9.5 , 10.5 , 11.5) :  +0.06189969,
    ( 9.5 , 11.5 , 11.5) :  -0.05530403,
    (10.0 , 10.0 , 10.0) :  -0.05055839,
    (10.0 , 10.0 , 11.0) :  +0.03616412,
    (10.0 , 10.0 , 12.0) :  +0.04520515,
    (10.0 , 11.0 , 11.0) :  +0.01420076,
    (10.0 , 11.0 , 12.0) :  -0.06058253,
    (10.0 , 12.0 , 12.0) :  +0.05391505
}


def shape_invariant_vector_q3(mesh_point,Jgex,selection_dim,intermediate_dim):
    """Generate matrix elements of Kumar-Cline q_3 operator.

    The shape invariants q_n are implicitly defined on page 698 of
    D. Cline, ARNPS 36, 683 (1986):

        q_3 = sqrt(5)*(QxQ)_0 ~ q^2

    If the requested dimensions are greater than the number of
    calculated eigenstates for any applicable J space, they will be
    adjusted downward accordingly.

    Global:
        g_6j_table_q3 (dict): table of "222" 6-j coefficients, labeled by the "bottom row"
            angular momenta (as floats) in canonical order

    Arguments:
        mesh_point (SpNCCIMeshPointData): mesh point containing E2 RMEs
        Jgex (tuple): (J,gex) labels for subspace on which to calculate Q-invariant
        selection_dim (int): requested dimension of final matrix
        intermediate_dim (int): requested dimension for intermediate resolutions of the identity

    Returns:
        (np.array): vector of Q-invariant values on given subspace
    """

    # set up target matrix
    (J,gex) = Jgex
    dim_J = min(mesh_point.num_eigenvalues[(J,gex)],selection_dim)
    target_matrix = np.zeros((dim_J,dim_J))

    # accumulate target matrix
    operator = "Qintr"
    prefactor = -math.sqrt(35/2)
    for J_bar_bar in mfdnres.am.product_angular_momenta(J,2):
        dim_J_bar_bar = min(mesh_point.num_eigenvalues[(J_bar_bar,gex)],intermediate_dim)
        matrix1 = mesh_point.get_rme_matrix(operator,((J,gex),(J_bar_bar,gex)))[:dim_J,:dim_J_bar_bar]

        # accumulate sum of matrix2*matrix3 products
        matrix23 = np.zeros((dim_J_bar_bar,dim_J))
        for J_bar in mfdnres.am.product_angular_momenta(J,2):
            if (not mfdnres.am.allowed_triangle(J_bar_bar,2,J_bar)):
                continue
            key = tuple(map(float,sorted((J,J_bar,J_bar_bar))))
            wigner_coefficient = g_6j_table_q3[key]
            coefficient = (
                (mfdnres.am.hat(J_bar)*mfdnres.am.hat(J_bar_bar))
                / mfdnres.am.hat(J)
                * wigner_coefficient
            )
            dim_J_bar = min(mesh_point.num_eigenvalues[(J_bar,gex)],intermediate_dim)
            matrix2 = mesh_point.get_rme_matrix(operator,((J_bar_bar,gex),(J_bar,gex)))[:dim_J_bar_bar,:dim_J_bar]
            matrix3 = mesh_point.get_rme_matrix(operator,((J_bar,gex),(J,gex)))[:dim_J_bar,:dim_J]
            matrix23 += prefactor * coefficient * np.dot(matrix2,matrix3)

        target_matrix += np.dot(matrix1,matrix23)

    return np.diag(target_matrix)

    
#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
