""" spncci.py

    

    TODO: finish neatening transition phases
         -- check GT conjugation phase
         -- define setter method for transitions, to only store in
            canonically decreasing (f<i) order

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    7/9/17 (mac): Extract

"""

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
    def get_basis (self):
        print('Need to implement')

    def get_rme(self,qnf,qni,op,Mj,default=np.nan):
        """

        !!!!!!!!!!!!!!!! WIP !!!!!!!!!!!!!!!!

            Arguments:
                qnf (tuple):   A set of quantum numbers in the form (J, g, n) that
                    represents the final state of the transition.
                qni (tuple):  A set of quantum numbers in the form (J, g, n) that
                    represents the initial state of the tranistion.
                op (string):  Specifies the operator.
                    Should be set to either 'M1' or 'E2'.
                Mj (float): (ADD DESCRIPTION HERE)
                default (float):  The value to be returned if elements are undefiend.
            Returned:
                values (numpy vector):  Contains the values of the reduced matrix elements
                    for different components of the operator.

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
        if (op not in {"M1","E2"}):
            raise ValueError("get_rme supports only M1 and E2 at present")

        if ((qnf,qni,op,Mj) in self.transitions):
            # available as direct entry
            values = np.array(self.transitions[(qnf,qni,op,Mj)])
        elif ((qni,qnf,op,Mj) in self.transitions):
            # available as reversed entry
            values = np.array(self.transitions[(qni,qnf,op,Mj)])
            (Ji,_,_) = qni
            (Jf,_,_) = qnf
            values *= (-1)**(Jf-Ji)
        else:
            # fall through to default
            if (op == "M1"):
                values = default*np.ones(4)
            elif (op == "E2"):
                values = default*np.ones(2)

        return values
 
    ########################################
    # Methods
    ########################################


#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
