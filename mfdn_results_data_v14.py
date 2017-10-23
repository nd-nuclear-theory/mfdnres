""" mfdn_results_data_v14.py

    Result storage and access for mfdn runs.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    10/6/17 (mac): Extract MFDnResultsData from res.py.
    10/23/17 (mac): Add get_radius accessor.
"""

import math

import numpy as np

import mfdnres.am
import mfdnres.results_data


#################################################
# MFDnResultsDataV14
#################################################

class MFDnResultsDataV14(mfdnres.results_data.ResultsData):
    """ Container for MFDn results -- legacy version from original mfdn v14 parsers.

    Deprecated after data hierarchy rethought for spncci parser.

    TODO clean up docstring

    TODO: finish neatening transition phases
        - check GT conjugation phase
        - define setter method for transitions, to only store in
           canonically decreasing (f<i) order
        - reduce states data member to sorted list

    Attributes:
        self.params (dict): INHERITED
        self.energies (dict): INHERITED

        self.states (dict): Maps from quantum number tuple to instance of MFDnStateData.
            A list of all states generated during a run.
            For MFDnResultsData, each state is identified by the tuple (J, g, n).  In
            SpNCCIResultsData, each state is identified by the tuple (hw, (J, gex, i)).
        self.properties (nested dictionary):  Maps from quantum number to properties dictionary.
            The outer key is the quantum number tuple of the form (J, g, n).
            The inner dictionary are the properties of the state specified by the key.
            The properties stored are J, g, n, and T.  
        self.transition (dictionary): Maps from (qn_final, qn_initial, type, Mj) to transition
            reduced matrix elements.
            The keys are tuples of the form (qn_final, qn_initial, type, Mj), where the
            quantum number values are tuples of the form (J, g, n).  Type is a string
            with one of the three values: 'Gt', 'M1', or 'E2'.  GT stands for Gamow-Teller,
            while M1 represent a dipole transition and E2 represents a quadrupole transition.
            The values of the dictionary are the transition reduced matrix elements, formatted
            as floats.
        self.tbo (nested dictionary): Maps from quantum number tuple to dictionary of two
            body observables.
            The main keys are quanum number tuples of the form (J, g, n).  The inner
            dictionaries contain the keys 'rp', 'rn', and 'r', as well as any other
            observables specified in the 'Other 2-body observables' section. The name of
            these observables are found as the file names for the TBME files. 
        self.moments (dictionary): Maps from (qn, type) to a list of moments.
            The keys are of the form ((J, g, n), type).  Type has three possible values:
            'M1EM', 'M1', or 'E2'.  The values of the dictionaries are the moments, as
            floats, associated with each set of quantum numbers and type.
        self.orbital_occupations (nested dictionary):  Maps from (J, g, n) to (n, l, j) to (np, nn).
            The main keys are quantum number tuples of the form (J, g, n).
            The inner keys are tuples of the form (n, l, j) (NOTE: MFDn version
            15 results files contain the value of 2*j instead of j).  The values are
            tuples of the form (np, nn).
    Accessors:
        get_levels: INHERITED
        get_energy: INHERITED
        get_property: Takes as arguments a set of quantum numbers
            and the property needed as a string.  If the set of quantum numbers and the property
            are both valid entries in the dictionary, the method returns the value.  Otherwise,
            it returns None and prints a message to the console.
        get_moment: Accessor for moments data by the tuple (qn, type). 
            Takes as arguments a set of quantum numbers and a type of moment, formatted
            as a string, called category.   If the tuple (qn, category) specifies a valid entry
            in the dictioanry self.moments, then the moments associated with that entry are
            returned.  If the tuple (qn, category) is not a valid entry, then None is returned. 
        get_tbo: Accessor two body operator data by quantum number tuple and operator name. 
            Takes as arguments a set of quantum numbers and an operator, formatted as a 
            string.  If the arguments qn and op are both valied entries in the 
            dictionary self.tbo, then this method returns the value associated with that entry.
            If either qn or op is invalid, the method returns None and prints a message to the
            console.
        get_orbital_occupation:  Accessor for (np, nn) tuple by (J, g, n) and (n, l, j).  
            It takes as arguments qn and inner.  qn is a set of quantum
            numbers.  inner can be a tuple of (n, l, j) or has a default value of None.               
            If only one argument, qn, is supplied, this method returns the dictionary
            associated with that set of quantum numbers.  If both qn and inner are specified,
            this method returns the value in the inner dictionary specified by both the
            (J, g, n) tuple and the (n, l, j) tuple.  If either qn or inner is an invalid
            entry, None is returned.
    Methods:
        has_rme:  Checks for reduced matrix element of transition operator.  
            Takes as arguments two sets of quantum numbers, the final and inital states of
            the transition, the type of transition as a string, and MJ.  The value returned
            indicates whether or not reduced matrix element (RME) of transition operator is
            availble, regardless of which direction it was calculated in the data set.
        get_rme:  Accessor for the reduced matric element of the transition operator.
            Takes as arguments two sets of quantum numbers, the final and inital states of
            the transition, the type of transition as a string, and MJ.  Retrieves reduced
            matrix element (RME) of transition operator, regardless of which direction it was
            calculated in the data set.
        get_rtp: Accessor for the reduced transition probability of transition operator. 
            Takes as arguments two sets of quantum numbers, the final and inital states of
            the transition, the type of transition as a string, and MJ. Retrieves reduced
            transition probability (RTP) of transition operator, regardless of which direction
            it was calculated in the data set.
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
        self.states = {}
        self.properties = {}
        self.transitions = {}
        self.tbo = {}
        self.moments = {}
        self.orbital_occupations = {}

    ########################################
    # Accessors
    ########################################        

    def get_property(self,qn,property):
        """
            Arguments:
                qn (tuple): A quantum number tuple of the form  (J, g, n)
                property (string):  The name of the property that is to be returned.
            Returned:
                value (varies):  The value of the property, given a valid qn and
                    property string.
                    The property specified by the argument 'propery', associated
                    with the set of quantum numbers specified by the argument 'qn'.  If a 
                    invalid entry is specifed by the arguments, value is set to None.

            Returns the value of the property, sepcified by the argument 'property', that is
            associated with the set of quantum numbers that are specified by the argument 'qn'.
            There are checks to make sure both 'property' and 'qn' are entries in the dictionaries.
            If either argument is not an entry, a message is printed to the console, and the 
            value of value is set to None.

        """
        if qn in self.properties:
            if property in self.properties[qn]:
                value = self.properties[qn][property]
            else:
                print (property, 'is not a valid property of', qn)
                print ('Returning None.')
                value = None
        else:
            print(qn, 'is not a valid set of quantum numbers.  Returning None.')
            value = None
        return value


    def get_moment (self, qn, category):
        """
            Arguments:
                qn (tuple):  A set of quantum numbers of the form (J, g, n).
                category (string):  Denotes the type of moments.
                    Values are either 'M1EM', 'M1', or 'E2'.
            Returned:
                value (varies): Moments data, given valid qn and category.
                    If the tuple of the form (qn, category) is a valid entry 
                    in self.moments, then value is set to the moments associated with the
                    (qn, tuple) pair.  If the tuple (qn, category) is invalid, value is set 
                    to None.

            If the tuple (qn, category) specifies a valid entry in the dictioanry self.moments, 
            then the moments associated with that entry are returned.  If the tuple (qn, category)
            is not a valid entry, then None is returned. 
        """
        if (qn, category) in self.moments:
            value = self.moments[(qn, category)]
        else:
            print(qn, 'and', category, 'do not make a valid entry in the self.moments dictionary.')
            print('Returning None.')
            value = None
        return value

    def get_tbo(self,qn,op):
        """
            Arguments:
                qn (tuple):  A set of quantum numbers in the form (J, g, n).
                op (string):  The name of the two body operator that is needed.
            Returned:
                value (varies): Two body operator data, given valid qn and op.
                    If qn is a valid set of quantum numbers and op is valid operator
                    associated with those quantum numbers, then value is set to the value
                    of the operator.  If either the set of quantum numbers of the operator
                    is not a valid entry in self.tbo, value is set to None.

            If the arguments qn and op are both valied entries in the dictionary self.tbo,
            then this method returns the value associated with that entry.  If either qn or
            op is invalid, the method returns None and prints a message to the console.
        """
        if qn in self.tbo:
            if op in self.tbo[qn]:
                value = self.tbo[qn][op]
            else:
                print(op, 'is not a valid operator for', qn,'.  Returning None.')
                value = None
        else:
            print(qn, 'is not a valid set of quantum numbers.  Returning None.')
            value = None
        return value

    def get_radius(self,radius_type,qn,default=np.nan):
        """Retrieve rms radius value.

        This is a special case of get_tbo, provided for compatibility
        with the results_data classes for other codes (mfdn_v15,
        spncci, ...).  We use the value from "Radii" two-body
        observable calculation (not the one-body radius also reported
        by MFDn in oscillator runs).

        Arguments:
           radius_type (str): radius type rp/rn/r
           qn (tuple): quantum numbers for state
           default (float,optional): default value to return for missing radius

        Returns
           (float): rms radius

        """

        # extract labels
        (J,gex,n) = qn

        rms_radius = self.get_tbo(qn,radius_type)

        return rms_radius

    def get_orbital_occupation (self, qn, inner=None):
        """
            Arguments:
                qn (tuple):  A set of quantum numbers of the form (J, g, n).
                inner (varies):  Defaults to None of the tuple (n, l, j) if supplied.
                    The default value is None.  If using the default value, all members of
                    the dictionary associated with the set of quantum numbers is returned.
                    inner can be set to a particular tuple (n, l, j) such that only the tuple
                    (np, nn) associated with that entry in the inner dictionary is returned.
            Returned:
                value (varies): (np, nn) if (n, l,j) and/or (J, g, n) are valid.
                    If inner and/or qn are both valid (depending on of the default value of
                    inner is used or not), value is set to the entry in the dictionary specified by the 
                    arguments.  If at least on of the arguments are invalid, value is set to None.

            If only one argument, qn, is supplied, this method returns the dictionary associated with
            that set of quantum numbers.  If both qn and inner are specified, this method returns
            the value in the inner dictionary specified by both the (J, g, n) tuple and the (n, l, j)
            tuple.  If either qn or inner is an invalid entry, None is returned.
        """
        if inner == None:
            if qn in self.orbital_occupations:
                value = self.orbital_occupations[qn]
            else:
                print(qn, 'is not a valid set of quantum numbers for self.orbital_occupations.')
                print('Returning None.')
                value = None
        else:
            if qn in self.orbital_occupations:
                if inner in self.orbital_occupations[qn]:
                    value = self.orbital_occupations[qn][inner]
                else:
                    print(inner, 'is not a valid key for self.orbital_occupations given',
                        qn, 'as quantum numbers.')
                    print('Returning None')
                    value = None
            else:
                print(qn, 'is not a valid set of quantum numbers for self.orbital_occupations.')
                print('Returning None.')
                value = None
        return value

    def has_rme(self,qnf,qni,op,Mj):
        """
            Arguments:
                qnf (tuple):  A set of quantum numbers in the form (J, g, n) that
                    represents the final state of the transition.
                qni (tuple): A set of quantum numbers in the form (J, g, n) that
                    represents the initial state of the tranistion.
                op (string):  Specifies the operator.  
                    Should be set to either 'M1' or 'E2'.
                Mj (float): (ADD DESCRIPTION HERE)
            Returned:
                available (boolean): Availiability of reduced matrix element for transition.
                    Its value is set depending on rather  or not the reduced matrix element
                    of the transition specified by the arguments is found.

            Indicates whether or not reduced matrix element (RME) of
            transition operator is availble,
            regardless of which direction it was calculated in the data set.
        """
        available = (
            ((qnf,qni,op,Mj) in self.transitions)
            or
            ((qni,qnf,op,Mj) in self.transitions)
        )

        return available

    def get_rme(self,qnf,qni,op,Mj,default=np.nan):
        """
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
 
    def get_rtp(self,qnf,qni,op,Mj,default=np.nan):
        """
            Arguments:
               qnf (tuple):  A set of quantum numbers in the form (J, g, n) that
                    represents the final state of the transition.
                qni (tuple):  A set of quantum numbers in the form (J, g, n) that
                    represents the initial state of the tranistion.
                op (string):  Specifies the operator.
                    Should be set to either 'M1' or 'E2'.
                Mj (float): (ADD DESCRIPTION HERE)
                default (float):  The value to be returned if elements are undefiend.
            Returned:
                values (numpy vector):  The values of the reduced transition probabilities 
                    for the different components of the operator specified in the arguments.

            Retrieves reduced transition probability (RTP) of transition
            operator, regardless of which direction it was calculated in
            the data set.

            Obtains RTP
  
               B(op;Ji->Jf) = (2Ji+1)^(-1) |<Jf||op||Ji>|^2

            from RME provided by get_rme.

            Caution: Note that the arguments to get_rtp are in the order
            qnf-qni, i.e., preserving the order used in the notation for
            the *matrix element*, not vice versa.

            Limitations: Only supports operators supported by get_rme.
        """
        (Ji,_,_) = qni
        values = 1/(2*Ji+1)*self.get_rme(qnf,qni,op,Mj)**2

        return values
 

##################################################
# MFDnStateData storage object
##################################################
class MFDnStateData(object):
    """Class to store results for single MFDn calculated state.

    DEPRECATED

    Attributes:
        qn (tuple): quantum numbers (J,g,n)
        properties (dict): quantum number-ish properties ("J", "g", "n", "T")
        energy (float): energy eigenvalue
        obo (dict): miscellaneous one-body observables, always calculated with state (one-body radii)

    Note that EM moments are stored elsewhere, not as part of the state data.

    """

    def __init__(self,qn,T,energy):
        """ Initialize MFDnStateData instance.

        Args:
            qn (tuple): quantum numbers (J,g,n)
            T (float): effective isospin
            energy (float): energy eigenvalue
        """

        # initialize data dictionaries
        self.qn = qn
        self.properties = {}
        (self.properties["J"], self.properties["g"], self.properties["n"]) = qn
        self.properties["T"] = T
        self.energy = energy
        self.obo = {}

#################################################
# test code
#################################################

if (__name__ == "__main__"):
    pass
