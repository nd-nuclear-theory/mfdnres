""" res.py -- import and retrieval of data from MFDn res file

    TODO: finish neatening transition phases
         -- check GT conjugation phase
         -- define setter method for transitions, to only store in
            canonically decreasing (f<i) order

    MORE TODO:
        - reduce states data member to sorted list
        - make data members protected/private

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    5/31/15 (mac): Initiated (as mfdn_res.py).
    6/5/15 (mac): Allow user-supplied res file parser.
    6/5/15 (mac): Restructure as subpackage.
    6/29/17 (jbutler): Added in inheritance for SpNCCI, updated documentation
    Last modified 6/29/17.

"""

import os

import numpy as np

import mfdnres.make_dict

#################################################
# parser registry                               #
#################################################

# global registration variables
res_format_parser = {}

def register_res_format(format_name,parser):
    """Register information for parsing res file.

    Args:
        format_name (str): name for res file format
        parser (callable): function for parsing file stream

    """

    res_format_parser[format_name] = parser


def read_file(filename,res_format,verbose=False):
    """
        Arguments:
            filename (string): Name of results file.
                The name of the results file to be parsed and analyzed.
            res_format (string):  Name of a parser.
                The parser to be used.  Must be registered in res_format_parser
            verbose (boolean): For debugging purposes.  Set to False by default.
        Returned:
            data_instances (list): Container for MFDnRunData/SpNCCIMeshPointData instances.
                A list of of the MFDnRunData or SpNCCIMeshPointData instances generated
                from the inputted results file.

        Parses the filename sepcified in the arguments using the parser also specified in the 
        arguments.  Returns a list of MFDnRunData/SpNCCIMeshPointData instances.  If MFDn results
        files are bring ana;yzed, only one instance should be returned.  If SpNCCI results files are
        being analyzed, one instance should be returned fro each mesh point in the results file.
    """
    # Holds the instances of MFDnRunData/SpNCCIMeshPointData 
    data_instance = []

    with open(filename, 'rt') as fin:
        # Check for MFDn results files
        if res_format == 'v14b06' or res_format == 'v15b00' or res_format == 'v14b05':
            data = MFDnRunData()
            res_format_parser[res_format](data, fin, verbose=verbose)
            data_instance.append(data)

        # Check for SpNCCI results files
        elif res_format == 'spncci':
            results = []
            for row in fin:
                results.append(row)
            # Makes the dictionary before invoking the "parser" so that make_dict is only
            # called once, no matter how many mesh points are analyzed
            results_dict, order = mfdnres.make_dict.make_dict_spncci(results)
            mesh_points = results_dict['Mesh']['hw']
            for x in mesh_points:
                # hw for a SpNCCIMeshPointData instance is defined in the constructor
                data = SpNCCIMeshPointData(x)
                # The arguments are different from the MFDn parser.  res_parser_spncci takes
                # the dictionary from make_dict_spncci as an argument instead of a file pointer.
                # Reduces run time since make_dict is not called for every mesh point
                res_format_parser[res_format](data, results_dict, verbose=verbose)
                data_instance.append(data)

    return data_instance        


#################################################
# BaseRunData                                   #
#################################################
class BaseData (object):
    """
        Notes:
            A set of quantum numbers for an MFDnRunData instance are of the
                form (J, g, n).  A set of quantum numbers for an instance of 
                SpNCCIRunData are of the form (J, gex, i).
            BaseRunData should not be invoked directly.  Only instances of its children,
                MFDnRunData and SpNCCIRunData, should be created.  BaseRunData contains
                attributes, acccessors, and methods that are common to both MFDnRunData and
                SpNCCIRunData.
        Attributes:
            self.params (a dictionary):  Container for properties of run. 
                Params holds various properties of the run, but the keys depend
                on rather it the run is MFDn of SpNCCI. There are only four entries
                in params for MFDnRunData: hw, Nmin, Nmax, and the tuple (Z, N).
                The entries in params for SpNCCIRunData are all the data stored under
                the headings 'Space', 'Interaction',  and also include the hw value for
                the run, which are currently nuclide, A, Nsigma0, Nsigmamax, N1v, Nmax,
                interaction, use_coulomb, and hw.
            self.energies (dictionary):  Maps from quantum number tuple to energy.
                The keys are the identifiers for a particualar state and the values
                are the ground state energy for that state.  For MFDnRunData, the keys
                are of the form (J, g, n) (or MFDnStateData.qn).  For SpNCCIRunData,
                they keys are of the form (J, gex, i).  
        Accessors:
            get_levels:  Accessor for all quantum numbers.
                Takes no arguments are returns a list of all quantum numbers produced
                by the run, sorted based on the energy associated with each set of quantum numbers.
            get_energy:  Accessor for energy by quantum number tuple.
                Takes as an argument a tuple of quantum numbers.  The the set of quantum numbers
                is valid, it returns the energy associated with those quantum numbers.  If the quantum
                numbers are not valid, it returns None and prints a message to the console.        
        Methods:
    """
    ########################################
    # Initializer                          #
    ########################################
    def __init__ (self):
        """
            Arguments:
                None.
            Returned:
                None.

            Initializes self.params, self.energies, self.states, and self.properties.  
        """
        self.params = {}
        self.energies = {}

    ########################################
    # Accessors                            #
    ########################################        
    def get_levels(self):
        """
            Arguments:
                None.
            Returned:
                qn_list (list): A list of quantum numbers sorted by their
                    associated energy

            Returns a list of quantum numbers ((J, g, n) for MFDn or (J, gex, i)
            for SpNCCI), sorted by the ground state energy associated with set of
             quantum numbers.
        """
        # Makes a list of unsorted quantum number tuples
        raw_qn_list = list(self.states.keys())
        # Sorts the quantum numbers based on their associated energy
        qn_list = sorted(raw_qn_list,key=(lambda qn : self.states[qn].energy))
        return qn_list

    def get_energies (self,qn):
        """
            Arguments:
                qn (tuple): A tuple of quantum numbers.
                    (J, g, n) for MFDn or (hw, (J, gex, i)) for SpNCCI.
                property (string):  The name of the property that is to be returned.
            Returned:
                value(varies): The energy for a valid set of quantum numbers.
                    If the set of quantum numbers is valid, value is set to the
                    ground state energy associated with the quantum numbers.  If the set of
                    quantum numbers is not valid, value is set to None. 

            Returns the ground state associated with the quantum numbers associated with 'qn', 
            if they are valid.  If they are not valid, it returns the None, and prints
            a message to the console.
        """
        # Check to be sure the quantum numbers supplied are in self.energies
        if qn in self.energies:
            value = self.energies[qn]
        else:
            print(qn, 'is not a valid set of quantum numbers.  Returning None.')
            value = None 
        return value

    ########################################
    # Methods                              #
    ########################################

#################################################
# SpNCCIMeshPointData (Child of BaseData)          #
#################################################
class SpNCCIMeshPointData (BaseData):
    """
        Child of BaseData
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
    # Initializer                          #
    ########################################
    def __init__ (self, hw):
        """
            Arguments;
                hw (float): h-bar omega.
                    The value of h-bar omega for the instance of SpNCCIMeshPointData.
            Returned:
                None.

            Initializes self.params, self.energies from BaseRunData. Also initializes
            self.decomposition, self.observables, self.spj_listing and self.baby_spncci_listing.
            Initialized self.hw to the value given in the arguments.
        """
        super().__init__()
        self.hw = float(hw)
        self.spj_listing = []
        self.baby_spncci_listing = []
        self.dimensions_by_omega = {}
        self.decomposition = {}
        self.observables = {}

    ########################################
    # Accessors                            #
    ########################################        
    def get_basis (self):
        print('Need to implement')

    ########################################
    # Methods                              #
    ########################################


#################################################
# MFDnRunData (Child of BaseData)            #
#################################################
class MFDnRunData (BaseData):
    """
        Child of BaseData
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
            self.energies (dictionary): Maps from quantum number tuple to energy.
                Inherited from BaseRunData.
                The keys are the identifiers for a particualar
                state and the values are the ground state energy for that state.  For
                MFDnRunData, the keys are of the form (J, g, n) (or MFDnStateData.qn).  
                For SpNCCIRunData, they keys are of the form (hw, (J, gex, i)) (or
                SpNCCIStateData.qn)).  
            self.states (dictionary): Maps from quantum number tuple to instance of MFDnStateData.
                A list of all states generated during a run.
                For MFDnRunData, each state is identified by the tuple (J, g, n).  In
                SpNCCIRunData, each state is identified by the tuple (hw, (J, gex, i)).
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
            get_levels:  Accessor for all quantum numbers. 
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
            get_rtp: Accessot for the reduced transition probability of transition operator. 
                Takes as arguments two sets of quantum numbers, the final and inital states of
                the transition, the type of transition as a string, and MJ. Retrieves reduced
                transition probability (RTP) of transition operator, regardless of which direction
                it was calculated in the data set.
    """
    ########################################
    # Initializer                          #
    ########################################
    def __init__ (self):
        """
            Arguments:
                None.
            Returned:
                None.

            Initializes self.params, self.energies from BaseRunData.  Also initializes
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
    # Accessors                            #
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


    ########################################
    # Methods                              #
    ######################################## 
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
# MFDnStateData storage object                   #
##################################################
class MFDnStateData(object):
    """Class to store results for single MFDn calculated state.

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
# register formats                              #
#################################################
# fails due to circularity
#
#      File ".\mfdnres\res.py", line 165, in <module>
#        from .formats import *
#      File ".\mfdnres\formats\res_parser_v14b06.py", line 282, in <module>
#        mfdnres.res.MFDnRunData.res_parser["v14b06"] = res_parser_v14b06
#    AttributeError: 'module' object has no attribute 'res'

## from .formats import *

#################################################
# test code                                     #
#################################################

if (__name__ == "__main__"):
    pass
