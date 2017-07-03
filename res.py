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


#################################################
# BaseRunData                                   #
#################################################
class BaseRunData (object):
    """
        Notes:
            A set of quantum numbers for an MFDnRunData instance are of the
                form (J, g, n).  A set of quantum numbers for an instance of 
                SpNCCIRunData are of the form (hw (J, gex, i)).
            BaseStateData should not be invoked directly.  Only instances of its children,
                MFDnStateData and SpNCCIStateData, should be created.  BaseStateData contains
                attributes, acccessors, and methods that are common to both MFDnStateData and
                SpNCCIStateData.
        Attributes:
            self.params: a dictionary.  Params holds various properties of the
                run, but the keys depend on rather it the run is MFDn of SpNCCI.
                There are only four entries in params for MFDnRunData: hw, Nmin,
                Nmax, and the tuple (Z, N).  The entries in params for SpNCCIRunData
                are all the data stored under the headings 'Space', 'Interaction', 
                and 'Mesh', which are currently nuclide, A, Nsigma0, Nsigmamax,
                N1v, Nmax, interaction, use_coulomb, and hw.
            self.energies: a dictionary.  The keys are the identifiers for a particualar
                state and the values are the ground state energy for that state.  For
                MFDnRunData, the keys are of the form (J, g, n) (or MFDnStateData.qn).  
                For SpNCCIRunData, they keys are of the form (hw, (J, gex, i)) (or
                SpNCCIStateData.qn)).  
            self.states: a dictionary.  A list of all states generated during a run.
                For MFDnRunData, each state is identified by the tuple (J, g, n).  In
                SpNCCIRunData, each state is identified by the tuple (hw, (J, gex, i)).
            self.properties:  a nested dictionary. The outer key is a tuple ((J, g, n) for
                MFDn and (hw, (J, gex, i)) for SpNCCI).  The inner dictionary are the properties
                of the state specified by the key. The properties stored for MFDn are J, g, 
                n, and T.  The properties stored by SpNCCI are J, gex, i, and hw.  
        Accessors:
            get_levels:  Takes no arguments are returns a list of all quantum numbers produced
                by the run, sorted based on the energy associated with each set of quantum numbers.
            get_property:  Takes as arguments a set of quantum numbers and the property needed as a
                string.  If the set of quantum numbers and the property are both valid entries in the
                dictionary, the method returns the value.  Otherwise, it returns None and prints a 
                message to the console.
            get_energy:  Takes as an argument a tuple of quantum numbers.  The the set of quantum numbers
                is valid, it returns the energy associated with those quantum numbers.  If the quantum
                numbers are not valid, it returns None and prints a message to the console.        
        Methods:
            read_file:  Takes as arguements a results file name, formatted as a string, and the name of 
                results file parser, also formatted as a string.  If the parser is a valid, register parser
                and the file name points to an existing file, this method makes a file pointer from the file 
                and sends it to the parser.
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
        self.states = {}
        self.properties = {}

    ########################################
    # Accessors                            #
    ########################################        
    def get_levels(self):
        """
            Arguments:
                None.
            Returned:
                qn_list: a list. A list of quantum numbers sorted by their
                    associated energy

            Returns a list of quantum numbers ((J, g, n) for MFDn or (hw, (J, gex, i))
            for SpNCCI), sorted by the ground state energy associated with set of
             quantum numbers.
        """
        raw_qn_list = list(self.states.keys())  # unsorted key list
        qn_list = sorted(raw_qn_list,key=(lambda qn : self.states[qn].energy))  # sort by energy
        return qn_list

    def get_property(self,qn,property):
        """
            Arguments:
                qn: a tuple.  (J, g, n) for MFDn or (hw, (J, gex, i)) for SpNCCI.
                    property: a string.  The name of the property that is to be returned.
            Returned:
                value: varies.  The property specified by the argument 'propery', associated
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

    def get_energies (self,qn):
        """
            Arguments:
                qn: a tuple.  (J, g, n) for MFDn or (hw, (J, gex, i)) for SpNCCI.
                    property: a string.  The name of the property that is to be returned.
            Returned:
                value: varies. If the set of quantum numbers is valid, value is set to the
                    ground state energy associated with the quantum numbers.  If the set of
                    quantum numbers is not valid, value is set to None. 

            Returns the ground state associated with the quantum numbers associated with 'qn', 
            if they are valid.  If they are not valid, it returns the None, and prints
            a message to the console.
        """
        if qn in self.energies:
            value = self.energies[qn]
        else:
            print(qn, 'is not a valid set of quantum numbers.  Returning None.')
            value = None 
        return value

    ########################################
    # Methods                              #
    ########################################
    def read_file(self,filename,res_format,verbose=False):
        """
            Arguments:
                filename: a string.  The name of the results file to be parsed and analyzed.
                res_format: a string.  The parser to be used.  Must be registered in res_format_parser
                verbose: a boolean. For debugging purposes.  Set to False by default.
            Returned:
                None.

            Parses the filename sepcified in the arguments using the parser also specified in the 
            arguments.
        """
        # Checks to make sure that the parser passed as an argument is a registered parser format
        if (res_format not in res_format_parser):
            raise ValueError("no parser registered for res file format {}".format(res_format))

        with open(filename,"rt") as fin:
            # Check to make sure the filename is valid
            try:
                res_format_parser[res_format](self,fin,verbose=verbose)
            except ValueError as err:
                print("Parsing error in file {} with format {}".format(filename,res_format))
                raise
        
        # For debugging purposes
       # if (verbose):
           # print("After import: states {}, moments {}, transitions {}".format(len(self.states),len(self.moments),len(self.transitions)))


#################################################
# SpNCCIRunData (Child of BaseRunData)          #
#################################################
class SpNCCIRunData (BaseRunData):
    """
        Child of BaseRunData
        Attributes:
            self.params: a dictionary.  Inherited from BaseRunData.
                Params holds various properties of the
                run, but the keys depend on rather it the run is MFDn of SpNCCI.
                There are only four entries in params for MFDnRunData: hw, Nmin,
                Nmax, and the tuple (Z, N).  The entries in params for SpNCCIRunData
                are all the data stored under the headings 'Space', 'Interaction', 
                and 'Mesh', which are currently nuclide, A, Nsigma0, Nsigmamax,
                N1v, Nmax, interaction, use_coulomb, and hw.
            self.energies: a dictionary.  Inherited from BaseRunData.
                The keys are the identifiers for a particualar
                state and the values are the ground state energy for that state.  For
                MFDnRunData, the keys are of the form (J, g, n) (or MFDnStateData.qn).  
                For SpNCCIRunData, they keys are of the form (hw, (J, gex, i)) (or
                SpNCCIStateData.qn)).  
            self.states: a dictionary.  Inherited from BaseRunData.
                A list of all states generated during a run.
                For MFDnRunData, each state is identified by the tuple (J, g, n).  In
                SpNCCIRunData, each state is identified by the tuple (hw, (J, gex, i)).
            self.properties:  a nested dictionary. Inherited from BaseRunData.
                The outer key is a tuple ((J, g, n) for
                MFDn and (hw, (J, gex, i)) for SpNCCI).  The inner dictionary are the properties
                of the state specified by the key. The properties stored for MFDn are J, g, 
                n, and T.  The properties stored by SpNCCI are J, gex, i, and hw.  
            self.spj_listing: a list of tuples. Stotes the information under the SpJ (listing) 
                data section.  Each tuple has the format (J, dim), where J is a float and dim is
                an int.
            self.baby_spncci_listing: a list of list.
            self.dimensions_by_omega: a dictionary
        Accessors:
            get_levels: Inherited from BaseRunData.  Takes no arguments are returns a list of all
                quantum numbers produced by the run, sorted based on the energy associated with
                each set of quantum numbers.
            get_property: Inherited from BaseRunData.  Takes as arguments a set of quantum numbers
                and the property needed as a string.  If the set of quantum numbers and the property
                are both valid entries in the dictionary, the method returns the value.  Otherwise,
                it returns None and prints a message to the console.
            get_energy: Inherited from BaseRunData.   Takes as an argument a tuple of quantum numbers.
                The the set of quantum numbers is valid, it returns the energy associated with those
                quantum numbers.  If the quantum numbers are not valid, it returns None and prints a
                message to the console.      
        Methods:
            read_file: Inherited from BaseRunData.  Takes as arguements a results file name, formatted
                as a string, and the name of results file parser, also formatted as a string.  If the
                parser is a valid, register parser and the file name points to an existing file, this
                method makes a file pointer from the file and sends it to the parser.
    """
    ########################################
    # Initializer                          #
    ########################################
    def __init__ (self):
        """
            Initializes self.params, self.energies, self.states, and self.properties from
            BaseRunData. Also initializes self.spj_listing and self.baby_spncci_listing.
        """
        super().__init__()
        self.spj_listing = []
        self.baby_spncci_listing = []
        self.dimensions_by_omega = {}

    ########################################
    # Accessors                            #
    ########################################        
    def get_basis (self):
        print('Need to implement')

    ########################################
    # Methods                              #
    ########################################


#################################################
# MFDnRunData (Child of BaseRunData)            #
#################################################
class MFDnRunData (BaseRunData):
    """
        Child of BaseRunData
        Attributes:
            self.params: a dictionary.  Inherited from BaseRunData.
                Params holds various properties of the
                run, but the keys depend on rather it the run is MFDn of SpNCCI.
                There are only four entries in params for MFDnRunData: hw, Nmin,
                Nmax, and the tuple (Z, N).  The entries in params for SpNCCIRunData
                are all the data stored under the headings 'Space', 'Interaction', 
                and 'Mesh', which are currently nuclide, A, Nsigma0, Nsigmamax,
                N1v, Nmax, interaction, use_coulomb, and hw.
            self.energies: a dictionary.  Inherited from BaseRunData.
                The keys are the identifiers for a particualar
                state and the values are the ground state energy for that state.  For
                MFDnRunData, the keys are of the form (J, g, n) (or MFDnStateData.qn).  
                For SpNCCIRunData, they keys are of the form (hw, (J, gex, i)) (or
                SpNCCIStateData.qn)).  
            self.states: a dictionary.  Inherited from BaseRunData.
                A list of all states generated during a run.
                For MFDnRunData, each state is identified by the tuple (J, g, n).  In
                SpNCCIRunData, each state is identified by the tuple (hw, (J, gex, i)).
            self.properties:  a nested dictionary. Inherited from BaseRunData.
                The outer key is a tuple ((J, g, n) for
                MFDn and (hw, (J, gex, i)) for SpNCCI).  The inner dictionary are the properties
                of the state specified by the key. The properties stored for MFDn are J, g, 
                n, and T.  The properties stored by SpNCCI are J, gex, i, and hw.  
            self.transition: a dictionary. The keys are tuples of the form (qn_final, qn_initial,
                type, Mj), where the quantum number values are tuples of the form (J, g, n).  Type
                is a string with one of the three values: 'Gt', 'M1', or 'E2'.  Gt stands for
                Gamow-Teller, while M1 represent a dipole transition and E2 represents a quadrupole
                transition.  The values of the dictionary are the transition reduced matriz elements,
                formatted as floats.
            self.tbo: a nested dictionary. The main keys are quanum number tuples of the
                form (J, g, n).  The inner dictionaries contain the keys 'rp', 'rn', and 'r',
                as well as any other observables specified in the 'Other 2-body observables' section.
                The name of these observables are found as the file names for the TBME files. 
            self.moments: a dictionary.  The keys are of the form ((J, g, n), type).  Type has
                three possible values: 'M1EM', 'M1', or 'E2'.  The values of the dictionary 
                are the moments, as floats, associated with each set of quantum numbers and type.
            self.orbital_occupations: a nested dictionary.  The main keys are quantum number
                tuples of the form (J, g, n).  The inner keys are tuples of the form (n, l, j) (NOTE:
                MFDn version 15 results files contain the value of 2*j instead of j).  The values are
                tuples of the form (np, nn).
        Accessors:
            get_levels: Inherited from BaseRunData.  Takes no arguments are returns a list of all
                quantum numbers produced by the run, sorted based on the energy associated with
                each set of quantum numbers.
            get_property: Inherited from BaseRunData.  Takes as arguments a set of quantum numbers
                and the property needed as a string.  If the set of quantum numbers and the property
                are both valid entries in the dictionary, the method returns the value.  Otherwise,
                it returns None and prints a message to the console.
            get_energy: Inherited from BaseRunData.   Takes as an argument a tuple of quantum numbers.
                The the set of quantum numbers is valid, it returns the energy associated with those
                quantum numbers.  If the quantum numbers are not valid, it returns None and prints a
                message to the console.
            get_moment:  Takes as arguments a set of quantum numbers and a type of moment, formatted
                as a string, called category.   If the tuple (qn, category) specifies a valid entry
                in the dictioanry self.moments, then the moments associated with that entry are
                returned.  If the tuple (qn, category) is not a valid entry, then None is returned. 
            get_tbo:  Takes as arguments a set of quantum numbers and an operator, formatted as a 
                string.  If the arguments qn and op are both valied entries in the 
                dictionary self.tbo, then this method returns the value associated with that entry.
                If either qn or op is invalid, the method returns None and prints a message to the
                console.
            get_orbital_occupation:  It takes as arguments qn and inner.  qn is a set of quantum
                numbers.  inner can be a tuple of (n, l, j) or has a default value of None.                             If only one argument, qn, is supplied, this method returns the dictionary
                associated with that set of quantum numbers.  If both qn and inner are specified,
                this method returns the value in the inner dictionary specified by both the
                (J, g, n) tuple and the (n, l, j) tuple.  If either qn or inner is an invalid
                entry, None is returned.
        Methods:
            read_file: Inherited from BaseRunData.  Takes as arguements a results file name, formatted
                as a string, and the name of results file parser, also formatted as a string.  If the
                parser is a valid, register parser and the file name points to an existing file, this
                method makes a file pointer from the file and sends it to the parser.
            has_rme:  Takes as arguments two sets of quantum numbers, the final and inital states of
                the transition, the type of transition as a string, and MJ.  The value returned
                indicates whether or not reduced matrix element (RME) of transition operator is
                availble, regardless of which direction it was calculated in the data set.
            get_rme:  Takes as arguments two sets of quantum numbers, the final and inital states of
                the transition, the type of transition as a string, and MJ.  Retrieves reduced
                matrix element (RME) of transition operator, regardless of which direction it was
                calculated in the data set.
            get_rtp:  Takes as arguments two sets of quantum numbers, the final and inital states of
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

            Initializes self.params, self.energies, self.states, and self.properties from
            BaseRunData.  Also initializes self.transitions, self.tbo, self.moments, and 
            self.orbital_occupations.
        """
        super().__init__()
        self.transitions = {}
        self.tbo = {}
        self.moments = {}
        self.orbital_occupations = {}
    ########################################
    # Accessors                            #
    ########################################        
    def get_moment (self, qn, category):
        """
            Arguments:
                qn: a tuple.  A set of quantum numbers of the form (J, g, n).
                category: a string.  Denotes the type of moments.  Values are either 
                    'M1EM', 'M1', or 'E2'.
            Returned:
                value: varies.  If the tuple of the form (qn, category) is a valid entry 
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
                qn: a tuple.  A set of quantum numbers in the form (J, g, n).
                op: a string.  The name of the two body operator that is needed.
            Returned:
                value: varies.  If qn is a valid set of quantum numbers and op is a 
                    valid operator associated with those quantum numbers, then value is
                    set to the value of the operator.  If either the set of quantum 
                    numbers of the operator is not a valid entry in self.tbo, value is set 
                    to None.

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
                qn: a tuple.  A set of quantum numbers of the form (J, g, n).
                inner: varies.  The default value is None.  If using the default value,
                    all members of the dictionary associated with the set of quantum numbers
                    is returned.  inner can be set to a particular tuple (n, l, j) such that only
                    the tuple (np, nn) associated with that entry in the inner dictionary is returned.
            Returned:
                value: varies. If inner and/or qn are both valid (dpending on of the default value of
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
                qnf: a tuple.   A set of quantum numbers in the form (J, g, n) that
                    represents the final state of the transition.
                qni: a tuple.  A set of quantum numbers in the form (J, g, n) that
                    represents the initial state of the tranistion.
                op: a string.  Specifies the operator.  Should be set to either 'M1' or
                    'E2'.
                Mj: a float. (ADD DESCRIPTION HERE)
            Returned:
                available: a boolean.  Its value is set depending on rather
                    or not the reduced matrix element of the transition specified
                    by the arguments is found.

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
               qnf: a tuple.   A set of quantum numbers in the form (J, g, n) that
                    represents the final state of the transition.
                qni: a tuple.  A set of quantum numbers in the form (J, g, n) that
                    represents the initial state of the tranistion.
                op: a string.  Specifies the operator.  Should be set to either 'M1' or
                    'E2'.
                Mj: a float. (ADD DESCRIPTION HERE)
                default: a float.  The value to be returned if elements are undefiend.
            Returned:
                values: a numpy vector.  Contains the values of the reduced matrix elements
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
               qnf: a tuple.   A set of quantum numbers in the form (J, g, n) that
                    represents the final state of the transition.
                qni: a tuple.  A set of quantum numbers in the form (J, g, n) that
                    represents the initial state of the tranistion.
                op: a string.  Specifies the operator.  Should be set to either 'M1' or
                    'E2'.
                Mj: a float. (ADD DESCRIPTION HERE)
                default: a float.  The value to be returned if elements are undefiend.
            Returned:
                values: a numpy vector.  The values of the reduced transition probabilities 
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
 

#################################################
# BaseStateData                                 #
#################################################
class BaseStateData (object):
    """  
        Notes:
            BaseStateData should not be invoked directly.  Only instances of its children,
                MFDnStateData and SpNCCIStateData, should be created.  BaseStateData contains
                attributes, acccessors, and methods that are common to both MFDnStateData and
                SpNCCIStateData.
        Attributes:
            self.properties: a dictionary.  Will contain properties of the state.  Only
                initialized by this class.
        Accessors:
            get_property: Takes as an argument a string of the property to be accessed.  If the 
                property is a valid key in self.properties, the method returns the value 
                associated with that key.  Otherwise, it returns None and prints a message to the 
                console.
        Methods:
    """
    ########################################
    # Initializer                          #
    ########################################
    def __init__(self):
        """
        Arguments:
                qn: a tuple.  A set of quantum numbers.  For MFDnStateData, it has the form 
                    (J, g, n).  For SpNCCIStateData, it has the form (hw, (J, gex, i)).
            Returned:
                None.

            Initializes self.properties and provides an accessor for the dictionary.
        """
        self.properties = {}


    ########################################
    # Accessors                            #
    ########################################        
    def get_property (prop):
        """
            Arguments:
                prop: a string.  The name of the property that is to be accessed.
            Returned:
                value: varies.  If prop is a valid key in self.properties, then value is
                    set to the associated value in the dictionary.  Otherwise, value is set to 
                    None.

            If prop is a valid key in the dictionary self.properties, this method returns the
            value of that property.  If prop is an invalid key, the method returns None and prints
            a message to the console.
        """
        if prop in self.properties:
            value = self.properties[prop]
        else:
            print(prop, 'is not a valid key in self.properties.  Returning None.')
            value = None
            return value
    

#################################################
# SpNCCIStateData (Child of BaseStateData)      #
#################################################
class SpNCCIStateData (BaseStateData):
    """  
        Attributes:
            self.properties: a dictionary.  Inherited by BaseStateData.  Initialized by 
                BaseStateData, but assigned values here based on the argument to the 
                initializer qn.  Keys in this dictionary should be 'hw', 'J', 'gex', and 
                'i'.
            self.hw: a float.  The value of hw that will be used to identify a instance of
                SpNCCIStateData.
            self.decomposition_Nex: a dictionary.
            self.decomposition_baby_spncci: a dictionary.
            self.observables: a dictionary.
        Accessors:
            get_property: Inherited from BaseStateData.  Takes as an argument a string
                of the property to be accessed.  If the property is a valid key in
                self.properties, the method returns the value associated with that key.
                Otherwise, it returns None and prints a message to the console.
            get_hw: Takes no arguments.  Returns the value of hw.
            get_decomposition_nex: Takes as an argument a value of J, as a float.  
                Returns the matrix from the data section 'Decompostion: Nex' associated
                with that J value.  Has a check to make sure J is a valid entry in 
                self.decomposition_nex.
            get_decomposition_baby_spncci: Takes as an argument a value of J, as a float.  
                Returns the matrix from the data section 'Decompostion: BabySpNCCI'
                associated with that J value.  Has a check to make sure J is a valid entry in 
                self.decomposition_baby_spncci.
            get_observable:  Takes three arguments: op (an string), qn_final (a tuple of
                the form (J_final, gex_final)), and qn_inital (a tuple of the form (J_initial, gex_initial)).
                If the combined tuple (op, qn_final, qn_inital) is a valid entry in 
                self.observables, then this accessor returns the matrix associated with that tuple.  If the 
                combined tuple is not a valid entry, then the method returns None.
        Methods:
    """
    ########################################
    # Initializer                          #
    ########################################
    def __init__(self, hw):
        """
            Arguments:
                hw: a float.  The hw value that will be used to identify a particular instance
                    of SpNCCIStateData.
            Returned:
                None.

            The BaseStateData initializer set.properties.  The SpNCCIStateData sets the
            value of self.hw to the hw supplied in the arguments.  It also makes 'hw' a
            property in self.properties.  SpNCCIStateData also initializes
            self.decomposition_nex, self.decomposition_baby_spncci, and self.obervables.
        """
        super().__init__()
        self.hw = hw
        self.properties['hw'] = hw
        self.decomposition_nex = {}
        self.decomposition_baby_spncci = {}
        self.observables = {}

    ########################################
    # Accessors                            #
    ########################################        
    def get_hw (self):
        """
            Arguments:
                None.
            Returned:
                value: a float.  The value of hw for the calling instance of SpNCCIStateData.

            Returns the value of hw for the instance of SpNCCIStateData that called the accessor.
        """
        value = self.hw
        return value

    def get_decomposition_nex (self, J):
        """
            Arguments:
                J: a float.  A valid J value for the instance of SpNCCIStateData.
            Returned:
                value: varies.  If J is a valid entry in self.decomposition_nex, value is
                    set to the list of lists associated with that J value.  If J is not a 
                    valid entry, value is set to None.

            If the argument, J, is a valid entry in self.decomposition_nex, this accessor returns
            the list of lists associated with that J value.  If J is not a valid entry, then None is
            returned. 
        """
        if J in self.decomposition:
           value = self.decomposition_nex[J]
        else:
           print(J, 'is not a valid entry in self.decomposition_nex.  Returning None.')
           value = None
        return value

    def get_decomposition_baby_spncci (self):
        """
            Arguments:
                J: a float.  A valid J value for the instance of SpNCCIStateData.
            Returned:
                value: varies.  If J is a valid entry in self.decomposition_baby_spncci, value is
                    set to the list of lists associated with that J value.  If J is not a 
                    valid entry, value is set to None.

            If the argument, J, is a valid entry in self.decomposition_baby_spncci, this accessor returns
            the list of lists associated with that J value.  If J is not a valid entry, then None is
            returned. 
        """
        if J in self.decomposition_baby_spncci:
            value = self.decomposition_baby_spncci[J]
        else:
            print(J, 'is not a valid entry in self.decomposition_baby_spncci.  Returning None.')
            return value
        return value

    def get_observable (self, op, qn_final, qn_initial):
        """
            Arguments:
                op: an int.  The operator name for the observable matrix.
                qn_final: a tuple.  The tuple is of the form (J, gex) and is the final state
                    of the observables matrix to be accessed.
                qn_inital:  a tuple.  The tuple is of the form (J, gex) and is the initial state
                    of the observables matrix to be accessed.
            Returned:
                value: varies.   If the combined tuple (op, qn_final, qn_initial) is a valid
                    entry in self.observables then value is set to the observable matrix associated with
                    that tuple.  If the combined tuple is not valid, then value is set to None.

            If the combined tuple made from the arguments, (op, qn_final, qn_inital), is a valid
            entry in self.observables, then the accessor returns the observables matrix associated with the tuple.
            If the combined tuple is not a valid entry, then the accessor returns None and prints a message to the 
            console.
        """
        if (op, qn_final, qn_initial) in self.observables:
            value = self.observables[(op, qn_final, qn_initial)]
        else:
            print('(', operator_index, qn_final, qn_initial, ')', 'is not a valid entry in self.observables.  Returning None.')
            value = None
        return None

#################################################
# MFDnStateData (Child of BaseStateData)        #
#################################################
class MFDnStateData (BaseStateData):
    """  
        Attributes:
            self.qn: a tuple. A set of quantum numbers. 
                For MFDnStateData, they are of the form (J, g, n) and for SpNCCIStateData,
                they are of the form (hw, (J, gex, i)).
            self.properties: a dictionary.  Inherited by BaseStateData.  Initialized by 
                BaseStateData, but assigned values here based on the argument to the 
                initializer qn.  Keys in this dictionary should be 'J', 'g', 'i', and 'T'.
            self.energy: a float.  The ground state energy
                of that state specified by the quantum numbers self.qn.
            self.obo: a dictionary.
        Accessors:
            get_qn:  Takes no arguments.
                Returns the tuple of quantum numbers for the particular
                instance of BaseStateData.
            get_energy:  Takes no arguments.
                Returns the ground state energy of the instance of BaseStateData.
            get_property: Inherited from BaseStateData.  Takes as an argument a string
                of the property to be accessed.  If the property is a valid key in
                self.properties, the method returns the value associated with that key.
                Otherwise, it returns None and prints a message to the console.
            get_obo:
        Methods:
    """
    ########################################
    # Initializer                          #
    ########################################
    def __init__(self,qn,T, E):
        """
            BasisStateData initializes self.properties.  MFDnStateStata sets the value of
            self.qn to the value of the argument, qn, and sets the value of self.energy to
            the value of the argument E.  MFDnStateData also initializes self.tbo.
            It also assigns entries to self.properties from the arguments qn and T (J, g, n, T).
        """
        super().__init__()
        self.energy = E
        self.qn = qn
        (self.properties["J"], self.properties["g"], self.properties["n"]) = qn
        self.properties['T'] = T
        self.obo = {}

    ########################################
    # Accessors                            #
    ########################################        
    def get_qn (self):
        """
            Arguments:
                None.
            Returned:
                value: a tuple.  The set of quantum numbers for the state.  

            Returns the tuple of quantum numbers for the particular instance of 
            BaseStateData.
        """
        value = self.qn
        return value

    def get_energy (self):
        """
            Arguments:
                None.
            Returned:
                value: a float.  The ground state energy for the instance of BaseStateData.

            Returns the ground state energy associated with the particular instance of 
            BaseStateData that calls the methdod.
        """
        value = self.energy
        return value


    def get_obo (self):
        print('Needs to be implemented')


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
