""" res.py -- import and retrieval of data from MFDn res file

    TODO: finish neatening transition phases
         -- check GT conjugation phase
         -- define setter method for transitions, to only store in
            canonically decreasing (f<i) order

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame
    5/31/15 (mac): Initiated (as mfdn_res.py).
    6/5/15 (mac): Allow user-supplied res file parser.
    6/5/15 (mac): Restructure as subpackage.
    Last modified 6/18/15.

"""

import os

import numpy as np

################################################################
# parser registry
################################################################

# global registration variables
res_format_parser = {}

def register_res_format(format_name,parser):
    """Register information for parsing res file.

    Args:
        format_name (str): name for res file format
        parser (callable): function for parsing file stream

    """

    res_format_parser[format_name] = parser


################################################################
# MFDnRunData storage object
################################################################

class MFDnRunData(object):
    """Class to store results of single MFDn run.

    Attributes:
        params (dict)
        states (dict) -- SUBJECT TO PHASEOUT
        properties (dict)
        energies (dict)
        tbo (dict)
        moments (dict)
        transitions (dict)
        orbital_occupations (dict)
        shell_occupations (dict)

    The expected keys in the params dictionary (and the corresponding expected value types) are:
        hw (float): hw value for run (from the run header)
        nuclide (tuple of int): tuple (Z,N) giving nuclide's proton and neutron numbers (from the run header)
        Nmin, Nmax (int): minimum and maximum numbers of HO quanta in calculation

    The expected keys in the states dictionary are (J,g,n) quantum
    number tuples.  The values are MFDnStateData objects.

    The expected keys in the properties dictionary are (J,g,n) quantum
    number tuples.  The values themselves are dictionaries of properties.

    The expected keys in the energies dictionary are (J,g,n) quantum
    number tuples.  The values are floats.

    The expected keys in the moments dictionary are ((J,g,n),type)
    tuples where

        type: "M1" for dipole, "E2" for quadrupole

    The values are tuples containing the moments as floats.

    The expected keys in the transitions dictionary are
    ((J,g,n)_final,(J,g,n)_initial,type) tuples where

        type: "GT" for Gamow-Teller, "M1" for dipole, "E2" for quadrupole

    The values are tuples containing the transition reduced matrix
    elements (RMEs) as floats.

    The expected keys in the orbital_occupations dictionary are
    (J,g,n) tuples.  The values are dictionaries of (n,l,j)->(np,nn).

    The expected keys in the shell_occupations dictionary are (J,g,n)
    tuples.  The values are dictionaries of N->(np,nn), where N=2n+l
    here is *zero* based (c.f., shell number N+1 in results file).

    The params and states dictionaries may or may not remain public in
    the future.  Use of accessors (e.g., get_energy) is preferred
    where possible.

    Methods:
        read_file(filename)
        get_levels() -- list of (J,g,n) quantum numbers
        get_property(qn,property) -- retrieves given property of state (e.g., isospin)
        get_energy(qn) -- retrieves energy
        has_rme(qnf,qni,op,Mj) -- indicates whether or not RME is available (in either direction)
        get_rme(qnf,qni,op,Mj) -- retrieves RME (in either direction, applying phase factor as needed)
        get_rtp(qnf,qni,op,Mj) -- retrieves RTP (in either direction, applying phase factor as needed)

    """

    ################################################################
    # constructor
    ################################################################

    def __init__(self):
        """ Initialize MFDnRunData instance.
        """

        # initialize data dictionaries
        self.params = {}
        self.states = {}

        # initialize dictionaries for lookup of data by qn
        self.properties = {}
        self.energies = {}
        self.tbo = {}
        self.moments = {}
        self.transitions = {}
        self.orbital_occupations = {}

    ################################################################
    # res file input
    ################################################################

    def read_file(self,filename,res_format,verbose=False):
        """Read result file data into MFDnRunData object.  If any
        data is duplicative of the old, the new data will overwrite
        the old.

        The parsing function should have the signature
        res_parser(self,fin,verbose).

        Args:
            filename (str): name of results file to read
            res_format (str): format of results file (i.e., MFDn version)
            verbose (bool): verbose output for debugging

        """

        if (res_format not in res_format_parser):
            raise ValueError("no parser registered for res file format {}".format(res_format))

        with open(filename,"rt") as fin:
            try:
                res_format_parser[res_format](self,fin,verbose=verbose)
            except ValueError as err:
                print("Parsing error in file {} with format {}".format(filename,res_format))
                raise

        if (verbose):
            print("After import: states {}, moments {}, transitions {}".format(len(self.states),len(self.moments),len(self.transitions)))

    ################################################################
    # accessors
    ################################################################

    def get_levels(self):
        """ Get list of (J,g,n) quantum numbers for all defined states.

        Returns:
            (list of tuple) : list of quantum numbers
        """

        raw_qn_list = list(self.states.keys())  # unsorted key list
        qn_list = sorted(raw_qn_list,key=(lambda qn : self.states[qn].energy))  # sort by energy
        return qn_list

    def get_property(self,qn,property):
        """Retrieves level energy by quantum numbers.

        Args:
            qn (tuple): state (J,g,n) quantum numbers
            property (string): dictionary key for desired property (e.g., "T" for isospin)

        Returns:
            (float): value

        """

        ## value = self.states[qn].properties[property]
        value = self.properties[qn][property]
        return value

    def get_energy(self,qn):
        """Retrieves level energy by quantum numbers.

        Args:
            qn (tuple): state (J,g,n) quantum numbers

        Returns:
            (float): energy

        """

        ## value = self.states[qn].energy
        value = self.energies[qn]
        return value

    def get_tbo(self,qn,op):
        """Retrieves two-body operator expectation value by quantum numbers and
        operator code.

        Arguments:
            qn (tuple): state (J,g,n) quantum numbers
            op (str): operator code

        Returns:
            (float): two-body operator expectation value
        """
        value = self.tbo[qn][op]
        return value

    def has_rme(self,qnf,qni,op,Mj):
        """Indicates whether or not reduced matrix element (RME) of
        transition operator is availble,
        regardless of which direction it was calculated in the data set.

        Args:
           qnf (tuple): final state (J,g,n) quantum numbers
           qni (tuple): final state (J,g,n) quantum numbers
           op (string): operator ("M1" or "E2")
           Mj (float): Mj of calculation for retrieval

        Returns:
            (bool): whether or not RME found

        """

        available = (
            ((qnf,qni,op,Mj) in self.transitions)
            or
            ((qni,qnf,op,Mj) in self.transitions)
        )

        return available

    def get_rme(self,qnf,qni,op,Mj,default=np.nan):
        """Retrieves reduced matrix element (RME) of transition operator,
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

        Args:
           qnf (tuple): final state (J,g,n) quantum numbers
           qni (tuple): final state (J,g,n) quantum numbers
           op (string): operator ("M1" or "E2")
           Mj (float): Mj of calculation for retrieval
           default (float): value to return if undefined (default: np.nan)

        Returns:
            (NumPy vector): values of RMEs for different components of operator

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
        """Retrieves reduced transition probability (RTP) of transition
        operator, regardless of which direction it was calculated in
        the data set.

        Obtains RTP

           B(op;Ji->Jf) = (2Ji+1)^(-1) |<Jf||op||Ji>|^2

        from RME provided by get_rme.

        Caution: Note that the arguments to get_rtp are in the order
        qnf-qni, i.e., preserving the order used in the notation for
        the *matrix element*, not vice versa.

        Limitations: Only supports operators supported by get_rme.

        Args:
           qnf (tuple): final state (J,g,n) quantum numbers
           qni (tuple): final state (J,g,n) quantum numbers
           op (string): operator ("M1" or "E2")
           Mj (float): Mj of calculation for retrieval
           default (float): value to return (for RME) if undefined (default: np.nan)

        Returns:
            (NumPy vector): values of RTPs for different components of operator

        """

        (Ji,_,_) = qni
        values = 1/(2*Ji+1)*self.get_rme(qnf,qni,op,Mj)**2

        return values


################################################################
# MFDnStateData storage object
################################################################

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


################################################################
# register formats
################################################################
# fails due to circularity
#
#      File ".\mfdnres\res.py", line 165, in <module>
#        from .formats import *
#      File ".\mfdnres\formats\res_parser_v14b06.py", line 282, in <module>
#        mfdnres.res.MFDnRunData.res_parser["v14b06"] = res_parser_v14b06
#    AttributeError: 'module' object has no attribute 'res'

## from .formats import *

################################################################
# test code
################################################################

if (__name__ == "__main__"):
    pass
