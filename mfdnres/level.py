""" NCCI level selection object.

    Mark A. Caprio
    University of Notre Dame

    - 07/31/22 (mac): Created, extracted from data.py.  Rename LevelSelector to Level.
    - 11/29/24 (mac): Define LevelQNNex.
"""

import numpy as np

from . import (
    analysis,
    data,
    ticks,
    tools,
)


################################################################
# Level interface
################################################################

class Level(object):
    """Select level from ResultsData spectrum.

    This is an interface class, not meant to be instantiated directly.

    For any daughter class, which we shall generically call Level,
    Level(args) constructs a level selector object.  Given a ResultsData
    object (or some daughter thereof), level_selector.select_level(results_data)
    selects a level based on the criteria specified via args, and returns the
    quantum number tuple (J,g,n) for the selected level (or None).

    Methods:

        level (ResultsData -> tuple): Retrieve level QN (or None) from given ResultsData

    Properties:

        descriptor_str (str): Text string describing level

        label_text (str): Formatted LaTeX text representing selected level

    """

    def __init__(self):
        """ Null initialize."""
        pass

    def select_level(self, results_data):
        """ Retrieve level.
        """
        return None

    @property
    def descriptor_str(self):
        """ Provide text string for use in descriptors.
        """
        return ""

    @property
    def label_text(self):
        """ Provide LaTeX label.
        """
        return ""


################################################################
# level selection by quantum numbers (trivial) -- (J,g,n)
################################################################

class LevelQN(Level):
    """Provides trivial level selector for given quantum numbers.

    Simply returns the given quantum numbers, unless the level is not found, in
    which case None is returned.

    """

    def __init__(self, qn):
        """ Initialize with given parameters.

        Arguments:

            qn (tuple): (J, g, n)
        """
        super().__init__()
        self._qn = qn

    def select_level(self, results_data):
        """ Retrieve level.
        """

        if self._qn not in results_data.levels:
            return None

        return self._qn

    @property
    def descriptor_str(self):
        """ Provide text string for use in descriptors.
        """
        text = data.qn_str(self._qn)
        ## text = "{:04.1f}-{:1d}-{:02d}".format(*self._qn)  # manual without reference to submodule data
        return text

    @property
    def label_text(self):
        """ Provide LaTeX label.
        """

        label = data.qn_text(self._qn)
        return label

################################################################
# level selection within isosopin suspace -- (J,g,n,T)
################################################################

class LevelQNT(Level):
    """Provides level selector within given isospin subspace.

    Returns level of given quantum numbers within given isospin subspace.

    Uses effective isospin from <T^2>, and bins by <T^2>, with the boundary
    between T1 and T2 occurring at [T1*(T1+1)+T2*(T2+1)]/2.  For the ideal case
    where states of isospin differing by unity undergo pure two-state mixing,
    this ensures that the crossover in identification happens at a mixing angle
    of 45 deg.

    """

    def __init__(self, qnT, debug=False):
        """ Initialize with given parameters.

        Arguments:
            qnT (tuple): (J, g, n_for_T, T)
        """
        super().__init__()
        self._qnT = qnT
        self._debug = debug
        if self._debug:
            print("Debug mode set for qnT {}".format(qnT))

    def select_level(self, results_data):
        """ Retrieve level.
        """

        # recover quantum numbers for sought level
        J, g, n_for_T, T = self._qnT
        if self._debug:
            print("Searching for qnT {}".format(self._qnT))

        # set up binning
        T_lower = T-1
        T_upper = T+1
        if T_lower < 0.:
            T_bound_lower = 0.  # assumes inclusive lower bound and positive definite result for comparison
        else:
            T_bound_lower = tools.effective_am((T_lower*(T_lower+1)+T*(T+1))/2)
        T_bound_upper = tools.effective_am((T_upper*(T_upper+1)+T*(T+1))/2)
            
        # scan for sought level
        current_n = 0
        current_n_for_T = 0
        while current_n_for_T < n_for_T:
            current_n +=1
            current_qn = (J,g,current_n)
            current_T = results_data.get_isospin(current_qn)
            if self._debug:
                print("T {} {} {}: {}".format(T,T_bound_lower,T_bound_upper,current_T))
            if np.isnan(current_T):
                return None  # ran out of levels
            if T_bound_lower <= current_T < T_bound_upper:
                current_n_for_T += 1
        return current_qn

    @property
    def descriptor_str(self):
        """ Provide text string for use in descriptors."""
        text = "{:04.1f}-{:1d}-{:02d}-T{:04.1f}".format(*self._qnT)
        return text

    @property
    def label_text(self):
        """ Provide LaTeX label.
        """

        J, g, n_for_T, T = self._qnT

        # format using J;T notation
        ## twice_T=int(2*T)
        ## T_str = "{}/2".format(twice_T) if twice_T % 2 else twice_T//2
        ## label = "{};{}".format(mfdnres.data.qn_text((J, g, n_for_T)), T_str)

        # format using J_{T=...} notation
        ## twice_J=int(2*J)
        ## J_str = "{}/2".format(twice_J) if twice_J % 2 else twice_J//2
        J_str = ticks.half_int_str(J)
        P_str = "+" if g==0 else "-"
        n_for_T_str = "{:d}".format(n_for_T)
        ## twice_T=int(2*T)
        ## T_str = "{}/2".format(twice_T) if twice_T % 2 else twice_T//2
        T_str = ticks.half_int_str(T)

        label = r"{{{}}}^{{{}}}_{{{};T={}}}".format(J_str,P_str,n_for_T_str,T_str)
        return label


################################################################
# level selection within 0hw/2hw subspace -- (J,g,n_for_Nex,Nex)
################################################################

class LevelQNNex(Level):
    """Provides level selector by Nex.

    Classifies states as nominally 0hw or 2hw based on simple P(Nex=0)/P(Nex=2)
    classification.  This counting is therefore not robust against severe mixing
    during crossings.

    """

    def __init__(self, qnNex, debug=False):
        """ Initialize with given parameters.

        Arguments:

            qnNex (tuple): (J, g, n_for_Nex, Nex)

        """
        super().__init__()
        self._qnNex = qnNex
        self._debug = debug
        if self._debug:
            print("Debug mode set for qnNex {}".format(qnNex))

    def select_level(self, results_data):
        """ Retrieve level.
        """

        # recover quantum numbers for sought level
        J, g, n_for_Nex, Nex = self._qnNex
        if not Nex in [0,2]:
            raise ValueError("Unexpected Nex ({})".format(Nex))
        if self._debug:
            print("  Searching for qnNex {}".format(self._qnNex))

        # scan for sought level
        current_n = 0
        current_n_for_Nex = 0
        while current_n_for_Nex < n_for_Nex:
            current_n +=1
            current_qn = (J,g,current_n)
            decomposition = results_data.get_decomposition("Nex", current_qn)
            if decomposition is None:  # ran out of levels
                return None
            Nex_ratio = decomposition[0]/decomposition[1]
            current_Nex = 0 if Nex_ratio>1 else 2
            if self._debug:
                print("  Nex_ratio {}: Nex {}".format(Nex_ratio, current_Nex))
            if current_Nex == Nex:
                current_n_for_Nex += 1
        if self._debug:
            print("  Selected: {}".format(current_qn))
        return current_qn

    @property
    def descriptor_str(self):
        """ Provide text string for use in descriptors."""
        text = "{:04.1f}-{:1d}-{:02d}-Nex{:1d}".format(*self._qnNex)
        return text

    @property
    def label_text(self):
        """ Provide LaTeX label.
        """

        J, g, n_for_Nex, Nex = self._qnNex

        # format using J_{...hw} notation
        J_str = ticks.half_int_str(J)
        P_str = "+" if g==0 else "-"
        n_for_Nex_str = "{:d}".format(n_for_Nex)
        Nex_str = "{:1d}".format(Nex)
        
        label = r"{{{}}}^{{{}}}_{{{};{}\hbar\omega}}".format(J_str,P_str,n_for_Nex_str,Nex_str)
        return label
    
    
################################################################
# level selection override
################################################################

class LevelOverride(Level):
    """Provides level selector which overrides another level selector for selected mesh points.

    Level descriptor and label are passed through untouched.

    Example:

        mfdnres.data.LevelOverride(
            mfdnres.data.LevelQN((0.0,0,1)),
            ("Nmax","hw"),
            {(10,15.): (0.0,0,2)},
            verbose = True,
        )

    """

    def __init__(self, base_level_selector, key_fields, qn_by_key, verbose=False):
        """Initialize with given parameters.

        Arguments:

            base_level_selector(mfdnres.data.Level): Level selector to
            use if no override applies

            key_fields (tuple of str): Names of parameters from which to construct key

            qn_by_key (dict): Mapping from key to (J,g,n) for overrides

        """
        super().__init__()
        self._base_level_selector = base_level_selector
        self._key_fields = key_fields
        self._qn_by_key = qn_by_key
        self._verbose = verbose

    def select_level(self, results_data):
        """ Retrieve level.
        """

        if self._verbose:
            print("  Mesh point {}".format(results_data.params))

        # recover quantum numbers for sought level
        key = analysis.extract_key(self._key_fields, results_data)
        if key in self._qn_by_key:
            qn = self._qn_by_key[key]
            if self._verbose:
                print("  Override {}".format(qn))
        else:
            qn = self._base_level_selector.select_level(results_data)
            if self._verbose:
                print("  Pass through {}".format(qn))
        
        # validate as existing level
        if qn not in results_data.levels:
            qn = None

        if self._verbose:
            print("key {} qn {}".format(key, qn))

        return qn

    @property
    def descriptor_str(self):
        """ Provide text string for use in descriptors."""
        text = self._base_level_selector.descriptor_str
        return text

    @property
    def label_text(self):
        """ Provide LaTeX label.
        """
        label = self._base_level_selector.label_text
        return label


