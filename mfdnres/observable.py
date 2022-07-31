""" NCCI observable extraction object.

    Mark A. Caprio
    University of Notre Dame

    - 07/26/22 (mac): Created, from stub code in data.py.
"""


import numpy as np
import pandas as pd

from . import (
    analysis,
    data,
)

###############################################################
# Observable interface
################################################################

class Observable(object):
    """Object to extract observable from ResultsData.

    This is an interface class.

    For any daughter class, which we shall generically call ObservableExtractor,
    ObservableExtractor(args) constructs an observable extractor object.  Given
    a ResultsData object (or some daughter thereof),
    observable_extractor.observable(results_data) returns the numerical
    observable data (or np.nan).

    Methods:

        data (list of ResultsData -> pd.DataFrame): Retrieve observable value
        mesh from given ResultsData mesh

        value (ResultsData -> float OR np.array): Retrieve observable value (or
        np.nan) from given ResultsData



    Properties:

        descriptor_str (str): Text string describing observable

        observable_label_text (str): Formatted LaTeX text representing
            observable, to be interpreted in math mode

        axis_label_text (str, str): Formatted LaTeX text representing axis label

            observable (str): observable label string, to be interpreted in math mode
            units (str): units string, to be interpreted in math mode, or None

    """

    def __init__(self):
        """ Null initialize."""
        pass

    def data(self, mesh_data, key_descriptor, verbose=False):
        """Extract data frame of observable values over mesh.

        To be overridden by child object if object does not identify with a
        single _nuclide and provide a value extractor.

        Arguments:

            key_descriptor (tuple of tuple, optional): dtype descriptor for key,
            e.g., (("Nmax", int), ("hw", float))

        """
        # select mesh down to given nuclide
        selector = {"nuclide": self._nuclide}
        mesh_data_selected = analysis.selected_mesh_data(mesh_data,selector)

        # generate hw table
        #
        # e.g., key_descriptor = (("Nmax",int),("hw",float))
        extractor = lambda results_data : self.value(results_data)  # is there a better way to return an instance method as a callable?
        table = analysis.make_obs_table(mesh_data_selected, key_descriptor, extractor, verbose=verbose)

        # structure table into DataFrame with compound index
        keys = [key for key, _ in key_descriptor]
        observable_data = pd.DataFrame(table).set_index(keys)
        
        return observable_data

    ## def value(self, results_data):
    ##     """ Extract observable value.
    ##     """
    ##     pass

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return ""

    @property
    def nuclide_label_text(self):
        """Formatted LaTeX text representing nuclide.

        To be overridden by child object if object does not identify with a
        single _nuclide.

        """
        try:
            return data.isotope(self._nuclide)
        except AttributeError:
            return ""

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        return ""

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_str = r""
        units_str = None
        return observable_str, units_str


################################################################
# legacy observable support (WIP)
################################################################

# Once this is implemented, the implementation in mfdnres.data can be cleaned
# up to simply wrap any traditional tuple with LegacyObservable.

class LegacyObservable(Observable):
    """ Provides support for legacy observable (with tuple identifier).

    """

    # WIP
    
    def __init__(self, specifier):
        """Initialize with given parameters.

        Arguments:

            specifier (tuple): legacy observable specifier
            (a.k.a. "nuclide_observable")

        """
        super().__init__()
        self._specifier = specifier

        
################################################################
# deduced observable: Difference
################################################################

class Difference(Observable):
    """ Observable extractor for difference of observables.

    """

    def __init__(self, observable1, observable2, observable_label_delimiters=None):
        """ Initialize with given parameters.

        Arguments:

            observable1, observable2 (Observable): first and second terms

        """
        super().__init__()
        self._arguments = observable1, observable2
        self._observable_label_delimiters = observable_label_delimiters

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        data_meshes = (
            self._arguments[0].data(mesh_data, key_descriptor, verbose=verbose),
            self._arguments[1].data(mesh_data, key_descriptor, verbose=verbose),
            )
        return data_meshes[0] - data_meshes[1]

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "diff"
        return r"{}_{}_{}".format(
            arithmetic_operation,
            self._arguments[0].descriptor_str,
            self._arguments[1].descriptor_str,
        )

    @property
    def nuclide_label_text(self):
        """Formatted LaTeX text representing nuclide.
        """
        nuclide_label_texts = self._arguments[0].nuclide_label_text, self._arguments[1].nuclide_label_text
        same_nuclide = nuclide_label_texts[0] == nuclide_label_texts[1]
        if same_nuclide:
            label = nuclide_label_texts[0]
        else:
            label = r"{}/{}".format(
                nuclide_label_texts[0],
                nuclide_label_texts[1],
            )
        return label

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_label_texts = self._arguments[0].observable_label_text, self._arguments[1].observable_label_text
        arithmetic_operation_sign = "-"
        observable_label_delimiters = self._observable_label_delimiters
        if observable_label_delimiters is None:
            observable_label_delimiters = (("",""),("",""))
        return r"{}{}{}{}{}{}{}".format(
            observable_label_delimiters[0][0],
            observable_label_texts[0],
            observable_label_delimiters[0][1],
            arithmetic_operation_sign,
            observable_label_delimiters[1][0],
            observable_label_texts[1],
            observable_label_delimiters[1][1],
        )

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        axis_label_texts = self._arguments[0].axis_label_text, self._arguments[1].axis_label_text
        same_axis = axis_label_texts[0] == axis_label_texts[1]
        if same_axis:
            return r"\Delta {}".format(axis_label_texts[0][0]), axis_label_texts[0][1]


################################################################
# deduced observable: Ratio
################################################################

class Ratio(Observable):
    """ Observable extractor for ratio of observables.

    """

    def __init__(self, observable1, observable2, observable_label_delimiters=None):
        """ Initialize with given parameters.

        Arguments:

            observable1, observable2 (Observable): first and second terms

        """
        super().__init__()
        self._arguments = observable1, observable2
        self._observable_label_delimiters = observable_label_delimiters

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        data_meshes = (
            self._arguments[0].data(mesh_data, key_descriptor, verbose=verbose),
            self._arguments[1].data(mesh_data, key_descriptor, verbose=verbose),
            )
        return data_meshes[0] / data_meshes[1]

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "ratio"
        return r"{}_{}_{}".format(
            arithmetic_operation,
            self._arguments[0].descriptor_str,
            self._arguments[1].descriptor_str,
        )

    @property
    def nuclide_label_text(self):
        """Formatted LaTeX text representing nuclide.
        """
        nuclide_label_texts = self._arguments[0].nuclide_label_text, self._arguments[1].nuclide_label_text
        same_nuclide = nuclide_label_texts[0] == nuclide_label_texts[1]
        if same_nuclide:
            label = nuclide_label_texts[0]
        else:
            label = r"{}/{}".format(
                nuclide_label_texts[0],
                nuclide_label_texts[1],
            )
        return label

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_label_texts = self._arguments[0].observable_label_text, self._arguments[1].observable_label_text
        arithmetic_operation_sign = "/"
        observable_label_delimiters = self._observable_label_delimiters
        if observable_label_delimiters is None:
            observable_label_delimiters = (("",""),("",""))
        return r"{}{}{}{}{}{}{}".format(
            observable_label_delimiters[0][0],
            observable_label_texts[0],
            observable_label_delimiters[0][1],
            arithmetic_operation_sign,
            observable_label_delimiters[1][0],
            observable_label_texts[1],
            observable_label_delimiters[1][1],
        )

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        axis_label_texts = self._arguments[0].axis_label_text, self._arguments[1].axis_label_text
        same_axis = axis_label_texts[0] == axis_label_texts[1]
        if same_axis:
            return r"\mathrm{Ratio}", None

################################################################
# deduced observable: FixSignTo
################################################################

class FixSignTo(Observable):
    """ Observable extractor for ratio of observables.

    """

    def __init__(self, observable1, observable2, observable_label_delimiters=None):
        """ Initialize with given parameters.

        Arguments:

            observable1, observable2 (Observable): first and second terms

        """
        super().__init__()
        self._arguments = observable1, observable2
        self._observable_label_delimiters = observable_label_delimiters

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        data_meshes = (
            self._arguments[0].data(mesh_data, key_descriptor, verbose=verbose),
            self._arguments[1].data(mesh_data, key_descriptor, verbose=verbose),
            )
        return data_meshes[0] * data_meshes[1].apply(np.sign, raw=True)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return self._arguments[0].descriptor_str

    @property
    def nuclide_label_text(self):
        """Formatted LaTeX text representing nuclide.
        """
        return self._arguments[0].nuclide_label_text

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        return self._arguments[0].observable_label_text

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        return self._arguments[0].axis_label_text


################################################################
# observable: Energy
################################################################

class Energy(Observable):
    """ Observable extractor for level energy.

    """

    def __init__(self, nuclide, level):
        """ Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            level (LevelSelector): level

        """
        super().__init__()
        self._nuclide = nuclide
        self._level = level

    def value(self, results_data):
        """ Extract observable.
        """
        qn = self._level.select_level(results_data)
        return results_data.get_energy(qn)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_str(self._nuclide),
            "energy",
            self._level.descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_text = r"E"
        level_text = self._level.label_text
        label = r"{}({})".format(observable_text,level_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"E"
        units_text = r"\mathrm{MeV}"
        return observable_text, units_text

    
################################################################
# observable: Isospin
################################################################

class Isospin(Observable):
    """ Observable extractor for level energy.

    """

    def __init__(self, nuclide, level):
        """ Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            level (LevelSelector): level

        """
        super().__init__()
        self._nuclide = nuclide
        self._level = level

    def value(self, results_data):
        """ Extract observable.
        """
        qn = self._level.select_level(results_data)
        return results_data.get_isospin(qn)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_str(self._nuclide),
            "isospin",
            self._level.descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_text = r"\bar{T}"
        level_text = self._level.label_text
        label = r"{}({})".format(observable_text,level_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"\bar{T}"
        units_text = None
        return observable_text, units_text

    
################################################################
# observable: LevelIndex
################################################################

class LevelIndex(Observable):
    """ Observable extractor for level index n.

    This observable provides a diagnostic on level selection.

    """

    def __init__(self, nuclide, level):
        """ Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            level (LevelSelector): level

        """
        super().__init__()
        self._nuclide = nuclide
        self._level = level

    def value(self, results_data):
        """ Extract observable.
        """
        qn = self._level.select_level(results_data)
        J, g, n = qn
        return n

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_str(self._nuclide),
            "n",
            self._level.descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_text = r"n"
        level_text = self._level.label_text
        label = r"{}({})".format(observable_text,level_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"n"
        units_text = None
        return observable_text, units_text
    
################################################################
# observable: Radius
################################################################

RADIUS_STR_BY_OPERATOR = {
    "rp" : r"r_p",
    "rn" : r"r_n",
    "r" : "r",
    "rp-ss" : r"r_{p,\mathrm{s.s.}}",
    "rn-ss" : r"r_{n,\mathrm{s.s.}}",
}

class Radius(Observable):
    """ Observable extractor for r.m.s. radius.

    """

    def __init__(self, nuclide, operator, level):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): identifier for radius operator, as accepted by
            MFDnResultsData.get_radius()

            level (LevelSelector): level

        """
        super().__init__()
        self._nuclide = nuclide
        self._operator = operator
        self._level = level

    def value(self, results_data):
        """ Extract observable.
        """
        qn = self._level.select_level(results_data)
        return results_data.get_radius(self._operator, qn)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_str(self._nuclide),
            "radius",
            self._operator,
            self._level.descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_text = RADIUS_STR_BY_OPERATOR[self._operator]
        level_text = self._level.label_text
        label = r"{}({})".format(observable_text,level_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"r"
        units_text = r"\mathrm{fm}"
        return observable_text, units_text

    
################################################################
# observable: RadiusSqr
################################################################

# TODO

# or replace with special ratio observables


################################################################
# observable: RadiusQuart
################################################################

# TODO

# or replace with special ratio observables


################################################################
# observable: Moment
################################################################


DIPOLE_ALIASES = {
    # sensible names for M1 components (a.k.a. "dipole terms")
    "Dlp": "M1lp",
    "Dln": "M1ln",
    "Dsp": "M1sp",
    "Dsn": "M1sn",
}

class Moment(Observable):
    """ Observable extractor for electromagnetic moment.

    """

    def __init__(self, nuclide, operator, level):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): identifier for electromagnetic operator, as accepted
            by MFDnResultsData.get_moment()

            level (LevelSelector): level

        """
        super().__init__()
        self._nuclide = nuclide
        self._operator = operator
        self._level = level

        # fix up "dipole term" names to standard multipole notation
        if self._operator in DIPOLE_ALIASES:
            self._operator = self._operator[self._operator]
                              
    def value(self, results_data):
        """ Extract observable.
        """
        qn = self._level.select_level(results_data)
        return results_data.get_moment(self._operator, qn)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_str(self._nuclide),
            "moment",
            self._operator,
            self._level.descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        if self._operator in {"M1", "M1lp","M1ln","M1sp","M1sn","M1l0","M1l1","M1s0","M1s1"}:
            observable_text = r"\mu_{{{}}}".format(self._operator[2:])
        elif self._operator in {"E2p","E2n","E20","E21","E2"}:
            observable_text = r"Q_{{{}}}".format(self._operator[2:])
        level_text = self._level.label_text
        label = r"{}({})".format(observable_text,level_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        if self._operator in {"M1", "M1lp","M1ln","M1sp","M1sn","M1l0","M1l1","M1s0","M1s1"}:
            observable_text = r"\mu"
            units_text = r"\mu_N"
        elif self._operator in {"E2p","E2n","E20","E21","E2"}:
            observable_text = r"Q"  ## r"eQ"
            units_text = r"\mathrm{fm}^{2}"  ## r"e\,\mathrm{fm}^{2}"
        return observable_text, units_text

    
################################################################
# observable: MomentSqr
################################################################

# TODO

# or replace with special ratio observables


################################################################
# observable: RME
################################################################

class RME(Observable):
    """ Observable extractor for electromagnetic reduced matrix element.

    """

    # TODO restore E0 and E1 functionality
    
    def __init__(self, nuclide, operator, levelf, leveli):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): identifier for electromagnetic operator, as accepted
            by MFDnResultsData.get_rme()

            level (LevelSelector): level

        """
        super().__init__()
        self._nuclide = nuclide
        self._operator = operator
        self._level_pair = levelf, leveli

        # fix up "dipole term" names to standard multipole notation
        if self._operator in DIPOLE_ALIASES:
            self._operator = self._operator[self._operator]
                              
    def value(self, results_data):
        """ Extract observable.
        """
        qn_pair = self._level_pair[0].select_level(results_data), self._level_pair[1].select_level(results_data)
        return results_data.get_rme(self._operator, qn_pair)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_str(self._nuclide),
            "rme",
            self._operator,
            self._level_pair[0].descriptor_str,
            self._level_pair[1].descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        if self._operator in {"M1", "M1lp","M1ln","M1sp","M1sn","M1l0","M1l1","M1s0","M1s1"}:
            ## observable_text = r"M1_{{{}}}".format(self._operator[-2:])
            observable_text = r"M_{{1{}}}".format(self._operator[2:])
        elif self._operator in {"E2p","E2n","E20","E21","E2"}:
            ## observable_text = r"Q_{{2{}}}".format(self._operator[2:])
            observable_text = r"E2_{{{}}}".format(self._operator[2:])
        level_pair_text = self._level_pair[0].label_text, self._level_pair[1].label_text
        label = r"\langle {} \Vert \mathcal{{M}}({}) \Vert {} \rangle".format(level_pair_text[0],observable_text,level_pair_text[1])  # <f|O|i> = i->f
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        if self._operator in {"M1", "M1lp","M1ln","M1sp","M1sn","M1l0","M1l1","M1s0","M1s1"}:
            ## observable_text = r"\langle M_1 \rangle"
            observable_text = r"\langle \mathcal{{M}}(M1) \rangle"
            units_text = r"\mu_N"
        elif self._operator in {"E2p","E2n","E20","E21","E2"}:
            ## observable_text = r"\langle Q_2 \rangle"
            observable_text = r"\langle \mathcal{{M}}(E2) \rangle"
            units_text = r"e\,\mathrm{fm}^{2}"
        return observable_text, units_text

    

################################################################
# observable: RTP
################################################################

class RTP(Observable):
    """Observable extractor for electromagnetic reduced transition probability.

    """

    # TODO restore E0 and E1 functionality
    
    def __init__(self, nuclide, operator, levelf, leveli):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): identifier for electromagnetic operator, as accepted
            by MFDnResultsData.get_rme()

            level (LevelSelector): level

        """
        super().__init__()
        self._nuclide = nuclide
        self._operator = operator
        self._level_pair = levelf, leveli

        # fix up "dipole term" names to standard multipole notation
        if self._operator in DIPOLE_ALIASES:
            self._operator = self._operator[self._operator]
                              
    def value(self, results_data):
        """ Extract observable.
        """
        qn_pair = self._level_pair[0].select_level(results_data), self._level_pair[1].select_level(results_data)
        return results_data.get_rtp(self._operator, qn_pair)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_str(self._nuclide),
            "rtp",
            self._operator,
            self._level_pair[0].descriptor_str,
            self._level_pair[1].descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        if self._operator in {"M1", "M1lp","M1ln","M1sp","M1sn","M1l0","M1l1","M1s0","M1s1"}:
            observable_text = r"M1_{{{}}}".format(self._operator[2:])
        elif self._operator in {"E2p","E2n","E20","E21","E2"}:
            observable_text = r"E2_{{{}}}".format(self._operator[2:])
        level_pair_text = self._level_pair[0].label_text, self._level_pair[1].label_text
        label = r"B({};{}\rightarrow{})".format(observable_text,level_pair_text[1],level_pair_text[0])  # <f|O|i> = i->f
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        if self._operator in {"M1", "M1lp","M1ln","M1sp","M1sn","M1l0","M1l1","M1s0","M1s1"}:
            observable_text = r"B(M1)"
            units_text = r"\mu_N^2"
        elif self._operator in {"E2p","E2n","E20","E21","E2"}:
            observable_text = r"B(E2)"
            units_text = r"e^2\,\mathrm{fm}^{4}"
        return observable_text, units_text

    
