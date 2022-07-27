""" NCCI observable extraction object.

    Mark A. Caprio
    University of Notre Dame

    - 07/26/22 (mac): Created, from stub code in data.py.
"""


import numpy as np

from . import (
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

        observable (ResultsData -> float OR np.array): Retrieve observable (or
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

    def observable(self, results_data):
        """ Extract observable.
        """
        return np.nan

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
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
# observable: energy
################################################################

class Energy(Observable):
    """ Provides observable extractor for level energy.
    """

    def __init__(self, nuclide, level):
        """ Initialize with given parameters.

        Arguments:
            level (LevelSelector): Level
        """
        super().__init__()
        self._nuclide = nuclide
        self._level = level

    def observable(self, results_data):
        """ Extract observable.
        """
        qn = self._level.select_level(results_data)
        return results_data.get_energy(qn)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_descriptor_str(self._nuclide),
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
