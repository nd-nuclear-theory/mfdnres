""" NCCI observable differencing.

    Mark A. Caprio
    University of Notre Dame

    - 07/07/23 (mac): Created.  Refactor differencing code from natorb_obs.
"""


import numpy as np
import pandas as pd

import mfdnres.data
import mfdnres.observable

################################################################
# Nmax shifting utility
################################################################

def Nmax_shifted(observable_data,shift):
    """ Shift Nmax values in observable_data.

    Data now found "at Nmax" was previously "at Nmax+shift".

    Example:

        last_step_data = Nmax_shifted(observable_data,0)-Nmax_shifted(observable_data,-2)

    Arguments:
        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

        shift (int): displacement to Nmax index

    Returns:
        (pd.DataFrame): shifted data multi-indexed by (Nmax,hw)
    """
    table = observable_data.reset_index()
    table["Nmax"] -= shift
    shifted_observable_data = pd.DataFrame(table).set_index(["Nmax","hw"])
    return shifted_observable_data


################################################################
# deduced observable: DifferenceRatio
################################################################

def Nmax_difference_ratio(observable_data):
    """Extract ratio between successive Nmax steps for observable.

    Provides step factor for geometric convergence (e.g., 0.5 for Zeno-esqe
    convergence).

    Arguments:
        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

    Returns:
        eta_data (pd.DataFrame): Nmax step ratio multi-indexed by (Nmax,hw)

    """

    Nmax_difference_ratio_data = (
        (Nmax_shifted(observable_data,0)-Nmax_shifted(observable_data,-2))
        /
        (Nmax_shifted(observable_data,-2)-Nmax_shifted(observable_data,-4))
    )
    # drop nan values
    Nmax_difference_ratio_data = mfdnres.data.hw_scan_drop_nan(Nmax_difference_ratio_data)
    return Nmax_difference_ratio_data


class DifferenceRatio(mfdnres.observable.Observable):
    """Observable extractor for difference ratio of observable w.r.t. Nmax.

    """

    def __init__(self, observable1, observable_label_delimiters=None):
        """ Initialize with given parameters.

        Arguments:

            observable1 (Observable): observable to be differenced

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the observable
            appearing in the power, e.g., ("[","]")

        """
        super().__init__()
        self._argument = observable1
        self._observable_label_delimiters = observable_label_delimiters
        if len(self.nuclide_set) == 1:
            self._nuclide = self.nuclide_set.pop()

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        observable_value_mesh = self._argument.data(mesh_data, key_descriptor, verbose=verbose)
        return Nmax_difference_ratio(observable_value_mesh)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "diff-ratio"

        # for single-nuclide observable, include nuclide at start of descriptor
        try:
            arithmetic_operation = "-".join([mfdnres.data.nuclide_str(self._nuclide), arithmetic_operation])
        except AttributeError:
            pass

        return r"{}_{}".format(
            arithmetic_operation,
            self._argument.descriptor_str,
        )

    @property
    def nuclide_label_text(self):
        """Formatted LaTeX text representing nuclide.
        """
        return self._argument.nuclide_label_text

    @property
    def nuclide_set(self):
        """Set of nuclides entering calculation of observable.
        """
        return self._argument.nuclide_set
    
    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_label_text = self._argument.observable_label_text
        observable_label_delimiters = self._observable_label_delimiters
        if observable_label_delimiters is None:
            observable_label_delimiters = ("[","]")
        return r"\eta {}{}{}".format(
            observable_label_delimiters[0],
            observable_label_text,
            observable_label_delimiters[1],
        )

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        axis_label_text = self._argument.axis_label_text
        # Note: Divide out unit inside log.
        return "\eta[{}]".format(axis_label_text[0]), None

    
################################################################
# deduced observable: LogDifference
################################################################

# TODO: extend sign handling on differences, so that can gracefully handle
# positive or negative (or both) differences

def Nmax_log_difference(observable_data):
    """Extract log of Nmax step for observable.

    Log is taken of *negative* of step, so value is defined if observable is
    *decreasing* with Nmax.

    Arguments:
        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

    Returns:
        eta_data (pd.DataFrame): Nmax step ratio multi-indexed by (Nmax,hw)

    """

    Nmax_log_diff_data = np.log(
        -(Nmax_shifted(observable_data,0)-Nmax_shifted(observable_data,-2))
    )
    # drop nan values
    Nmax_log_diff_data = mfdnres.data.hw_scan_drop_nan(Nmax_log_diff_data)
    return Nmax_log_diff_data


class LogDifference(mfdnres.observable.Observable):
    """Observable extractor for logarithmic difference of observable w.r.t. Nmax.

    """

    def __init__(self, observable1, observable_label_delimiters=None):
        """ Initialize with given parameters.

        Arguments:

            observable1 (Observable): observable to be differenced

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the observable
            appearing in the power, e.g., ("[","]")

        """
        super().__init__()
        self._argument = observable1
        self._observable_label_delimiters = observable_label_delimiters
        if len(self.nuclide_set) == 1:
            self._nuclide = self.nuclide_set.pop()

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        observable_value_mesh = self._argument.data(mesh_data, key_descriptor, verbose=verbose)
        return Nmax_log_difference(observable_value_mesh)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "log-diff"

        # for single-nuclide observable, include nuclide at start of descriptor
        try:
            arithmetic_operation = "-".join([mfdnres.data.nuclide_str(self._nuclide), arithmetic_operation])
        except AttributeError:
            pass

        return r"{}_{}".format(
            arithmetic_operation,
            self._argument.descriptor_str,
        )

    @property
    def nuclide_label_text(self):
        """Formatted LaTeX text representing nuclide.
        """
        return self._argument.nuclide_label_text

    @property
    def nuclide_set(self):
        """Set of nuclides entering calculation of observable.
        """
        return self._argument.nuclide_set
    
    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_label_text = self._argument.observable_label_text
        observable_label_delimiters = self._observable_label_delimiters
        if observable_label_delimiters is None:
            observable_label_delimiters = ("[","]")
        return r"\log - \Delta {}{}{}".format(
            observable_label_delimiters[0],
            observable_label_text,
            observable_label_delimiters[1],
        )

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        axis_label_text = self._argument.axis_label_text
        # Note: Could divide out unit inside log.
        return "\log - \Delta [{}]".format(axis_label_text[0]), None
