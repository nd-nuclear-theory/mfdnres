""" NCCI observable differencing.

    Mark A. Caprio
    University of Notre Dame

    - 07/07/23 (mac): Created.  Refactor differencing code from natorb_obs.
    - 07/08/23 (mac): Add RelativeDifference observable.
    - 03/04/24 (mac): Add ExponentialExtrapolation observable.
    - 03/21/24 (mac): Extract Nmax_shifted() utility to observable.py.
    - 09/14/24 (mac): Add ExcitationEnergyExponentialExtrapolationSeparate observable.
"""


import numpy as np
import pandas as pd

import mfdnres.data
import mfdnres.observable

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
        (mfdnres.observable.Nmax_shifted(observable_data,0)-mfdnres.observable.Nmax_shifted(observable_data,-2))
        /
        (mfdnres.observable.Nmax_shifted(observable_data,-2)-mfdnres.observable.Nmax_shifted(observable_data,-4))
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
# deduced observable: RelativeDifference
################################################################

def Nmax_relative_difference(observable_data):
    """Extract "logarithmic difference", i.e., relative Nmax step, for observable.

    Arguments:
        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

    Returns:
        eta_data (pd.DataFrame): Nmax step ratio multi-indexed by (Nmax,hw)

    """

    Nmax_relative_difference_data = (
        (mfdnres.observable.Nmax_shifted(observable_data,0)-mfdnres.observable.Nmax_shifted(observable_data,-2))
        /
        ((mfdnres.observable.Nmax_shifted(observable_data,0)+mfdnres.observable.Nmax_shifted(observable_data,-2)) / 2)
    )
    # drop nan values
    Nmax_relative_difference_data = mfdnres.data.hw_scan_drop_nan(Nmax_relative_difference_data)
    return Nmax_relative_difference_data


class RelativeDifference(mfdnres.observable.Observable):
    """Observable extractor for relative difference (or "logarithmic difference") of
    observable w.r.t. Nmax.

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
        return Nmax_relative_difference(observable_value_mesh)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "relative-diff"

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
        ## return r"\Delta {}{}{}/{}{}{}".format(
        ##     observable_label_delimiters[0],
        ##     observable_label_text,
        ##     observable_label_delimiters[1],
        ##     observable_label_delimiters[0],
        ##     observable_label_text,
        ##     observable_label_delimiters[1],
        ## )
        return r"\Delta_{{\mathrm{{rel}}}} {}{}{}".format(
            observable_label_delimiters[0],
            observable_label_text,
            observable_label_delimiters[1],
        )

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        axis_label_text = self._argument.axis_label_text
        return "\Delta_{{\mathrm{{rel}}}} [{}]".format(axis_label_text[0]), None


################################################################
# deduced observable: LogOfDifference
################################################################

def Nmax_log_of_difference(observable_data):
    """Extract log of Nmax step for observable.

    Log is taken of absolute value of difference, so value is defined regardless
    of whether observable is increasing or decreasing with Nmax.

    Arguments:
        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

    Returns:
        eta_data (pd.DataFrame): Nmax step ratio multi-indexed by (Nmax,hw)

    """

    # Prevoiusly: Log was taken of *negative* of step, so value is defined if observable is
    # *decreasing* with Nmax.
    ## Nmax_log_of_difference_data = np.log(
    ##     -(mfdnres.observable.Nmax_shifted(observable_data,0)-mfdnres.observable.Nmax_shifted(observable_data,-2))
    ## )

    Nmax_log_of_difference_data = np.log(
        np.fabs(mfdnres.observable.Nmax_shifted(observable_data,0)-mfdnres.observable.Nmax_shifted(observable_data,-2))
    )
    # drop nan values
    Nmax_log_of_difference_data = mfdnres.data.hw_scan_drop_nan(Nmax_log_of_difference_data)
    return Nmax_log_of_difference_data


class LogOfDifference(mfdnres.observable.Observable):
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
        return Nmax_log_of_difference(observable_value_mesh)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "log-of-diff"

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
        return r"\log \vert \Delta {}{}{} \vert".format(
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
        return r"\log \vert \Delta [{}] \vert".format(axis_label_text[0]), None
    
################################################################
# deduced observable: ExponentialExtrapolation
################################################################

def Nmax_extrapolation_exp3(observable_data):
    """Obtain three-point exponential extrapolation for observable.

    Uses "generalized Zeno's paradox" formula...

        X_infinity = X(N) - eta(N)/[eta(N)-1]*Delta(N)

        Delta(N) = X(N) - X(N-2)

        eta(N) = Delta(N)/Delta(N-2)


    Arguments:
        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

    Returns:
        extrapolated_data (pd.DataFrame): extrapolated data multi-indexed by (Nmax,hw)

    """

    Delta_data = mfdnres.observable.Nmax_shifted(observable_data, 0) - mfdnres.observable.Nmax_shifted(observable_data, -2)
    eta_data = Delta_data / mfdnres.observable.Nmax_shifted(Delta_data, -2)
    extrapolated_data = observable_data - eta_data / (eta_data - 1) * Delta_data
    
    # drop nan values
    extrapolated_data = mfdnres.data.hw_scan_drop_nan(extrapolated_data)
    return extrapolated_data


class ExponentialExtrapolation(mfdnres.observable.Observable):
    """Observable extractor for 3-point eponential extrapolation of
    observable w.r.t. Nmax.

    """

    def __init__(self, observable1, observable_label_delimiters=None):
        """ Initialize with given parameters.

        Arguments:

            observable1 (Observable): observable to be extrapolated

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
        return Nmax_extrapolation_exp3(observable_value_mesh)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "extrapolation-exp3"

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
        return r"{}{}{} \mathrm{{exp3}}".format(
            observable_label_delimiters[0],
            observable_label_text,
            observable_label_delimiters[1],
        )

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        axis_label_text = self._argument.axis_label_text
        return axis_label_text


################################################################
# deduced observable: ExcitationEnergyExponentialExtrapolationSeparate
################################################################

class ExcitationEnergyExponentialExtrapolationSeparate(mfdnres.observable.Observable):
    """Observable extractor for 3-point eponential extrapolation of excitation
    energy, done as difference of independently extrapolated absolute energies.

    Results are equivalent to, but labeling different from, those of

        mfdnres.observable.Difference(
            mfdnres.observable_differencing.ExponentialExtrapolation(
                mfdnres.observable.Energy(
                    nuclide,
                    level,  # e.g., mfdnres.level.LevelQN(qn_ex),
                )
            ),
            mfdnres.observable_differencing.ExponentialExtrapolation(
                mfdnres.observable.Energy(
                    nuclide,
                    reference_level,  # e.g., mfdnres.level.LevelQN(qn),
                )
            ),
        ),

    """

    def __init__(self,  nuclide, level, reference_level, *, Nmax_shift=0, label_as_difference=False):
        """ Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            level (LevelSelector): level

            reference_level (LevelSelector): reference ("ground state") level for energy difference

            Nmax_shift (int, optional): shift in Nmax to apply to excited level,
            for cross-parity excitation energies (e.g., Nmax_shift=1 to
            calculated "Nmax+1" excited level energy relative to "Nmax"
            reference level energy, reported at "Nmax")

            label_as_difference (bool, optional): whether or not to show reference level in observable label

        """
        super().__init__()
        self._nuclide = nuclide
        self._level = level
        self._reference_level = reference_level
        self._Nmax_shift = Nmax_shift
        self._label_as_difference = label_as_difference

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        level_energy_value_mesh = mfdnres.observable.Energy(self._nuclide, self._level).data(mesh_data, key_descriptor, verbose=verbose)
        level_energy_value_mesh = mfdnres.observable.Nmax_shifted(level_energy_value_mesh, self._Nmax_shift)
        reference_level_energy_value_mesh = mfdnres.observable.Energy(self._nuclide, self._reference_level).data(mesh_data, key_descriptor, verbose=verbose)
        return Nmax_extrapolation_exp3(level_energy_value_mesh) - Nmax_extrapolation_exp3(reference_level_energy_value_mesh)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "extrapolation-exp3-sep"  # flag that extrapolation is separate for each level's energy

        descriptor_str = "-".join([
            mfdnres.data.nuclide_str(self._nuclide),
            "energy-ex",
            self._level.descriptor_str,
            self._reference_level.descriptor_str,
        ])

        return r"{}_{}".format(
            arithmetic_operation,
            descriptor_str,
        )

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        level_text = self._level.label_text
        reference_level_text = self._reference_level.label_text
        if self._label_as_difference:
            observable_text = r"E"
            label = r"{}({})-{}({})".format(observable_text,level_text,observable_text,reference_level_text)
        else:
            observable_text = r"E_x"
            label = r"{}({})".format(observable_text,level_text)
        observable_text = "{} \mathrm{{exp3'}".format(observable_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        if self._label_as_difference:
            observable_text = r"\Delta E"
        else:
            observable_text = r"E_x"
        units_text = r"\mathrm{MeV}"
        return observable_text, units_text

    
