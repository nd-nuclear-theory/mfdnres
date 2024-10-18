""" NCCI observable extraction object.

    Mark A. Caprio
    University of Notre Dame

    - 07/26/22 (mac): Created, from stub code in data.py.
    - 10/16/22 (mac): Provide Pow and E2 dimensionless ratio observables.
    - 03/21/23 (mac): Provide secondary axis labels for E2 dimensionless ratio observables.
    - 04/25/23 (mac): Provide ME observable for scalar operators.
    - 08/15/23 (mac): Simplify secondary axis labels for dimensionless ratio observables to omit power.
    - 09/07/23 (mac): Provide conversion dictionaries for matter/proton/neutron observable tag.
    - 09/16/23 (mac): Refactor dimensionless ratio observables to observable_ratio.py.
    - 01/16/24 (zz): Add support for other scalars with user defined label text in ME.
    - 03/21/24 (mac): Add option Nmax_shift to ExcitationEnergy for cross-parity energy differencing.
    - 05/09/24 (mac/zz): Return np.nan for observable value if level is missing (instead of crashing).
    - 08/23/24 (mac): Add wrapper observable OverrideLabels.
    - 09/27/24 (mac): Extend Nmax_shifted() to work with generic index sets.
"""


import numpy as np
import pandas as pd

from . import (
    am,
    analysis,
    data,
)

################################################################
# conversion dictionaries for observable tag ("m"/"p"/"n")
################################################################

NUCLEON_NUMBER_STR_BY_OBSERVABLE_TAG = {
    "m": "A",
    "p": "Z",
    "n": "N",
}
E2_OPERATOR_BY_OBSERVABLE_TAG = {
    "m": "E20",
    "p": "E2p",
    "n": "E2n",
}
RADIUS_OPERATOR_BY_OBSERVABLE_TAG = {
    "m": "r",
    "p": "rp",
    "n": "rn",
}

def nucleon_number_by_observable_tag(nuclide):
    """ Provide dictionary giving nucleon number for mass/proton/neutron observable.

    Arguments:

        nuclide (tuple): (Z, N)

    Returns:

        (dict): observable tag ("m", "p", "n") -> (A, Z, N)
    """

    A = sum(nuclide)
    Z, N = nuclide

    nucleon_number_by_observable_tag = {
        "m": A,
        "p": Z,
        "n": N,
    }
    return nucleon_number_by_observable_tag


################################################################
# Nmax shifting utility
################################################################

def Nmax_shifted(observable_data, shift, *, shift_index_name="Nmax"):
    """Shift Nmax values in observable_data.

    Data now found "at Nmax" was previously "at Nmax+shift".

    Example:

        Delta_data = Nmax_shifted(observable_data, 0) - Nmax_shifted(observable_data, -2)

    Arguments:
        observable_data (pd.DataFrame): data multi-indexed by (Nmax,hw)

        shift (int): displacement to Nmax index

        shift_index_name (str, optional): name of multiindex component to shift
        (defaults to "Nmax")

    Returns:
        (pd.DataFrame): shifted data multi-indexed by (Nmax,hw)

    """
    index_names = observable_data.index.names
    table = observable_data.reset_index()
    table[shift_index_name] -= shift
    shifted_observable_data = pd.DataFrame(table).set_index(index_names)
    return shifted_observable_data


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

        nuclide_set (set): Set of nuclides entering calculation of observable

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

        This default data extractor is provided for "simple" observable objects,
        which identify with a single _nuclide and provide a value extractor.  It
        is to be overridden by any child object which does not follow this
        model, e.g., one which carries out an arithmetic operation on the entire
        mesh.

        The verbose argument provides for debugging of missing observable values
        on mesh.  This argument can be passed through, e.g.,
        mfdnres.data.make_hw_scan_data.

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
        #
        # TODO: Replace use of analysis.make_obs_table with stripped-down tabulation code.
        table = analysis.make_obs_table(mesh_data_selected, key_descriptor, self.value, verbose=verbose)

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
    def nuclide_set(self):
        """Set of nuclides entering calculation of observable.

        This information may be used by calling code in determining, e.g., what
        Nmax values to use in plotting.

        To be overridden by child object if object does not identify with a
        single _nuclide.

        """
        try:
            return {self._nuclide}
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
# deduced observable: OverrideLabels
################################################################

class OverrideLabels(Observable):
    """ Override axis label text

    """

    def __init__(self, observable, observable_label_text=None, axis_label_text=None):
        """Initialize with given parameters.

        Arguments:

            observable (Observable): observable providing value

            observable_label_text (str, optional): observable label text

            axis_label_text (tuple, optional): axis label text

        """
        super().__init__()
        self._observable = observable
        self._observable_label_text = observable_label_text
        self._axis_label_text = axis_label_text

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        return self._observable.data(mesh_data, key_descriptor, verbose=verbose)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return self._observable.descriptor_str

    @property
    def nuclide_label_text(self):
        """Formatted LaTeX text representing nuclide.
        """
        return self._observable.nuclide_label_text

    @property
    def nuclide_set(self):
        """Set of nuclides entering calculation of observable.
        """
        return self._observable.nuclide_set

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_label_text = self._observable_label_text
        if observable_label_text is None:
            observable_label_text = self._observable.observable_label_text
        return observable_label_text

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        axis_label_text = self._axis_label_text
        if axis_label_text is None:
            axis_label_text = self._observable.axis_label_text
        return axis_label_text


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

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the first/second observable
            appearing in the ratio, e.g., (("[","]"),("[","]"))

        """
        super().__init__()
        self._arguments = observable1, observable2
        self._observable_label_delimiters = observable_label_delimiters
        if len(self.nuclide_set) == 1:
            self._nuclide = self.nuclide_set.pop()

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        observable_value_meshes = (
            self._arguments[0].data(mesh_data, key_descriptor, verbose=verbose),
            self._arguments[1].data(mesh_data, key_descriptor, verbose=verbose),
            )
        return observable_value_meshes[0] - observable_value_meshes[1]

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "diff"
        # for single-nuclide observable, include nuclide at start of descriptor
        try:
            arithmetic_operation = "-".join([data.nuclide_str(self._nuclide), arithmetic_operation])
        except AttributeError:
            pass
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
    def nuclide_set(self):
        """Set of nuclides entering calculation of observable.
        """
        return self._arguments[0].nuclide_set | self._arguments[1].nuclide_set

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
        axis_label_texts_same = axis_label_texts[0] == axis_label_texts[1]
        if not axis_label_texts_same:
            raise ValueError("Expect observables with identical axis labels in Difference (found {} and {})".format(axis_label_texts[0], axis_label_texts[1]))
        return r"\Delta {}".format(axis_label_texts[0][0]), axis_label_texts[0][1]

################################################################
# deduced observable: Sum
################################################################

class Sum(Observable):
    """ Observable extractor for sum of observables.

    """

    def __init__(self, observable1, observable2, observable_label_delimiters=None):
        """ Initialize with given parameters.

        Arguments:

            observable1, observable2 (Observable): first and second terms

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the first/second observable
            appearing in the ratio, e.g., (("[","]"),("[","]"))

        """
        super().__init__()
        self._arguments = observable1, observable2
        self._observable_label_delimiters = observable_label_delimiters
        if len(self.nuclide_set) == 1:
            self._nuclide = self.nuclide_set.pop()

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        observable_value_meshes = (
            self._arguments[0].data(mesh_data, key_descriptor, verbose=verbose),
            self._arguments[1].data(mesh_data, key_descriptor, verbose=verbose),
            )
        return observable_value_meshes[0] + observable_value_meshes[1]

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "sum"
        # for single-nuclide observable, include nuclide at start of descriptor
        try:
            arithmetic_operation = "-".join([data.nuclide_str(self._nuclide), arithmetic_operation])
        except AttributeError:
            pass
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
    def nuclide_set(self):
        """Set of nuclides entering calculation of observable.
        """
        return self._arguments[0].nuclide_set | self._arguments[1].nuclide_set

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_label_texts = self._arguments[0].observable_label_text, self._arguments[1].observable_label_text
        arithmetic_operation_sign = "+"
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
        axis_label_texts_same = axis_label_texts[0] == axis_label_texts[1]
        if not axis_label_texts_same:
            raise ValueError("Expect observables with identical axis labels in Difference (found {} and {})".format(axis_label_texts[0], axis_label_texts[1]))
        return r"\Sigma {}".format(axis_label_texts[0][0]), axis_label_texts[0][1]


################################################################
# deduced observable: Ratio
################################################################

class Ratio(Observable):
    """ Observable extractor for ratio of observables.

    """

    def __init__(self, observable1, observable2, observable_label_delimiters=None):
        """Initialize with given parameters.

        Arguments:

            observable1, observable2 (Observable): first and second terms

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the first/second observable
            appearing in the ratio, e.g., (("[","]"),("[","]"))

        """
        super().__init__()
        self._arguments = observable1, observable2
        self._observable_label_delimiters = observable_label_delimiters
        if len(self.nuclide_set) == 1:
            self._nuclide = self.nuclide_set.pop()

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        observable_value_meshes = (
            self._arguments[0].data(mesh_data, key_descriptor, verbose=verbose),
            self._arguments[1].data(mesh_data, key_descriptor, verbose=verbose),
            )
        return observable_value_meshes[0] / observable_value_meshes[1]

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "ratio"

        # for single-nuclide observable, include nuclide at start of descriptor
        try:
            arithmetic_operation = "-".join([data.nuclide_str(self._nuclide), arithmetic_operation])
        except AttributeError:
            pass

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
    def nuclide_set(self):
        """Set of nuclides entering calculation of observable.
        """
        return self._arguments[0].nuclide_set | self._arguments[1].nuclide_set

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
        enforce_axis_label_match = False
        if (enforce_axis_label_match):
            axis_label_texts = self._arguments[0].axis_label_text, self._arguments[1].axis_label_text
            axis_label_texts_same = axis_label_texts[0] == axis_label_texts[1]
            if not axis_label_texts_same:
                raise ValueError("Expect observables with identical axis labels in Ratio (found {} and {})".format(axis_label_texts[0], axis_label_texts[1]))
        return r"\mathrm{Ratio}", None


################################################################
# deduced observable: Power
################################################################

class Power(Observable):
    """ Observable extractor for power of observable.

    TODO: In observable label, provide some sort of pass-through mechanism to
    obtain r^4() instead of r()^4.

    """

    def __init__(self, observable1, power, observable_label_delimiters=None):
        """ Initialize with given parameters.

        Arguments:

            observable1 (Observable): observable to be exponentiated

            power (int): power to which to raise observable

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the observable
            appearing in the power, e.g., ("[","]")

        """
        super().__init__()
        self._argument = observable1
        self._power = power
        self._observable_label_delimiters = observable_label_delimiters
        if len(self.nuclide_set) == 1:
            self._nuclide = self.nuclide_set.pop()

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        observable_value_mesh = self._argument.data(mesh_data, key_descriptor, verbose=verbose)
        return observable_value_mesh**self._power

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        arithmetic_operation = "pow"

        # for single-nuclide observable, include nuclide at start of descriptor
        try:
            arithmetic_operation = "-".join([data.nuclide_str(self._nuclide), arithmetic_operation])
        except AttributeError:
            pass

        return r"{}_{:g}_{}".format(
            arithmetic_operation,
            self._power,
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
        try:
            # use "nice" power label if provided by observable, e.g., "r^4(...)"
            # rather than generic "r(...)^4"
            observable_label_text = self._argument.observable_power_label_text(self._power)
            return observable_label_text
        except AttributeError:
            # fall back on generic power notation
            observable_label_text = self._argument.observable_label_text
            observable_label_delimiters = self._observable_label_delimiters
            if observable_label_delimiters is None:
                observable_label_delimiters = ("","")
            return r"{}{}{}^{:d}".format(
                observable_label_delimiters[0],
                observable_label_text,
                observable_label_delimiters[1],
                self._power,
            )

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        axis_label_text = self._argument.axis_label_text
        if self._power==1/2:
            power_text = r"1/2"
        else:
            power_text = "{:d}".format(self._power)
        return "{}^{{{:s}}}".format(axis_label_text[0], power_text), "{}^{:s}".format(axis_label_text[1], power_text)


################################################################
# deduced observable: FixSignTo
################################################################

class FixSignTo(Observable):
    """Fix sign of one observable by multiplying it by sign of another observable
    (e.g., to fix signs from arbitrary phases of initial and final states).

    As a special case, takes the absolute value for a single observable.

    """

    def __init__(self, observable1, observable2=None, observable_label_delimiters=None):
        """Initialize with given parameters.

        Arguments:

            observable1 (Observable): observable providing value

            observable2 (Observable, optional): observable providing sign; if
                omitted, observable1 is used, which amounts to taking
                abs(observble1)

            observable_label_delimiters (tuple, optional): left/right delimiter
                pair to put around the label for the observable,
                e.g., (r"\vert",r"\vert")

        """
        super().__init__()
        if observable2 is None:
            observable2 = observable1
        self._arguments = observable1, observable2
        self._observable_label_delimiters = observable_label_delimiters

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        observable_value_meshes = (
            self._arguments[0].data(mesh_data, key_descriptor, verbose=verbose),
            self._arguments[1].data(mesh_data, key_descriptor, verbose=verbose),
            )
        return observable_value_meshes[0] * observable_value_meshes[1].apply(np.sign, raw=True)

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
    def nuclide_set(self):
        """Set of nuclides entering calculation of observable.
        """
        return self._arguments[0].nuclide_set | self._arguments[1].nuclide_set

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_label_delimiters = self._observable_label_delimiters
        if observable_label_delimiters is None:
            observable_label_delimiters = ("","")
        return r"{}{}{}".format(
            observable_label_delimiters[0],
            self._arguments[0].observable_label_text,
            observable_label_delimiters[1],
        )

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
        if qn is None:
            return np.nan
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
# observable: ExcitationEnergy
################################################################

class ExcitationEnergy(Observable):
    """Observable extractor for difference of energies in single nucleus, labeled as "excitation energy".

    Note that a difference of energies may alternatively be obtained as a
    Difference of two Energy observables.  But then the descriptor text (used in
    data filenames) is that for a generic difference, as is automatic labeling
    (which will appear as delta E), while the labeling for ExcitationEnergy is
    likely more appropriate for presentation of the excitation energy.

    For an energy difference between levels of different parity, the Nmax_shift
    option may be used.

    """

    def __init__(self, nuclide, level, reference_level, *, Nmax_shift=0, label_as_difference=False):
        """Initialize with given parameters.

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
        ## super().__init__(Energy(nuclide, level), Energy(nuclide, reference_level))  # if implement as child of Difference rather than freshly calculating value()
        self._nuclide = nuclide
        self._level = level
        self._reference_level = reference_level
        self._Nmax_shift = Nmax_shift
        self._label_as_difference = label_as_difference

    ## def value(self, results_data):
    ##     """ Extract observable.
    ##     """
    ##     qn = self._level.select_level(results_data)
    ##     reference_qn = self._reference_level.select_level(results_data)
    ##     return results_data.get_energy(qn) - results_data.get_energy(reference_qn)

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """
        level_energy_value_mesh = Energy(self._nuclide, self._level).data(mesh_data, key_descriptor, verbose=verbose)
        print("E unshifted {}".format(level_energy_value_mesh))
        level_energy_value_mesh = Nmax_shifted(level_energy_value_mesh, self._Nmax_shift)
        print("E shifted {}".format(level_energy_value_mesh))
        reference_level_energy_value_mesh = Energy(self._nuclide, self._reference_level).data(mesh_data, key_descriptor, verbose=verbose)
        print("Eref {}".format(reference_level_energy_value_mesh))
        return level_energy_value_mesh - reference_level_energy_value_mesh
    
    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_str(self._nuclide),
            "energy-ex",
            self._level.descriptor_str,
            self._reference_level.descriptor_str,
        ])

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


################################################################
# observable: Isospin
################################################################

class Isospin(Observable):
    """Observable extractor for level isospin.

    This is the isosopin as calculated and then printed to low precision in the
    mfdn results.

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
        if qn is None:
            return np.nan
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
        if qn is None:
            return np.nan
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
        if qn is None:
            return np.nan
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

    def observable_power_label_text(self, power):
        """ Formatted LaTeX text representing observable raised to power.

        Arguments:

            power (int): exponent
        """
        observable_text = RADIUS_STR_BY_OPERATOR[self._operator]
        level_text = self._level.label_text
        label = r"{}^{:d}({})".format(observable_text,power,level_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"r"
        units_text = r"\mathrm{fm}"
        return observable_text, units_text


################################################################
# observable: RadiusPower -- OMIT in favor of deduced Power
################################################################

## class RadiusPower(Observable):
##     """ Observable extractor for power of r.m.s. radius.
##
##     """
##
##     def __init__(self, nuclide, operator, level, power):
##         """Initialize with given parameters.
##
##         Arguments:
##
##             nuclide (tuple): (Z, N)
##
##             operator (str): identifier for radius operator, as accepted by
##             MFDnResultsData.get_radius()
##
##             level (LevelSelector): level
##
##             power (int): power to which to raise radius
##
##         """
##         super().__init__()
##         self._nuclide = nuclide
##         self._operator = operator
##         self._level = level
##         self._power = power
##
##     def value(self, results_data):
##         """ Extract observable.
##         """
##         qn = self._level.select_level(results_data)
##         return results_data.get_radius(self._operator, qn)**self._power
##
##     @property
##     def descriptor_str(self):
##         """ Text string describing observable.
##         """
##         return "-".join([
##             data.nuclide_str(self._nuclide),
##             "radius-pow{:d}".format(self._power),
##             self._operator,
##             self._level.descriptor_str,
##         ])
##
##     @property
##     def observable_label_text(self):
##         """ Formatted LaTeX text representing observable.
##         """
##         observable_text = RADIUS_STR_BY_OPERATOR[self._operator]
##         level_text = self._level.label_text
##         label = r"{}^{:d}({})".format(observable_text,self._power,level_text)
##         return label
##
##     @property
##     def axis_label_text(self):
##         """ Formatted LaTeX text representing axis label.
##         """
##         observable_text = r"r^{:d}".format(self._power)
##         units_text = r"\mathrm{{fm}}^{:d}".format(self._power)
##         return observable_text, units_text


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
        if qn is None:
            return np.nan
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

    def observable_power_label_text(self, power):
        """ Formatted LaTeX text representing observable raised to given power.

        Arguments:

            power (int): exponent
        """
        if self._operator in {"M1", "M1lp","M1ln","M1sp","M1sn","M1l0","M1l1","M1s0","M1s1"}:
            observable_text = r"\mu_{{{}}}".format(self._operator[2:])
        elif self._operator in {"E2p","E2n","E20","E21","E2"}:
            observable_text = r"Q_{{{}}}".format(self._operator[2:])
        level_text = self._level.label_text
        if power==1/2:
            power_text = "1/2"
        else:
            power_text = "{:d}".format(power)
        label = r"{}^{{{:s}}}({})".format(observable_text,power_text,level_text)
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
# observable: RME
################################################################

class RME(Observable):
    """ Observable extractor for electromagnetic reduced matrix element.

    """

    # TODO restore E1 functionality; consider combining E0/E1/E2 implementation

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
        if (qn_pair[0] is None) or (qn_pair[1] is None):
            return np.nan
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
            observable_text = r"M1_{{{}}}".format(self._operator[2:])
        elif self._operator in {"E2p","E2n","E20","E21","E2"}:
            ## observable_text = r"Q_{{2{}}}".format(self._operator[2:])
            observable_text = r"E2_{{{}}}".format(self._operator[2:])
        elif self._operator in {"E0p","E0n","E00","E01","E0"}:
            ## observable_text = r"Q_{{2{}}}".format(self._operator[2:])
            observable_text = r"E0_{{{}}}".format(self._operator[2:])
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
        elif self._operator in {"E0p","E0n","E00","E01","E0"}:
            ## observable_text = r"\langle Q_2 \rangle"
            observable_text = r"\langle \mathcal{{M}}(E0) \rangle"
            units_text = r"e\,\mathrm{fm}^{2}"
        return observable_text, units_text


################################################################
# observable: ME -- scalar observable only
################################################################

class ME(Observable):
    """ Observable extractor for electromagnetic non-reduced matrix element.

    Applies only to scalar (E0) observable by default.

    Also provides support for other scalars with user defined observable and axis label text.

    """

    def __init__(self, nuclide, operator, levelf, leveli, *, label_text = None, rank="ob"):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): identifier for electromagnetic operator, as accepted
            by MFDnResultsData.get_rme()

            level (LevelSelector): level

            label_text (str, (str, str)): overide for observable and axis label text
            (observable_label_text, axis_label_text)

        """
        super().__init__()
        self._nuclide = nuclide
        self._operator = operator
        self._level_pair = levelf, leveli
        self._label_text = label_text
        self._rank = rank

    def value(self, results_data):
        """ Extract observable.
        """
        qn_pair = self._level_pair[0].select_level(results_data), self._level_pair[1].select_level(results_data)
        if (qn_pair[0] is None) or (qn_pair[1] is None):
            return np.nan
        # assert(self._operator in {"E0p","E0n","E00","E01","E0"})
        if self._label_text is None:
            assert(self._operator in {"E0p","E0n","E00","E01","E0"})
        return results_data.get_me(self._operator, qn_pair, rank = self._rank)

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            data.nuclide_str(self._nuclide),
            "me",
            self._operator,
            self._level_pair[0].descriptor_str,
            self._level_pair[1].descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        if self._label_text is not None:
            observable_text = self._label_text[0]
        elif self._operator in {"E0p","E0n","E00","E01","E0"}:
            ## observable_text = r"Q_{{2{}}}".format(self._operator[2:])
            observable_text = r"E0_{{{}}}".format(self._operator[2:])
            observable_text = r"\mathcal{{M}}({})".format(observable_text)
        level_pair_text = self._level_pair[0].label_text, self._level_pair[1].label_text
        label = r"\langle {} \vert {} \vert {} \rangle".format(level_pair_text[0],observable_text,level_pair_text[1])  # <f|O|i> = i->f
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        if self._label_text is not None:
            observable_text = self._label_text[1][0]
            units_text = self._label_text[1][1]
        elif self._operator in {"E0p","E0n","E00","E01","E0"}:
            ## observable_text = r"\langle Q_2 \rangle"
            observable_text = r"\langle \mathcal{{M}}(E0) \rangle"
            units_text = r"e\,\mathrm{fm}^{2}"
        return observable_text, units_text

################################################################
# observable: RTP
################################################################

class RTP(Observable):
    """Observable extractor for electromagnetic reduced transition probability.

    """

    # TODO restore E1 functionality; consider combining E0/E1/E2 implementation

    def __init__(self, nuclide, operator, levelf, leveli):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): identifier for electromagnetic operator, as accepted
            by MFDnResultsData.get_rme()

            levelf (LevelSelector): final level

            leveli (LevelSelector): final initial level

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
        if (qn_pair[0] is None) or (qn_pair[1] is None):
            return np.nan
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
        elif self._operator in {"E0p","E0n","E00","E01","E0"}:
            observable_text = r"E0_{{{}}}".format(self._operator[2:])
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
        elif self._operator in {"E0p","E0n","E00","E01","E0"}:
            observable_text = r"B(E0)"
            units_text = r"e^2\,\mathrm{fm}^{4}"
        return observable_text, units_text
