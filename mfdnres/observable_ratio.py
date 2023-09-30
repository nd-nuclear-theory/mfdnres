"""NCCI dimensionless ratio observables.

    Mark A. Caprio
    University of Notre Dame

    - 09/07/23 (mac): Created.  Refactor dimensionless ratio observables from
        observable.py.

"""

import numpy as np

import mfdnres.data
import mfdnres.observable

################################################################
# deduced observable: RatioBE2Q2
################################################################
        
class RatioBE2Q2(mfdnres.observable.Ratio):
    """Observable extractor for dimensionless ratio B(E2)/(e^2Q^2).

    """

    def __init__(self, observable1, observable2, observable_label_delimiters=(("",""),("[e^2","]"))):
        """Initialize with given parameters.

        Arguments:

            observable1, observable2 (Observable): first and second terms

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the first/second observable
            appearing in the ratio, e.g., (("[","]"),("[","]"))

        """
        if not (isinstance(observable1, mfdnres.observable.RTP) and isinstance(observable2, mfdnres.observable.Moment)):
            raise ValueError("Unexpected observable types in RatioBE2Q2 (found {} and {})".format(observable1, observable2))
        super().__init__(observable1, mfdnres.observable.Power(observable2, 2), observable_label_delimiters)

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        ##return r"B(E2)/(e^2Q^2)", None
        return r"B(E2)/(eQ)^2", None

    @property
    def secondary_axis_label_text(self):
        """ Formatted LaTeX text representing axis label for secondary "calibrated" axis.
        """

        calibrator_label = "Q"  # redefined form "[via Q]"
        return r"{}~({})~~[\mathrm{{via}}~{}]".format(*self._arguments[0].axis_label_text,calibrator_label)

    
################################################################
# deduced observable: RatioBE2r4
################################################################

# TODO enable analogous E0 ratio support

class RatioBE2r4(mfdnres.observable.Ratio):
    """Observable extractor for dimensionless ratio B(E2)/(e^2r^4).

    """

    def __init__(self, observable1, observable2, observable_label_delimiters=(("",""),("[e^2","]"))):
        """Initialize with given parameters.

        Arguments:

            observable1, observable2 (Observable): first and second terms

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the first/second observable
            appearing in the ratio, e.g., (("[","]"),("[","]"))

        """
        if not (isinstance(observable1, mfdnres.observable.RTP) and isinstance(observable2, mfdnres.observable.Radius)):
            raise ValueError("Unexpected observable types in RatioBE2r4 (found {} and {})".format(observable1, observable2))
        super().__init__(observable1, mfdnres.observable.Power(observable2, 4), observable_label_delimiters)

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        ## axis_label_texts = self._arguments[0].axis_label_text, self._arguments[1].axis_label_text
        # e.g., ('B(E2)', 'e^2\\,\\mathrm{fm}^{4}') and ('r^4', '\\mathrm{fm}^4'))
        
        # Would need more info than available from observable objects to
        # specialize as, e.g., "$B(E2)/(e^2r_p^4)$".
        
        return r"B(E2)/(e^2r^4)", None

    @property
    def secondary_axis_label_text(self):
        """ Formatted LaTeX text representing axis label for secondary "calibrated" axis.
        """

        calibrator_label = "r"  # redefined form "[via r]"
        return r"{}~({})~~[\mathrm{{via}}~{}]".format(*self._arguments[0].axis_label_text,calibrator_label)

################################################################
# deduced observable: RatioQr2
################################################################

class RatioQr2(mfdnres.observable.Ratio):
    """ Observable extractor for dimensionless ratio Q/r^2.

    """

    def __init__(self, observable1, observable2, observable_label_delimiters=None):
        """Initialize with given parameters.

        Arguments:

            observable1, observable2 (Observable): first and second terms

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the first/second observable
            appearing in the ratio, e.g., (("[","]"),("[","]"))

        """
        if not (
                isinstance(observable1, mfdnres.observable.Moment)
                and isinstance(observable2, mfdnres.observable.Radius)
        ):
            raise ValueError("Unexpected observable types in RatioQr2 (found {} and {})".format(observable1, observable2))
        super().__init__(observable1, mfdnres.observable.mfdnres.observable.Power(observable2, 2), observable_label_delimiters)

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
       
        return r"Q/r^2", None

    @property
    def secondary_axis_label_text(self):
        """ Formatted LaTeX text representing axis label for secondary "calibrated" axis.
        """
        
        ## calibrator_label = *self._arguments[1].axis_label_text  # original form "[via r^2]"
        calibrator_label = "r"  # redefined form "[via r]"
        return r"{}~({})~~[\mathrm{{via}}~{}]".format(*self._arguments[0].axis_label_text,calibrator_label)

    
################################################################
# deduced observable: RatiorQ12
################################################################

class RatiorQ12(mfdnres.observable.Ratio):
    """ Observable extractor for dimensionless ratio r/[Q^(1/2)].

    """

    def __init__(self, observable1, observable2, observable_label_delimiters=None):
        """Initialize with given parameters.

        Caution: Arguments are still given in the order (Q, r), even with the
        ratio "flipped" so r is in the numerator.

        Arguments:

            observable1, observable2 (Observable): first and second terms

            observable_label_delimiters (tuple, optional): left/right delimiter
            pairs to put around the labels for the first/second observable
            appearing in the ratio, e.g., (("[","]"),("[","]"))

        """
        if not (isinstance(observable1, mfdnres.observable.Moment) and isinstance(observable2, mfdnres.observable.Radius)):
            raise ValueError("Unexpected observable types in RatiorQ12 (found {} and {})".format(observable1, observable2))
        super().__init__(observable2, mfdnres.observable.Power(observable1, 1/2), observable_label_delimiters)

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
       
        return r"r/Q^{1/2}", None


    @property
    def secondary_axis_label_text(self):
        """ Formatted LaTeX text representing axis label for secondary "calibrated" axis.
        """
        
        ## calibrator_label = *self._arguments[1].axis_label_text  # original form "[via r^2]"
        calibrator_label = "Q"
        print(r"{}~({})~~[\mathrm{{via}}~{}]".format(*self._arguments[0].axis_label_text,calibrator_label))
        return r"{}~({})~~[\mathrm{{via}}~{}]".format(*self._arguments[0].axis_label_text,calibrator_label)

    
################################################################
# deduced observable: BetaFromRatioQr2
################################################################

class BetaFromRatioQr2(mfdnres.observable.Observable):
    """ Observable extractor for beta deformation.

    """

    def __init__(self, nuclide, observable_tag, level, J, K, axis_label_has_observable_tag=True):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            observable_tag (str): identifier tag for beta observable ("p", "n", or "m")

            level (LevelSelector): level

            J (float): J quantum number for level (assumed unique across mesh)

            K (float): K quantum number for level (assumed unique across mesh)

            axis_label_has_observable_tag (bool): whether to include subscript
                on axis label (as beta_{m,p,n})

        """
        super().__init__()
        self._nuclide = nuclide
        self._observable_tag = observable_tag
        self._level = level
        self._J = J
        self._K = K
        self._axis_label_has_observable_tag = axis_label_has_observable_tag

    def data(self, mesh_data, key_descriptor, verbose=False):
        """ Extract data frame of observable values over mesh.
        """

        nuclide = self._nuclide
        observable_tag = self._observable_tag
        level = self._level
        J = self._J
        K = self._K
        
        # retrieve nucleon number
        A = sum(nuclide)
        Z, N = nuclide
        nucleon_number = mfdnres.observable.nucleon_number_by_observable_tag(nuclide)[observable_tag]


        # select E2 moment
        e2_operator = mfdnres.observable.E2_OPERATOR_BY_OBSERVABLE_TAG[observable_tag]
        moment_observable = mfdnres.observable.Moment(nuclide, e2_operator, level)
        
        # select radius
        radius_operator = mfdnres.observable.RADIUS_OPERATOR_BY_OBSERVABLE_TAG[observable_tag]
        radius_observable = mfdnres.observable.Radius(nuclide, radius_operator, level)

        # deduce ratio
        ratio_observable = RatioQr2(moment_observable, radius_observable)

        # calculate mesh
        ratio_mesh = ratio_observable.data(mesh_data, key_descriptor, verbose=verbose)
        
        # convert to beta
        prefactor = (J+1)*(2*J+3)/(3*K**2-J*(J+1)) * np.sqrt(np.pi/5)/nucleon_number
        beta_mesh = prefactor * ratio_mesh

        print("ratio_mesh {}".format(ratio_mesh))
        print("beta_mesh {}".format(beta_mesh))
        
        return beta_mesh

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            mfdnres.data.nuclide_str(self._nuclide),
            "beta-from-ratio-q-rsqr",
            self._observable_tag,
            self._level.descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_text = r"\beta_{{{}}}".format(self._observable_tag)
        level_text = self._level.label_text
        label = r"{}({})".format(observable_text,level_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        if self._axis_label_has_observable_tag:
            observable_text = r"\beta_{{{}}}".format(self._observable_tag)
        else:
            observable_text = r"\beta"

        units_text = None
        return observable_text, units_text
    

