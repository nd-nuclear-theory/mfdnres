"""Intrinsic QxQ_0 observables calculated from lab-frame two-body observables.

    Mark A. Caprio
    University of Notre Dame

    - 07/24/22 (mac): Created, as q_invariant_intr.py.
    - 10/19/22 (mac): Upgrade from observable tuple to observable object.

"""

import numpy as np

import mfdnres
import mfdnres.am
import mfdnres.data
import mfdnres.observable

import ncci.utils

################################################################
# Q-invariant c.m. correction
################################################################

def get_QxQ_0_intr_rme(results_data, observable_operator, qn_pair, default=np.nan):
    """ Extract QxQ_0-type two-body observable expectation value with analytic c.m. correction.

    See "Q-invariant rederivation" pencilwork page 11.
    """

    # parse quantum numbers
    qnf, qni = qn_pair
    Ji, gi, ni = qnf
    Jf, gf, nf = qni
    if (Ji, gi) != (Jf, gf):
        return default
    J, g = Ji, gi

    # retrieve nuclide
    nuclide = results_data.params["nuclide"]
    A = sum(nuclide)
    Z, N = nuclide
    Tz = 1/2*(Z-N)
    
    # retrieve rmes
    #
    # These are computed with the "E2" normalization for the quadrupole operator.
    QxQ_0_e2 = results_data.get_rme("QxQ_0",qn_pair,rank="tb")
    QpxQp_0_e2 = results_data.get_rme("QpxQp_0",qn_pair,rank="tb")
    QnxQn_0_e2 = results_data.get_rme("QnxQn_0",qn_pair,rank="tb")
    DivxDiv_0 = results_data.get_rme("DivxDiv_0",qn_pair,rank="tb")

    # convert to moment normalization for quadrupole operator
    moment_prefactor = 16*np.pi/5
    QxQ_0 = moment_prefactor*QxQ_0_e2
    QpxQp_0 = moment_prefactor*QpxQp_0_e2
    QnxQn_0 = moment_prefactor*QnxQn_0_e2

    # calculate oscillator length
    hw = results_data.params["hw"]
    b = ncci.utils.oscillator_length(hw)
    b_cm = 1/np.sqrt(A)*b

    # carry out cm correction
    #
    # QxQ_0_intr_obs represents intrinsic QxQ_0, QpxQp_0, or QnxQn_0, depending
    # on given observable_operator.
    QxQ_0_cm = 3*np.sqrt(5)*b_cm**4
    xxx0_cm = -np.sqrt(3)/2*b_cm**2
    delta_f_i = 1 if nf==ni else 0
    if observable_operator == "E20":
        QxQ_0_intr_obs = QxQ_0 - delta_f_i * QxQ_0_cm
    elif observable_operator in ["E2p", "E2n"]:
        alpha = +1 if observable_operator=="E2p" else -1
        middle_term = 2*np.sqrt(5)*xxx0_cm*(DivxDiv_0 + delta_f_i*4*Tz**2*xxx0_cm)
        final_term = -1/4*(1-4*alpha*Tz-4*Tz**2)*delta_f_i*QxQ_0_cm
        QxQ_0_intr_obs = QpxQp_0 + middle_term + final_term

    ## print(QxQ_0_intr_obs)
    return QxQ_0_intr_obs
    
def get_QxQ_0_intr_expectation(results_data, observable_operator, qn, default=np.nan):
    """ Extract QxQ_0-type two-body observable expectation value with analytic c.m. correction.

    See "Q-invariant rederivation" pencilwork page 11.
    """

    # parse quantum number
    J, g, n = qn

    rme = get_QxQ_0_intr_rme(results_data, observable_operator, (qn,qn), default=default)
    expectation_value = 1/mfdnres.am.hat(J)*rme

    return expectation_value

def get_Q0_intr(results_data,observable_operator,qn):
    """ Retrieve intrinsic quadrupole moment from QxQ_0 two-body observables.

    See "Q-invariant rederivation" pencilwork page 11.
    """

    QxQ_0_intr_obs = get_QxQ_0_intr_expectation(results_data,observable_operator,qn)
    Q0_obs = np.sqrt(np.sqrt(5)*QxQ_0_intr_obs)

    ## print("  Q0_obs {}".format(Q0_obs))
    return Q0_obs

################################################################
# Deformation deduced from QxQ_0
################################################################

def get_beta_intr(results_data,observable_operator,qn):
    """ Deduce intrinsic beta from intrinsic QxQ_0 and radius observables.
    """
    
    (J, g, n) = qn

    # retrieve nuclide
    nuclide = results_data.params["nuclide"]
    A = sum(nuclide)
    Z, N = nuclide

    # extract intrinsic quadrupole moment
    Q0 = get_Q0_intr(results_data,observable_operator,qn)

    # convert to beta
    if observable_operator == "E20":
        nucleon_number = A
        radius = results_data.get_radius("r",qn)
    elif observable_operator == "E2p":
        nucleon_number = Z
        radius = results_data.get_radius("rp",qn)
    elif observable_operator == "E2n":
        nucleon_number = N
        radius = results_data.get_radius("rn",qn)

    beta = Q0/(np.sqrt(5/np.pi)*nucleon_number*(radius**2))
    return beta

################################################################
# intrinsic Q0 and beta "observables" -- LEGACY tuple observable
################################################################

# Q0-intr (intrinsic quadrupole moment derived from lab-frame QxQ_0 with CM correction)

def expectation_Q0_intr_extractor(nuclide,observable_operator,observable_qn_list):

    def extractor(results_data):
        # extract qn
        ##print("extractor {}".format([nuclide,observable_operator,observable_qn_list]))
        try:
            resolved_qn_list = tuple([mfdnres.data.resolve_qn(results_data, qn) for qn in observable_qn_list])
            qn, = resolved_qn_list
            if qn is None:
                return np.nan
            value = get_Q0_intr(results_data,observable_operator,qn)
        except Exception as err:
            print("expectation_Q0_intr_extractor for {} {} {} encounters exception {}".format(nuclide,observable_operator,observable_qn_list,err))
            
        return value
    
    return extractor

def expectation_Q0_intr_observable_label(nuclide,observable_operator,observable_qn_list):
    observable_str = r"Q_{{0,{}}}".format(observable_operator[2:])
    qn_str = mfdnres.data.resolve_qn_text(observable_qn_list[0])
    label = r"{}({})".format(observable_str,qn_str)
    return label

def expectation_Q0_intr_axis_label(nuclide,observable_operator,observable_qn_list):
    observable_str = r"Q_0"
    units_str = r"e\,\mathrm{fm}^{2}"
    return observable_str, units_str

mfdnres.data.register_observable(
    "Q0-intr",
    mfdnres.data.Observable(expectation_Q0_intr_extractor, expectation_Q0_intr_observable_label, expectation_Q0_intr_axis_label)
)

# beta-intr (deformation derived from lab-frame QxQ_0 with CM correction)

def expectation_beta_intr_extractor(nuclide,observable_operator,observable_qn_list):

    def extractor(results_data):
        ##print("extractor {}".format([nuclide,observable_operator,observable_qn_list]))
        try:
            # extract qn
            resolved_qn_list = tuple([mfdnres.data.resolve_qn(results_data, qn) for qn in observable_qn_list])
            qn, = resolved_qn_list
            if qn is None:
                return np.nan
            value = get_beta_intr(results_data,observable_operator,qn)
        except Exception as err:
            print("expectation_beta_intr_extractor for {} {} {} encounters exception {}".format(nuclide,observable_operator,observable_qn_list,err))
    
        return value
    
    return extractor

def expectation_beta_intr_observable_label(nuclide,observable_operator,observable_qn_list):
    observable_str = r"\beta_{{{}}}".format(observable_operator[2:])
    qn_str = mfdnres.data.resolve_qn_text(observable_qn_list[0])
    label = r"{}({})".format(observable_str,qn_str)
    return label

def expectation_beta_intr_axis_label(nuclide,observable_operator,observable_qn_list):
    observable_str = r"\beta"
    units_str = None
    return observable_str, units_str

mfdnres.data.register_observable(
    "beta-intr",
    mfdnres.data.Observable(expectation_beta_intr_extractor, expectation_beta_intr_observable_label, expectation_beta_intr_axis_label)
)

################################################################
# intrinsic Q0 and beta observables
################################################################

OBSERVABLE_TAG_STR_BY_OPERATOR = {
    "E2p" : r"p",
    "E2n" : r"n",
    "E20" : r"",
}

class IntrinsicQuadrupoleMoment(mfdnres.observable.Observable):
    """ Observable extractor for intrinsic quadrupole moment Q0.

    Derived from lab-frame QxQ_0 with CM correction.

    """

    def __init__(self, nuclide, operator, level):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): identifier for electromagnetic operator ("E2p",
            "E2n", or "E20")

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
        value = get_Q0_intr(results_data, self._operator,qn)
        return value

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            mfdnres.data.nuclide_str(self._nuclide),
            "Q0",
            self._operator,
            self._level.descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_text = r"Q_{{0{}}}".format(OBSERVABLE_TAG_STR_BY_OPERATOR[self._operator])
        level_text = self._level.label_text
        label = r"{}({})".format(observable_text, level_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"Q_0"
        units_text = r"e\,\mathrm{fm}^{2}"
        return observable_text, units_text


class IntrinsicBeta(mfdnres.observable.Observable):
    """ Observable extractor for beta deformation.

    Derived from lab-frame QxQ_0 with CM correction and radius.

    """

    def __init__(self, nuclide, operator, level):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): identifier for electromagnetic operator ("E2p",
            "E2n", or "E20")

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
        value = get_beta_intr(results_data, self._operator, qn)
        return value

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            mfdnres.data.nuclide_str(self._nuclide),
            "beta",
            self._operator,
            self._level.descriptor_str,
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        observable_text = r"\beta_{{{}}}".format(OBSERVABLE_TAG_STR_BY_OPERATOR[self._operator])
        level_text = self._level.label_text
        label = r"{}({})".format(observable_text,level_text)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"\beta"
        units_text = None
        return observable_text, units_text
    
################################################################
# main
################################################################

def main():
    pass

if __name__ == "__main__":
    main()
