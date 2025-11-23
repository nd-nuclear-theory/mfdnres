"""Observables for two-state mixing, as diagnosed by the transition matrix element.

    For mixing angle relations, see, e.g., McCoy et al., PLB, DOI:
    10.1016/j.physletb.2024.138870.

    Mark A. Caprio
    University of Notre Dame

    - 10/18/23 (mac): Created, extracting code from 10be-shape_obs.py.
    - 07/21/24 (mac): Add extraction of mixing angle from fragmentation.

"""

import numpy as np

import mfdnres.data
import mfdnres.level
import mfdnres.observable


################################################################
# mixing angle calculation
################################################################

def mixing_angle_from_trans(results_data, operator, qn_pair):
    """Mixing angle extracted from transition between states.

    Assumes vanishing E0 (or E2) between the unmixed states, which yields

        tan(2*theta) = 2 * M'_12 / (M'_22 - M'_11)

    Arguments:

        results_data (MFDnResultsData): results data

        operator (str): operator type ("E0p", "E0n", "E00", "E2p", "E2n", "E20")

        qn_pair (tuple): pair (qn_1, qn_2) of (J,g,n) quantum numbers for states

    Returns:

        theta (float): magnitue of mixing angle in degrees

    """

    qn_1, qn_2 = qn_pair
    M11 = results_data.get_rme(operator, (qn_1, qn_1))
    M12 = results_data.get_rme(operator, (qn_1, qn_2))
    M22 = results_data.get_rme(operator, (qn_2, qn_2))

    theta = np.abs(1/2*np.arctan(2*M12/(M22-M11)))
    return theta


def mixing_angle_from_fragmentation(results_data, operator, qn_pair, qn_i):
    """Mixing angle extracted from fragmentation of transition over states.

    Assumes vanishing transition to the unmixed second state, which yields

        tan(theta) = abs(M'_2 / M'_1)

    Arguments:

        results_data (MFDnResultsData): results data

        operator (str): operator type ("E0p", "E0n", "E00", "E2p", "E2n", "E20")

        qn_pair (tuple): pair (qn_1, qn_2) of (J,g,n) quantum numbers for states

        qn_i (tuple): (J,g,n) for initial state

    Returns:

        theta (float): magnitue of mixing angle in degrees

    """

    qn_1, qn_2 = qn_pair
    M1 = results_data.get_rme(operator, (qn_1, qn_i))
    M2 = results_data.get_rme(operator, (qn_2, qn_i))

    theta = np.abs(np.arctan(M2/M1))
    return theta


def mixing_angle(results_data, operator, qn_pair, qn_i):
    """Mixing angle dispatch function.

    Returns mixing angle either from transition between states or from
    fragmentation of transition to states, according to value of argument qn_i.

    Arguments:

        results_data (MFDnResultsData): results data

        operator (str): operator type ("E0p", "E0n", "E00", "E2p", "E2n", "E20")

        qn_pair (tuple): pair (qn_1, qn_2) of (J,g,n) quantum numbers for states

        qn_i (tuple, optional): (J,g,n) for initial state

    Returns:

        theta (float): magnitue of mixing angle in degree

    """
    if qn_i is None:
        theta = mixing_angle_from_trans(results_data, operator, qn_pair)
    else:
        theta = mixing_angle_from_fragmentation(results_data, operator, qn_pair, qn_i)
    return theta

################################################################
# mixing observables
################################################################


class MixingObservable(mfdnres.observable.Observable):

    """Observable extractor interface class.

    Provides common __init__ method for mixing problem.

    """

    def __init__(self, nuclide, operator, subspace, qn_i=None):
        """ Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): operator type ("E0p", "E0n", "E00")

            subspace (tuple): (J,g) quantum numbers

            qn_i (tuple, optional): (J,g,n) for initial state (for branching analysis)

        """
        super().__init__()

        if operator not in {"E0p", "E0n", "E00", "E2p", "E2n", "E20"}:
            raise ValueError("inappropriate operator code for mixing calculation")
        
        self._nuclide = nuclide
        self._operator = operator
        self._operator_family = operator[:-1]  # "E0" or "E2"
        self._operator_flavor = operator[-1]  # "p", "n", or "0"
        self._J, self._g = subspace
        self._qn_pair = ((*subspace, 1), (*subspace, 2))
        self._qn_i = qn_i
       
        
class MixingAngle(MixingObservable):
    """Observable extractor for magnitude of mixing angle (converted to deg).

    """

    def value(self, results_data):
        """ Extract observable.
        """
        theta = mixing_angle(results_data, self._operator, self._qn_pair, self._qn_i)
        theta_deg = theta*180/np.pi
        return theta_deg

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            "mixing-angle",
            self._operator,
            "J{:1.0f}".format(self._J),
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        qn_text = mfdnres.data.qn_text((self._J,self._g,None), show_parity=True, show_index=False)
        label = r"\theta_{{{}}}~[\mathrm{{from}}~{}_{}]".format(qn_text, self._operator_family, self._operator_flavor)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"\theta"
        units_text = r"\mathrm{deg}"
        return observable_text, units_text

    
class MixingAdmixture(MixingObservable):
    """Observable extractor for admixture norm (sin^2 theta).

    """

    def value(self, results_data):
        """ Extract observable.
        """
        theta = mixing_angle(results_data, self._operator, self._qn_pair, self._qn_i)
        admixture = np.sin(theta)**2
        return admixture

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            "mixing-admixture",
            self._operator,
            "J{:1.0f}".format(self._J),
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        ## label = r"\sin^2\theta({:1.0f}^+_{{0\hbar\omega}}/{:1.0f}^+_{{2\hbar\omega}})_{}".format(self._J, self._J, self._operator_flavor)
        ##label = r"\sin^2\theta_{{{:1.0f},{}}}".format(self._J, self._operator_flavor)
        qn_text = mfdnres.data.qn_text((self._J,self._g,None), show_parity=True, show_index=False)
        label = r"\sin^2\theta_{{{}}}~[\mathrm{{from}}~{}_{}]".format(qn_text, self._operator_family, self._operator_flavor)        
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"\sin^2\theta"
        units_text = None
        return observable_text, units_text

    
class MixingMatrixElement(MixingObservable):
    """ Observable extractor for magnitude of mixing matrix element.

    """

    def value(self, results_data):
        """ Extract observable.
        """
        theta = mixing_angle(results_data, self._operator, self._qn_pair, self._qn_i)
        r = 1/2*np.tan(2*theta)
        qn_1, qn_2 = self._qn_pair
        energy_difference = results_data.get_energy(qn_2) - results_data.get_energy(qn_1)
        V = r*(1+4*r**2)**(-1/2)*energy_difference
        return V

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            "mixing-matrix-element",
            self._operator,
            "J{:1.0f}".format(self._J),
        ])

    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        ## label = r"\langle {:1.0f}^+_{{0\hbar\omega}} \vert V \vert {:1.0f}^+_{{2\hbar\omega}} \rangle_{}".format(self._J, self._J, self._operator_flavor)
        ## label = r"V_{{{:1.0f},{}}}".format(self._J, self._operator_flavor)
        ## P_str = "+" if self._g==0 else "-"
        qn_text = mfdnres.data.qn_text((self._J,self._g,None), show_parity=True, show_index=False)

        label = r"V_{{{}}}~[\mathrm{{from}}~{}_{}]".format(qn_text, self._operator_family, self._operator_flavor)
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"\langle V \rangle"
        units_text = r"\mathrm{MeV}"
        return observable_text, units_text
    

class MixingEnergyDifference(mfdnres.observable.ExcitationEnergy):
    """Observable extractor for mixing energy difference.

    """

    def __init__(self, nuclide, subspace):
        """ Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            subspace (tuple): (J,g) quantum numbers

        """
        self._nuclide = nuclide
        self._J, self._g = subspace
        self._qn_pair = ((*subspace, 1), (*subspace, 2))
        super().__init__(
            nuclide,
            mfdnres.level.LevelQN(self._qn_pair[1]), mfdnres.level.LevelQN(self._qn_pair[0]),
            label_as_difference=True,
        )
        
    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            "mixing-energy-difference",
            "J{:1.0f}".format(self._J),
        ])
    



################################################################
# main
################################################################

def main():
    pass

if __name__ == "__main__":
    main()
