"""Observables for two-state mixing, as diagnosed by the transition matrix element.

    Mark A. Caprio
    University of Notre Dame

    - 10/18/23 (mac): Created, extracting code from 10be-shape_obs.py.

"""

import numpy as np

import mfdnres.level
import mfdnres.observable


################################################################
# E0 (or E2) mixing analysis
################################################################

def mixing_angle_from_trans(results_data, operator, qn_pair):
    """Extract magnitude of mixing angle from E0 (or E2) between states.

    Assumes vanishing E0 (or E2) between the unmixed states, which yields

        tan(2*theta) = 2 M'_12 / (M'_22 - M'_11)

    Arguments:

        results_data (MFDnResultsData): results data

        operator (str): operator type ("E0p", "E0n", "E00", "E2p", "E2n", "E20")

        qn_pair (tuple): pair (qn_1, qn_2) of (J,g,n) quantum numbers for states

    Returns:

        theta (float): magnitue of mixing angle in degree

    """

    qn_1, qn_2 = qn_pair
    M11 = results_data.get_rme(operator, (qn_1, qn_1))
    M12 = results_data.get_rme(operator, (qn_1, qn_2))
    M22 = results_data.get_rme(operator, (qn_2, qn_2))

    theta = np.abs(1/2*np.arctan(2*M12/(M22-M11)))
    return theta


class MixingObservable(mfdnres.observable.Observable):
    """Observable extractor interface class.

    Provides common __init__ method for E0 mixing problem.

    """

    def __init__(self, nuclide, operator, subspace):
        """ Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            operator (str): operator type ("E0p", "E0n", "E00")

            subspace (tuple): (J,g) quantum numbers

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
        
        
class MixingAngle(MixingObservable):
    """Observable extractor for magnitude of mixing angle (converted to deg).

    """

    def value(self, results_data):
        """ Extract observable.
        """
        theta = mixing_angle_from_trans(results_data, self._operator, self._qn_pair)
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
        P_str = "+" if self._g==0 else "-"
        label = r"\theta_{{{:1.0f}^{}}}~[\mathrm{{from}}~{}_{}]".format(self._J, P_str, self._operator_family, self._operator_flavor)
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
        theta = mixing_angle_from_trans(results_data, self._operator, self._qn_pair)
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
        P_str = "+" if self._g==0 else "-"
        label = r"\sin^2\theta_{{{:1.0f}^{}}}~[\mathrm{{from}}~E0_{}]".format(self._J, P_str, self._operator_flavor)        
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
        theta = mixing_angle_from_trans(results_data, self._operator, self._qn_pair)
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
        P_str = "+" if self._g==0 else "-"
        label = r"V_{{{:1.0f}^{}}}~[\mathrm{{from}}~E0_{}]".format(self._J, P_str, self._operator_flavor)
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
