"""Spectroscopic amplitude and spectroscopic factor observables.

    Mark A. Caprio
    University of Notre Dame

    - 07/06/25 (mac): Created.

"""

import numpy as np

import mfdnres
import mfdnres.am
import mfdnres.data
import mfdnres.ncci
import mfdnres.observable
import mfdnres.ticks
import mfdnres.tools


class SpectroscopicAmplitude(mfdnres.observable.Observable):
    """ Observable extractor for spectroscopic amplitude.

    """

    def __init__(self, nuclide, delta_nuclide, levelf, leveli, orbital):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            delta_nuclide (tuple): (delta_Z, delta_N) identifier for
            spectroscopic amplitude/factor.

            levelf (LevelSelector): Level selector for final level.

            leveli (LevelSelector): Level selector for initial level.

            orbital (tuple): Orbital (n,l,j) for spectroscopic amplitude.

        """
        super().__init__()
        self._nuclide = nuclide  # DEPRECATED
        self._delta_nuclide = delta_nuclide
        self._level_pair = levelf, leveli
        self._orbital = orbital

        # deduce final nuclide
        final_nuclide = (nuclide[0]+delta_nuclide[0], nuclide[1]+delta_nuclide[1])
        self._nuclide_pair = (final_nuclide, nuclide)
        
    def value(self, results_data):
        """ Extract observable.
        """
        # TODO (mac): Resolve level selection on differeing initial and final results data.
        ## qn_pair = self._level_pair[0].select_level(results_data_final), self._level_pair[1].select_level(results_data_initial)
        qn_pair = self._level_pair[0], self._level_pair[1].select_level(results_data)  # INTERIM
        if (qn_pair[0] is None) or (qn_pair[1] is None):
            return np.nan
        amplitudes = results_data.get_spectroscopic_amplitudes(self._delta_nuclide, qn_pair)
        if amplitudes is None:
            return np.nan
        amplitude = amplitudes.get(self._orbital)
        if amplitude is None:
            return np.nan
        return amplitude

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            mfdnres.tools.nuclide_str(self._nuclide_pair[0]),
            mfdnres.tools.nuclide_str(self._nuclide_pair[1]),
            "spamp",
            ##self._operator,
            mfdnres.tools.qn_str(self._level_pair[0]),  ## self._level_pair[0].descriptor_str,  # interim
            self._level_pair[1].descriptor_str,
            "{:d}-{:d}-{:.1f}".format(*self._orbital),
        ])

    @property
    def nuclide_label_text(self):
        """Formatted LaTeX text representing nuclide.
        """
        nuclide_label_texts = data.isotope(self._nuclide_pair[0]), data.isotope(self._nuclide_pair[1])
        label = r"{}\rightarrow{}".format(
            nuclide_label_texts[1],
            nuclide_label_texts[0],
        )
        return label
    
    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        ## observable_text = r"Q_{{0{}}}".format(OBSERVABLE_TAG_STR_BY_OPERATOR[self._operator])
        ## level_text = self._level.label_text
        ## label = r"{}({})".format(observable_text, level_text)
        n, l, j = self._orbital
        label = "a({:d},{:d},{:s})".format(n,l,mfdnres.ticks.half_int_str(j))
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"a"
        units_text = None
        return observable_text, units_text


class SpectroscopicFactor(mfdnres.observable.Observable):
    """Observable extractor for spectroscopic factor.

    Calculations follow equations (1)-(4) of Sargsyan et al., PRC 108, 054303
    (2023), defined for a creation operator from the mass A-1 system to the mass
    A system.  Under our convention of storing a transition with the initial
    nuclide, this means A-1=sum(nuclide).

    """

    def __init__(self, nuclide, delta_nuclide, levelf, leveli, lj):
        """Initialize with given parameters.

        Arguments:

            nuclide (tuple): (Z, N)

            delta_nuclide (tuple): (delta_Z, delta_N) identifier for
            spectroscopic amplitude/factor.

            levelf (LevelSelector): Level selector for final level.

            leveli (LevelSelector): Level selector for initial level.

            lj (tuple): Angular momentum channel (l,j) for spectroscopic factor.

        """
        super().__init__()
        self._nuclide = nuclide  # DEPRECATED
        self._delta_nuclide = delta_nuclide
        self._level_pair = levelf, leveli
        self._lj = lj

        # deduce final nuclide
        final_nuclide = (nuclide[0]+delta_nuclide[0], nuclide[1]+delta_nuclide[1])
        self._nuclide_pair = (final_nuclide, nuclide)

    def value(self, results_data):
        """ Extract observable.
        """

        if self._delta_nuclide not in {(+1,0), (0,+1)}:
            raise ValueError("Unexpected delta_nuclide.  Only nucleon addition presently supported by implemented formulas.")
        
        # TODO (mac): Resolve level selection on differeing initial and final results data.
        ## qn_pair = self._level_pair[0].select_level(results_data_final), self._level_pair[1].select_level(results_data_initial)
        qn_pair = self._level_pair[0], self._level_pair[1].select_level(results_data)  # INTERIM
        if (qn_pair[0] is None) or (qn_pair[1] is None):
            return np.nan
        amplitudes = results_data.get_spectroscopic_amplitudes(self._delta_nuclide, qn_pair)
        if amplitudes is None:
            return np.nan

        # accumulate spectroscopic factor sum
        mass_ratio = sum(self._nuclide_pair[0])/sum(self._nuclide_pair[1])
        Jf, _, _ = qn_pair[0]
        am_factor = 1/(2*Jf+1)
        S = 0
        for orbital, amplitude in amplitudes.items():
            n, l, j = orbital
            if (l, j) != self._lj:
                continue
            N = 2*n+l
            S += mass_ratio**N * am_factor * amplitude**2
        return S

    @property
    def descriptor_str(self):
        """ Text string describing observable.
        """
        return "-".join([
            mfdnres.tools.nuclide_str(self._nuclide_pair[0]),
            mfdnres.tools.nuclide_str(self._nuclide_pair[1]),
            "sf",
            ##self._operator,
            mfdnres.tools.qn_str(self._level_pair[0]),  ## self._level_pair[0].descriptor_str,  # interim
            self._level_pair[1].descriptor_str,
            "{:d}-{:.1f}".format(*self._lj),
        ])

    @property
    def nuclide_label_text(self):
        """Formatted LaTeX text representing nuclide.
        """
        nuclide_label_texts = data.isotope(self._nuclide_pair[0]), data.isotope(self._nuclide_pair[1])
        label = r"{}\rightarrow{}".format(
            nuclide_label_texts[1],
            nuclide_label_texts[0],
        )
        return label
    
    @property
    def nuclide_set(self):
        """Set of nuclides entering calculation of observable.
        """
        return set(self._nuclide_pair)
    
    @property
    def observable_label_text(self):
        """ Formatted LaTeX text representing observable.
        """
        l, j = self._lj
        label = "S({:d},{:s})".format(l,mfdnres.ticks.half_int_str(j))
        return label

    @property
    def axis_label_text(self):
        """ Formatted LaTeX text representing axis label.
        """
        observable_text = r"S"
        units_text = None
        return observable_text, units_text

    
################################################################
# main
################################################################

def main():
    pass

if __name__ == "__main__":
    main()
