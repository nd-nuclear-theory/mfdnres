""" ncci.py -- NCCI identities and definitions

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    05/31/19 (mac): Created.
    07/12/19 (mac): Add augment_params_with_parity().
    06/20/23 (mac): Pull in oscillator length functions from mcscript-ncci/utils.py.
    08/01/24 (mac): Add augment_params_with_Nmax().
    09/18/24 (mac): Fix missing constants from mcscript-ncci/constants.py.

"""

################################################################
# oscillator length calculations
################################################################

# Taken from mcscript-ncci/utils.py and mcscript-ncci/constants.py.

k_hbar_c  = 197.326_980_4     # (hbar c) in MeV fm [1,2]
k_mp_csqr = 938.272_088_16    # proton mass in MeV/c^2 [1,2]
k_mn_csqr = 939.565_420_52    # neutron mass in MeV/c^2 [1,2]
k_mN_csqr = (k_mp_csqr+k_mn_csqr)/2  # (m_N c^2) in MeV [1]

def oscillator_length(hw):
    """Calculate oscillator length for given oscillator frequency.

    b(hw) = (hbar c)/[(m_N c^2) (hbar omega)]^(1/2)

    Arguments:
        hw (numeric): hbar omega in MeV

    Returns:
        (float): b in fm
    """
    return k_hbar_c/math.sqrt(k_mN_csqr*hw)


def hw_from_oscillator_length(b):
    """Calculate oscillator frequency for given oscillator length.

    hw(b) = (hbar c)^2/[(m_N c^2) (b^2)]

    Arguments:
        b (numeric): oscillator length in fm

    Returns:
        (float): hbar omega in MeV
    """
    return k_hbar_c**2/(k_mN_csqr*b**2)


################################################################
# shell model N0
################################################################

def N0_for_nuclide(nuclide):
    """ Calculate quanta in lowest oscillator configuration for given nuclide.

    Natural parity grade can then be obtained as N0_for_nuclide(nuclide)%2.

    Inspired by spncci lgi::Nsigma0ForNuclide.

    Arguments:
        nuclide (tuple): (Z,N) for nuclide
    
    Returns:
        N0 (int): number of quanta
    """

    # each major shell eta=2*n+l (for a spin-1/2 fermion) contains (eta+1)*(eta+2) substates

    N0 = 0;
    for species_index in (0,1):
        num_particles = nuclide[species_index]
        eta=0
        while(num_particles>0):
            # add contribution from particles in shell
            shell_degeneracy = (eta+1)*(eta+2)
            num_particles_in_shell = min(num_particles,shell_degeneracy)
            N0 += num_particles_in_shell*eta

            # discard particles in shell
            num_particles -= num_particles_in_shell

            # move to next shell
            eta += 1

    return N0

def augment_params_with_Nmax(results_data):
    """Postprocess mesh point to add Nmax as parameter, with default value 0.

    Meant for traditional shell model runs in a harmonic oscillator valence shell.

    Arguments:

        results_data (ResultsData): results data object to augment

    """
    results_data.params.setdefault("Nmax", 0)

def augment_params_with_parity(results_data):
    """Postprocess mesh point to add parity as parameter, based on Nmax and nuclide.

    The parity parameter is already normally available in the results from mfdn
    but not after merging with results from obscalc-ob.

    The natural parity is defined by the number of oscillator quanta N0, as each
    oscillator quantum carries negative parity.

    Note that the grade can be extracted then as

        g = mfdnres.am.parity_grade(results_data.params["parity"])

    Arguments:

        results_data (ResultsData): results data object to augment

    """

    g = (results_data.params["Nmax"]+N0_for_nuclide(results_data.params["nuclide"]))%2
    results_data.params.update({
        "parity" : (-1)**g
    })


################################################################
# test code
################################################################

if (__name__=="__main__"):

    # test N0_for_nuclide
    for nuclide in [(1,1),(2,2),(3,3),(4,3),(3,4),(8,8),(9,8)]:
        N0 = N0_for_nuclide(nuclide)
        print("nuclide {} N0 {}".format(nuclide,N0))
