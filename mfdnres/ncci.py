""" ncci.py -- NCCI identities and definitions

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    05/31/19 (mac): Created.
    07/12/19 (mac): Add augment_params_with_parity().

"""

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

def augment_params_with_parity(results_data):
    """Postprocess mesh point to add parity as parameter, based on Nmax and nuclide.

    The parity parameter is already normally available in the results from mfdn
    but not after merging with results from obscalc-ob.

    Arguments:
        results_data (ResultsData): results data object to augment
    """

    results_data.params.update({
        "parity":(-1)**(results_data.params["Nmax"]+N0_for_nuclide(results_data.params["nuclide"]))
    })


################################################################
# test code
################################################################

if (__name__=="__main__"):

    # test N0_for_nuclide
    for nuclide in [(1,1),(2,2),(3,3),(4,3),(3,4),(8,8),(9,8)]:
        N0 = N0_for_nuclide(nuclide)
        print("nuclide {} N0 {}".format(nuclide,N0))
