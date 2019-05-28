""" band.py -- rotational band analysis

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    - 07/30/17 (mac): Extract band analysis functions from analysis.py (originated
         6/2/15) to band.py.
    - 04/27/18 (mac): Rename parameter Mj to M.
    - 10/26/18 (mac): Update band energy fitting to accommodate mfdnres mesh reorganization.
        + Update configuration file format (need M value for energy retrieval, need band parity).
        + Revise energy fitting to take energies as dictionary.
    - 12/12/18 (mac): Reimplement BandDefinition derived data structures as Python properties.
    - 04/25/19 (mac): Remove hard-coded workaround for early SpNCCI operator implementation in
      write_network_table, and add transition_operators optional argument.
"""

import math
import configparser # for band file

import numpy as np

################################################################
# band configuration definition
################################################################

class BandDefinition(object):
    """Stores parameters of a band (or an arbitrary set of levels).

    Attributes from config file:

        [band]
        K (float): K quantum number
        g (int): parity grade [P=(-)^g]
        levels (list of tuple): list of (J,g,n) quantum numbers
        M (dict): dictionary mapping J -> M

        [trans]
        signs (dict): dictionary mapping (J,M) -> sigma

        [fit]
        J_list_energy (list of float): J values to use for energies in energy fit
        J_Q (float): J value to use for normalization of E2 intrinsic matrix element
        J_list_M1_moment (list of float): J values to use for M1 moments in M1 fit
        J_list_M1_trans (list of float): J values to use for initial levels for M1 transitions in M1 fit

    Attributes derived from these:
        members (dict): dictionary mapping J -> (J,g,n)
        J_values (list of float): sorted list of J values

    While the BandDefinition class is primarily intended to represent
    traditional bands for rotational analysis, it can also be used to
    represent a generalized set of states for certain purposes.  For
    traditional rotational analysis, the given levels must be unique
    by J (this is true for any analysis which uses the members
    dictionary to look up states by J).  However, code can be written
    which does not rely on this property (e.g., for compiling a list
    of transitions from a set of levels, for a "network diagram", one
    might only make use of the attributes levels and M).

    """

    ################################################################
    # initialization (and configuration file input)
    ################################################################

    def __init__(self,filename=None):
        """Initialize attributes from file.

        Configuration files follow Python configparser syntax.

        Changes for 2017-2018 mfdnres restucturing: Move M listing from [trans]
        to [band].  Add parity field to [band] and remove parity from level qn.
        Wave function signs may move out of trans section, as well, when return
        to transition analysis with preprocessor-based transition handling.

        Example band configuration file contents:

            [band]
            K = 0.5
            g = 1
            levels =
                0.5 1 1
                1.5 1 1
                2.5 1 1
                3.5 1 1
                4.5 1 2
            M =
                0.5 0.5
                1.5 0.5
                2.5 0.5
                3.5 0.5
                4.5 0.5

            [trans]
            signs =
                0.5 0.5 +1
                1.5 0.5 +1
                2.5 0.5 +1
                3.5 0.5 +1
                4.5 0.5 +1

            [fit]
            J_list_energy = 0.5 1.5 3.5
            J_Q = 1.5
            J_list_M1_moment = 0.5 1.5 2.5 3.5
            J_list_M1_trans = 0.5 1.5 2.5 3.5


        Args:
            filename (str): filename of band config file (default: None)

        """

        # quit if no filename
        # it is up to user to initialize any necessary fields
        if (filename is None):
            return

        # read configuration
        config = configparser.ConfigParser()
        config.read_file(open(filename))

        # read band composition parameters
        self.K = None
        self.g = None
        self.levels = []
        self.M = {}

        if (config.has_section("band")):

            # parse K
            if (config.has_option("band","K")):
                self.K = float(config["band"]["K"])

            # parse parity
            if (config.has_option("band","g")):
                self.g = int(config["band"]["g"])

            # parse multiline list of levels
            if (config.has_option("band","levels")):
                levels_string = config["band"]["levels"]
                levels_string_lines = levels_string.strip().split("\n")
                for line in levels_string_lines:
                    entries = line.split()
                    qn = (float(entries[0]),int(entries[1]),int(entries[2]))
                    self.levels.append(qn)

            # read band member M values for level energies
            if (config.has_option("band","M")):
                M_string = config["band"]["M"]
                M_string_lines = M_string.strip().split("\n")
                for line in M_string_lines:
                    entries = line.split()
                    (J, M) = (float(entries[0]),float(entries[1]))
                    self.M[J] = M

        # read transition parameters
        self.signs = {}
        if (config.has_section("trans")):

            # read band member signs
            if (config.has_option("trans","signs")):
                signs_string = config["trans"]["signs"]
                signs_string_lines = signs_string.strip().split("\n")
                for line in signs_string_lines:
                    entries = line.split()
                    (J,M,sigma) = (float(entries[0]),float(entries[1]),float(entries[2]))
                    self.signs[(J,M)] = sigma

        # read fit parameters
        self.J_list_energy = []
        self.J_Q = None
        self.J_list_M1_moment = []
        self.J_list_M1_trans = []
        if (config.has_section("fit")):
            if (config.has_option("fit","J_list_energy")):
                J_list_energy_string = config["fit"]["J_list_energy"]
                self.J_list_energy = list(map(float,J_list_energy_string.split()))
            if (config.has_option("fit","J_Q")):
                self.J_Q = float(config["fit"]["J_Q"])
            if (config.has_option("fit","J_list_M1_moment")):
                J_list_M1_moment_string = config["fit"]["J_list_M1_moment"]
                self.J_list_M1_moment = list(map(float,J_list_M1_moment_string.split()))
            if (config.has_option("fit","J_list_M1_trans")):
                J_list_M1_trans_string = config["fit"]["J_list_M1_trans"]
                self.J_list_M1_trans = list(map(float,J_list_M1_trans_string.split()))

    ################################################################
    # derived properties
    ################################################################

    # It is desirable to derive these listings on the fly, instead of creating
    # them in _init_, so that they reflect any subsequent changes to self.levels.

    @property
    def members(self):
        """ Dictionary of band members indexed by J. """
        members_by_J = {}
        for qn in self.levels:
            (J,g,n) = qn
            members_by_J[J] = qn
        return members_by_J

    @property
    def J_values(self): 
        """ Sorted list of J values. """
        return sorted(self.members.keys())

    ################################################################
    # container-like access
    ################################################################

    def __iter__(self):
        """Make band iterable, allowing iteration over full list of state quantum numbers.

        A BandDefinition object can therefore be used in many places
        where a list of levels may be used, e.g., as the levels
        argument to write_level_table.

        Duplication of J values is allowed.

        Returns:
            (iter) : iterator over levels list

        """

        return iter(self.levels)


################################################################
# output tabulation: in-band moment and transition data
################################################################

def write_band_table(results,filename,band_definition,fields=None,default=np.nan):
    """Writes level energy, moment, and in-band transition data.

    With default arguments, recovers behavior of write_level_table.

    Each output line is written in the form

        "seq"=0 J p n T Eabs <"E2_moments"> <"M1_moments"> <"E2_transitions_dJ1"> <"E2_transitions_dJ2"> <"M1_transitions_dJ1">

    The field seq is included for historical reasons but written as a
    dummy zero.

    Any missing or undefined transition RMEs are written as a NaN (or
    as the numerical value given by the argument default).

    Entries within a field are written in the order RME(lp,ln,sp,sn)
    for M1 or RME(p,n) for E2.  Note that the electromagnetic M1
    moment is *not* included in the list of moments, even though it
    was in the original tabulation format for berotor.  The
    electromagnetic M1 RME can be recovered by taking the linear
    combination with the standard gyromagnetic ratios as coefficients.

    Some of these fields may optionally be omitted.

    The given levels are assumed to be unique by J.

    RME lookup is done with lower-J state as final state (presumes
    lower-J state served as a reference state in the MFDn
    calculation).

    Hint: If all the transitions come out as missing, have you set the correct M value?

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        band_definition (BandDefinition): band definition
        fields (set of str): fields to include, else all fields written if None (default: None)
        default (float): value to use for missing numerical entries (default: np.nan)

    """

    # resolve special values of fields argument
    if (fields is None):
        fields = {"E2_moments","M1_moments","E2_transitions_dJ1","E2_transitions_dJ2","M1_transitions_dJ1"}

    # assemble table lines
    value_format = " {:9.4f}"  # format string for numerical values
    lines = []
    for J in band_definition.J_values:

        # determine level
        qn = band_definition.members[J]

        # initial state data
        line = "{:1d} {:6.3f} {:1d} {:2d} {:6.3f} {:8.3f}".format(
            0,
            results.get_property(qn,"J"),
            results.get_property(qn,"g"),
            results.get_property(qn,"n"),
            results.get_property(qn,"T"),
            results.get_energy(qn)
        )

        # loop over moment fields
        #    (field,op,entries)
        field_definitions = [
            ("E2_moments","E2",2),
            ("M1_moments","M1",4)
        ]
        for (field,op,entries) in field_definitions:
            if (field in fields):
                values = default*np.ones(entries)
                if ((qn,op) in results.moments):
                    # values exist
                    values = np.array(results.moments[(qn,op)])
                line += (entries*value_format).format(*values)

        # loop over transition fields
        #    (field,op,dJ,entries)
        field_definitions = [
            ("E2_transitions_dJ1","E2",1,2),
            ("E2_transitions_dJ2","E2",2,2),
            ("M1_transitions_dJ1","M1",1,4)
        ]
        for (field,op,dJ,entries) in field_definitions:
            if (field in fields):
                Ji = J
                Jf = Ji-dJ
                qni = qn
                values = default*np.ones(entries)
                if ((Jf in band_definition.members) and (Ji in band_definition.M)):
                    # final level and appropriate M calculation are defined
                    M = band_definition.M[Ji]
                    qnf = band_definition.members[Jf]
                    values = results.get_rme(qnf,qni,op,M,default=default)
                    if ((Ji,M) in band_definition.signs) and ((Jf,M) in band_definition.signs):
                        # phases are defined for these states at the required M
                        values *= band_definition.signs[(Ji,M)]*band_definition.signs[(Jf,M)]
                line += (entries*value_format).format(*values)

        # finalize line
        line += "\n"
        lines.append(line)

    # write to file
    with open(filename,"wt") as fout:
        fout.writelines(lines)

################################################################
# output tabulation: RME network from band members
################################################################

def write_network_table(
        results,filename,band_definition,
        energy_cutoff=None,
        transition_operators=["Qpintr","Qnintr","Qintr"]
):
    """Writes table of E2 RMEs

    WARNING: Currently adapted for spncci use.  Must generalize to recover MFDn
    use.  Currently puts mass Q operator in place of both proton and neutron
    operators.

    Data format:

        J_f gex_f n_f Eabs_f ; J_i p_i n_i Eabs_i ; |RME_p| |RME_n|

    Legacy data format (2015-era) -- included isospin, swapped initial/final label order:

        J_i gex_i n_i T_i Eabs_i ; J_f p_f n_f T_f Eabs_f ; |RME_p| |RME_n|

    Tabulation is of J-descending transitions only.  This may be
    generalized in the future.

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        band_definition (BandDefinition): band providing set of initial levels
           (and M values)
        energy_cutoff (float,optional): energy cutoff to limit output size
        transition_operators (list of str, optional): operator identifiers for transition operators to tabulate
    """

    # assemble table lines
    value_format = " {:9.4f}"  # format string for numerical values
    lines = []

    # loop over initial states in band
    for qni in band_definition.levels:

        # loop over all final states in results
        for qnf in results.levels:

            # test for transition to include in tabulation
            #   -- J-descending
            #   -- M defined for initial level
            #   -- transition available in calculations
            (Ji,gexi,ni) = qni
            (Jf,gex,nf) = qnf
            allowed_sense = (Jf < Ji)
            if (not allowed_sense):
                continue

            # MFDn:
            ## op = "E2"
            ## M = band_definition.M.get(Ji,None)
            ## available = results.get_rme(qnf,qni,op,M)

            # retrieve values
            rme_values = [
                results.get_rme(operator,(qnf,qni))
                for operator in transition_operators
            ]

            # short circuit if all values are nans
            all_nan = True
            for is_nan in map(np.isnan,rme_values):
                all_nan = all_nan and is_nan
            if (all_nan):
                continue

            # initiate line
            line = ""

            # final state data
            energy = results.get_energy(qnf)
            if ((energy_cutoff is not None) and (energy>energy_cutoff)):
                continue
            line += "{qn[0]:4.1f} {qn[1]:1d} {qn[2]:2d} {energy:8.3f}    ".format(
                qn=qnf,
                energy=energy
            )

            # initial state data
            energy = results.get_energy(qni)
            if ((energy_cutoff is not None) and (energy>energy_cutoff)):
                continue
            line += "{qn[0]:4.1f} {qn[1]:1d} {qn[2]:2d} {energy:8.3f}    ".format(
                qn=qni,
                energy=energy
            )

            # value data
            ## values = np.abs(results.get_rme(qnf,qni,op,M))
            num_values = len(rme_values)
            line += (num_values*value_format).format(*rme_values)

            # finalize line
            line += "\n"
            lines.append(line)

    # write to file
    with open(filename,"wt") as fout:
        fout.writelines(lines)


################################################################
# band fitting
################################################################

def band_fit_energy(band_definition,level_energy_dict,verbose=False):
    """Obtain band energy parameters, by fitting band energies.

    Args:
        band_definition (BandDefinition): band providing set of initial levels
        level_energy_dict (dict): level energies, as mapping J->energy

    Returns:
        (np.array 1x3): band energy parameters (E0,A,a)

    """

    # construct energy vector
    b = np.array([
        level_energy_dict[J]
        for J in band_definition.J_list_energy
    ])
    if (verbose):
        print("J:",band_definition.J_list_energy)
        print("Members:",band_definition.members)
        print("Energies:",b)

    # construct coefficient matrix
    K = band_definition.K
    if (K == 0.5):
        A = np.array([
            [1,J*(J+1),(-1)**(J+1/2)*(J+1/2)]
            for J in band_definition.J_list_energy
        ])
    else:
        A = np.array([
            [1,J*(J+1)]
            for J in band_definition.J_list_energy
        ])

    # solve system
    parameters = np.linalg.lstsq(A,b)[0]
    if (verbose):
        print("Parameters:",parameters)

    # upgrade parameter vector to three entry
    #   zero pads for parameter a if not already present
    if (parameters.shape == (2,)):
        parameters = np.append(parameters,[0],axis=0)

    # convert last coefficient to a = c3/c2
    parameters[2] /= parameters[1]

    return parameters

def band_fit_E2(results,band_definition):
    """Extracts band E2 moments for normalization.

    CAVEAT: Requires updating to properly handle mesh with distinct M values.

    Args:
        results (MFDnRunData): results object
        band_definition (BandDefinition): band definition

    Returns:
        ((np.array 1x2)): band intrinsic quadrupole moments (Q0p,Q0n)

    """

    K = band_definition.K
    J = band_definition.J_Q
    qn = band_definition.members[J]

    factor = (3*K**2-J*(J+1))/((J+1)*(2*J+3))
    if (factor == 0):
        raise ValueError("attempt to extract Q0 in case where Q(J) factor is 0")
    Q_values = np.array(results.moments.get((qn,"E2"),np.nan))
    Q0_values = Q_values/factor

    return Q0_values

def band_fit_M1(results,band_definition,verbose=False):
    """Extracts band M1 fit parameters.

    CAVEAT: Requires updating to properly handle mesh with distinct M values.

    Args:
        results (MFDnRunData): results object
        band_definition (BandDefinition): band definition

    Returns:
        (np.array 3x4): band M1 fit parameters, or all np.nan for blatantly undersized systems

    This is the array of least squares solution vectors to the linear
    systems, arranged as

       [[a0,a0,a0,a0],
        [a1,a1,a1,a1],
        [a2,a2,a2,a2]]

    where successive columns are for Dlp, Dln, Dsp, Dsn, respectively,
    and a2 is set to 0 if K is not 1/2.

    Derived from code in rotor_fit.py (mac, 5/3/14).
    """

    # helper functions
    def f_moment(K,J):
        """ Generate coefficient matrix row for moments.
        """
        f0 = J
        f1 = K/(J+1)
        f2 = (-1)**(J-0.5) / (2*math.sqrt(2)) * (2*J+1)/(J+1)
        if (K == 0):
            # case K = 0
            coefficients = [f0]
        elif (K == 1/2):
            # case K = 1/2
            coefficients = [f0,f1,f2]
        else:
            # case K generic
            coefficients = [f0,f1]

        return coefficients
    def f_trans(K,J):
        """ Generate coefficient matrix row for transitions.
        """
        f0 = 0
        f1 = -math.sqrt(3/(4*math.pi)) * math.sqrt((J**2-K**2)/J)
        f2 = (-1)**(J-0.5) / math.sqrt(2) * f1
        if (K == 0):
            # case K = 0
            coefficients = [f0]
        elif (K == 1/2):
            # case K = 1/2
            coefficients = [f0,f1,f2]
        else:
            # case K generic
            coefficients = [f0,f1]

        return coefficients

    # setup
    K = band_definition.K

    # accumulate moment entries
    A_moment = []
    b_moment = []
    for J in band_definition.J_list_M1_moment:
        A_moment.append(f_moment(K,J))
        b_moment.append(results.moments[(band_definition.members[J],"M1")])

    # accumulate transition entries
    A_trans = []
    b_trans = []
    for J in band_definition.J_list_M1_trans:
        A_trans.append(f_trans(K,J))
        Ji = J
        Jf = J - 1
        M = band_definition.M[Ji]
        values = np.array(results.get_rme(band_definition.members[Jf],band_definition.members[Ji],"M1",M))
        values *= band_definition.signs[(Ji,M)]*band_definition.signs[(Jf,M)]
        b_trans.append(values)

    # combine moment and transition arrays
    A = np.array(A_moment+A_trans,float)
    b = np.array(b_moment+b_trans,float)
    if (verbose):
        print("J_list_M1_moment:",band_definition.J_list_M1_moment)
        print("J_list_M1_trans:",band_definition.J_list_M1_trans)
        print("Coefficient matrix")
        print(A)
        print("Ordinate matrix")
        print(b)

    # trap blatantly insufficient system
    if ( not (
            ((K==0) and (len(b)>=1))
            or
            ((K==1/2) and (len(b)>=3))
            or
            ((K>1/2) and (len(b)>=2))
    ) ):
            parameters = np.nan*np.ones((3,4))
            return parameters

    # solve system
    parameters = np.linalg.lstsq(A,b)[0]

    # upgrade parameter matrix to three rows
    #   zero pads for parameter a2 if not already present
    if (parameters.shape == (1,4)):
        parameters = np.append(parameters,[[0,0,0,0],[0,0,0,0]],axis=0)
    elif (parameters.shape == (2,4)):
        parameters = np.append(parameters,[[0,0,0,0]],axis=0)
    if (verbose):
        print("Parameter matrix")
        print(parameters)

    return parameters

def write_band_fit_parameters(results,filename,band_definition,fields=None,verbose=False):
    """Writes band fit parameters.

    The output contains one line, of the form:

        <"energy"> <"E2"> <"M1">

    where the fields are

        <"energy">: E0 A a
        <"E2">: Q0p Q0n
        <"M1"> alp0 alp1 alp2 ; aln0 ... ; asp0 ... ; asn0 ...

    Some of these fields may optionally be omitted.

    Args:
        results (MFDnRunData): results object
        filename (str): output filename
        band_definition (BandDefinition): band definition
        fields (set of str): fields to include, else all fields written if None (default: None)

    The fields argument is mostly useful to drop EM band parameters (fields={"energy"})if
    attention is only being paid to energy fitting.

    """

    # formatting configuration
    value_format = " {:9.4f}"  # format string for numerical values

    # resolve special values of fields agument
    if (fields is None):
        fields = {"energy","E2","M1"}

    # set up line
    line = ""

    # energy parameters
    if ("energy" in fields):
        parameters = band_fit_energy(results,band_definition,verbose)
        line += (3*value_format).format(*parameters)

    # E2 parameters
    if ("E2" in fields):
        parameters = band_fit_E2(results,band_definition)
        line += (2*value_format).format(*parameters)

    # M1 parameters
    if ("M1" in fields):
        parameters = band_fit_M1(results,band_definition,verbose)
        flat_parameters = parameters.T.flatten()
        line += (12*value_format).format(*flat_parameters)

    # finalize line
    line += "\n"

    # write to file
    with open(filename,"wt") as fout:
        fout.write(line)

if (__name__ == "__main__"):
    pass
