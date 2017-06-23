""" analysis.py -- analysis tools for processing MFDn results file

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame

    6/2/15 (mac): Initiated (as mfdn_analysis.py).
    6/5/15 (mac): Restructure as subpackage.
    Last modified 7/6/15.

"""

import os
import math
import glob
import configparser # for band file
import numpy as np

# intra-package references
import mfdnres.res
import mfdnres.descriptor
import mfdnres.nonlinear

################################################################
# global constants
################################################################

# electromagnetic g factors (glp,gln,gsp,gsn)
MU_EM = np.array([1,0,5.586,-3.826])

################################################################
# load results files into database
################################################################

def import_res_files(
        data,filename,key_fields,filename_format,res_format,
        base_path=None,verbose=False
):
    """ Imports set of results files into given dictionary with customizable key tuple.

    Args:
        data (dict) : container for resulting MFDnRunData objects (key contents determined
            by key_fields argument)
        filename_list (str, list of str) : names of files or directories to process
        key_fields (tuple of str) : tuple of run descriptor keys to use to generate key to results
        filename_format (str) : determines which parser to use to parse run descriptor from filename
        res_format (str) : determines which parser to use to parse MFDn res file
        base_path (str) : base path to preappend to all filenamse (default: None)
        verbose (bool): verbose output (default: False)

    """

    # expand filename specification
    # upgrade single item to list
    if (type(filename) == list):
        raw_filename_list = filename
    else:
        raw_filename_list = [filename]
    # prepend path
    if (base_path is not None):
        qualified_filename_list = [
            os.path.join(base_path,filename)
            for filename in raw_filename_list
        ]
    else:
        qualified_filename_list = raw_filename_list
    # expand directory names to explicit lists of files
    filename_list = []
    for filename in qualified_filename_list:
        if (filename[-4:] == ".res"):
            # individual file
            filename_list.append(filename)
        else:
            # file prefix
            prefix_filename_list = glob.glob(filename + "*.res")
            filename_list.extend(prefix_filename_list)
            # directory
            directory_filename_list = glob.glob(os.path.join(filename,"*.res"))
            filename_list.extend(directory_filename_list)

    # process files
    for filename in filename_list:

        # determine key from filename
        if (verbose):
            print("Reading:",filename)
        res_file_info = mfdnres.descriptor.parse_res_filename(filename,filename_format=filename_format)
        ##key = tuple(map((lambda x : res_file_info[x]),key_fields))
        key = tuple([res_file_info[key] for key in key_fields])
        if (verbose):
            print(res_file_info["descriptor"],key)

        # generate new results object if key not already present
        if (key not in data):
            data[key] = mfdnres.res.MFDnRunData()

        # read data
        data[key].read_file(filename,res_format=res_format)

################################################################
# output tabulation: basic levels
################################################################

def write_level_table(results,filename,levels=None):
    """Writes basic level data.

        seq J p n T Eabs

    The field seq is included for historical reasons but written as a
    dummy zero.

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        levels (list of tuples): levels to include or None for all (default: None)

    """

    # determine list of levels to use
    if (levels is None):
        levels = results.get_levels()

    # assemble lines
    lines = []
    for qn in levels:
        line = (
            "{:1d} {:6.3f} {:1d} {:2d} {:6.3f} {:8.3f}"
            "\n"
        ).format(
            0,
            results.get_property(qn,"J"),
            results.get_property(qn,"g"),
            results.get_property(qn,"n"),
            results.get_property(qn,"T"),
            results.get_energy(qn)
        )
        lines.append(line)

    # write to file
    with open(filename,"wt") as fout:
        fout.writelines(lines)

################################################################
# output tabulation: radii
################################################################

def write_radii_table(results,filename,levels=None):
    """Writes radius data.

        seq J p n T rp rn r

    The field seq is included for historical reasons but written as a
    dummy zero.

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        levels (list of tuples): levels to include or None for all (default: None)

    """

    # determine list of levels to use
    if (levels is None):
        levels = results.get_levels()

    # assemble lines
    lines = []
    for qn in levels:
        line = (
            "{:1d} {:6.3f} {:1d} {:2d} {:6.3f} {:8.3f} {:8.3f} {:8.3f}"
            "\n"
        ).format(
            0,
            results.get_property(qn,"J"),
            results.get_property(qn,"g"),
            results.get_property(qn,"n"),
            results.get_property(qn,"T"),
            results.get_tbo(qn,"rp"),
            results.get_tbo(qn,"rn"),
            results.get_tbo(qn,"r")
        )
        lines.append(line)

    # write to file
    with open(filename,"wt") as fout:
        fout.writelines(lines)

################################################################
# output tabulation: angular momentum
################################################################

def effective_am(J_sqr):
    """  Converts mean square angular momentum to effective angular momentum.

    Args:
        J_sqr (float): value representing <J.J>

    Returns:
        (float): effective J
    """

    J = (math.sqrt(4*J_sqr+1)-1)/2
    return J

def write_level_am_table(results,filename,levels=None,default=np.nan):
    """Writes basic level data.

        seq J p n T Eabs

    The field seq is included for historical reasons but written as a
    dummy zero.

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        levels (list of tuples): levels to include or None for all (default: None)
        default (float): value to use for missing numerical entries (default: np.nan)

    """

    # determine list of levels to use
    if (levels is None):
        levels = results.get_levels()

    # assemble lines
    lines = []
    for qn in levels:
        # basic properties
        line = ("{:1d} {:6.3f} {:1d} {:2d} {:6.3f} {:8.3f}").format(
            0,
            results.get_property(qn,"J"),
            results.get_property(qn,"g"),
            results.get_property(qn,"n"),
            results.get_property(qn,"T"),
            results.get_energy(qn)
        )

        # angular momenta
        tbo = results.tbo[qn]
        value_sqr_list = [tbo.get("L2",default),tbo.get("Sp2",default),tbo.get("Sn2",default),tbo.get("S2",default)]
        value_list = list(map(effective_am,value_sqr_list))
        value_format = " {:6.3f}"
        line += (4*value_format).format(*value_list)

        # finalize line
        line += "\n"
        lines.append(line)

    # write to file
    with open(filename,"wt") as fout:
        fout.writelines(lines)


################################################################
# band configuration definition
################################################################

class BandDefinition(object):
    """Stores parameters of a band (or an arbitrary set of levels).

    Attributes from config file:

        [band]
        K (float): K quantum number
        levels (list of tuple): list of (J,g,n) quantum numbers

        [trans]
        Mj (dict): dictionary mapping J -> Mj
        signs (dict): dictionary mapping (J,Mj) -> sigma

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
    might only make use of the attributes levels and Mj).

    """

    ################################################################
    # initialization (and configuration file input)
    ################################################################

    def __init__(self,filename=None):
        """Initialize attributes from file.

        Configuration files follow Python configparser syntax.

        Example band configuration file contents:

            [band]
            K = 0.5
            levels =
                0.5 0 1
                1.5 0 1
                2.5 0 1
                3.5 0 1
                4.5 0 2

            [trans]
            Mj =
                0.5 0.5
                1.5 0.5
                2.5 0.5
                3.5 0.5
                4.5 0.5
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

        # read band
        self.K = None
        self.levels = []
        if (config.has_section("band")):

            # parse K
            if (config.has_option("band","K")):
                self.K = float(config["band"]["K"])

            # parse multiline list of levels
            if (config.has_option("band","levels")):
                levels_string = config["band"]["levels"]
                levels_string_lines = levels_string.strip().split("\n")
                for line in levels_string_lines:
                    entries = line.split()
                    qn = (float(entries[0]),int(entries[1]),int(entries[2]))
                    self.levels.append(qn)

        self.Mj = {}
        self.signs = {}
        if (config.has_section("trans")):

            # read band member Mj values for transitions
            if (config.has_option("trans","Mj")):
                Mj_string = config["trans"]["Mj"]
                Mj_string_lines = Mj_string.strip().split("\n")
                for line in Mj_string_lines:
                    entries = line.split()
                    (J, Mj) = (float(entries[0]),float(entries[1]))
                    self.Mj[J] = Mj

            # read band member signs
            if (config.has_option("trans","signs")):
                signs_string = config["trans"]["signs"]
                signs_string_lines = signs_string.strip().split("\n")
                for line in signs_string_lines:
                    entries = line.split()
                    (J,Mj,sigma) = (float(entries[0]),float(entries[1]),float(entries[2]))
                    self.signs[(J,Mj)] = sigma

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

        # derived
        # construct dictionary of band members indexed by J
        self.members = {}
        for qn in self.levels:
            (J,g,n) = qn
            self.members[J] = qn
        # construct sorted list of J values
        self.J_values = sorted(self.members.keys())

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

def write_band_table(results,filename,band,fields=None,default=np.nan):
    """Writes level energy, moment, and in-band transition data.

    With default arguments, recovers behavior of write_level_table.

    Each output line is written in the form

        seq J p n T Eabs <"E2_moments"> <"M1_moments"> <"E2_transitions_dJ1"> <"E2_transitions_dJ2"> <"M1_transitions_dJ1">

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

    Hint: If all the transitions come out as missing, have you set the correct Mj value?

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        band (BandDefinition): band definition
        fields (set of str): fields to include, else all fields written if None (default: None)
        default (float): value to use for missing numerical entries (default: np.nan)

    """

    # resolve special values of fields agument
    if (fields is None):
        fields = {"E2_moments","M1_moments","E2_transitions_dJ1","E2_transitions_dJ2","M1_transitions_dJ1"}

    # assemble table lines
    value_format = " {:9.4f}"  # format string for numerical values
    lines = []
    for J in band.J_values:

        # determine level
        qn = band.members[J]

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
                if ((Jf in band.members) and (Ji in band.Mj)):
                    # final level and appropriate Mj calculation are defined
                    Mj = band.Mj[Ji]
                    qnf = band.members[Jf]
                    values = results.get_rme(qnf,qni,op,Mj,default=default)
                    if ((Ji,Mj) in band.signs) and ((Jf,Mj) in band.signs):
                        # phases are defined for these states at the required Mj
                        values *= band.signs[(Ji,Mj)]*band.signs[(Jf,Mj)]
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

def write_network_table(results,filename,band):
    """Writes table of E2 RMEs

        Ji pi ni Ti Eabsi ; Jf pf nf Tf Eabsf ; |RME_p| |RME_n|

    Tabulation is of J-descending transitions only.  This may be
    generalized in the future.

    Args:
        results (MFDnRunData): results object containing the levels
        filename (str): output filename
        band (BandDefinition): band providing set of initial levels
           (and Mj values)

    """

    # assemble table lines
    value_format = " {:9.4f}"  # format string for numerical values
    lines = []

    # loop over initial states in band
    for qni in band.levels:

        # loop over all final states in results
        for qnf in results.get_levels():

            # test for transition to include in tabulation
            #   -- J-descending
            #   -- Mj defined for initial level
            #   -- transition available in calculations
            (Ji,_,_) = qni
            (Jf,_,_) = qnf
            allowed_sense = (Jf < Ji)
            op = "E2"
            Mj = band.Mj.get(Ji,None)
            available = results.has_rme(qnf,qni,op,Mj)
            if (not (allowed_sense and available)):
                continue

            # initiate line
            line = ""

            # initial state data
            qn = qni
            line += "{:6.3f} {:1d} {:2d} {:6.3f} {:8.3f}".format(
                results.get_property(qn,"J"),
                results.get_property(qn,"g"),
                results.get_property(qn,"n"),
                results.get_property(qn,"T"),
                results.get_energy(qn)
            )

            # final state data
            qn = qnf
            line += "{:6.3f} {:1d} {:2d} {:6.3f} {:8.3f}".format(
                results.get_property(qn,"J"),
                results.get_property(qn,"g"),
                results.get_property(qn,"n"),
                results.get_property(qn,"T"),
                results.get_energy(qn)
            )

            # value data
            values = np.abs(results.get_rme(qnf,qni,op,Mj))
            entries = 2
            line += (entries*value_format).format(*values)

            # finalize line
            line += "\n"
            lines.append(line)

    # write to file
    with open(filename,"wt") as fout:
        fout.writelines(lines)

################################################################
# extrapolation
################################################################

def extrapolate_energies_exp(Nmax_values,E_values,c2_guess=0.3,verbose=False):
    """Obtain exponentially extrapolated energy.

    Args:
        Nmax_values (NumPy vector): vector of Nmax values for the fit
        E_values (NumPy vector): vector of energy values for the fit
        c2_guess (float, optional): assumed exponential decay constant for initial linear fit
        verbose (bool, optional): whether or not to display fit progress

    Returns:
        (float): extrapolated energy value

    Example:

        >>> Nmax_values = np.array([6,8,10])
        >>> E_values = np.array([-34.849,-37.293,-38.537])
        >>> mfdnres.analysis.extrapolate_energies_exp(Nmax_values,E_values,verbose=True)

        Exponential fit
        Nmax values: [ 6  8 10]
        E values: [-34.849 -37.293 -38.537]
        Linear fit assuming c2 = 0.3
        (c0,c1): (-40.15798483275664, 32.030180719893316)
        Nonlinear fit returns: (array([-39.82661333,  37.74536425,   0.33765202]), 2)
    """

    if (verbose):
        print("Exponential fit")
        print("Nmax values:",Nmax_values)
        print("E values:",E_values)

    # do preliminary linear fit
    # based on typical assumed exponential decay constant
    A = np.array([
        [1,math.exp(-c2_guess*Nmax)]
        for Nmax in Nmax_values
    ])
    (c0_guess,c1_guess) = np.linalg.lstsq(A,E_values)[0]
    if (verbose):
        print("Linear fit assuming c2 = {}".format(c2_guess))
        print("(c0,c1):",(c0_guess,c1_guess))

    # do nonlinear fit
    c_guess = np.array([c0_guess,c1_guess,c2_guess])
    fit = mfdnres.nonlinear.fit(mfdnres.nonlinear.model_exp,Nmax_values,E_values,c_guess)
    if (verbose):
        print("Nonlinear fit returns:",fit)
    (c_values,_) = fit
    (c0,_,_) = c_values

    return c0

################################################################
# band fitting
################################################################

def band_fit_energy (results,band,verbose=False):
    """Fits band energies.

    Args:
        results (MFDnRunData): results object
        band (BandDefinition): band definition
        extrapolation (string, optional): keyword for stored energy extrapolation
            to retrieve (or None to use the unextrapolated energy)

    Returns:
        (np.array 1x3): band energy parameters (E0,A,a)

    """

    # construct energy vector
    b = np.array([
        results.get_energy(band.members[J])
        for J in band.J_list_energy
    ])
    if (verbose):
        print("J:",band.J_list_energy)
        print("Members:",band.members)
        print("Energies:",b)

    # construct coefficient matrix
    K = band.K
    if (K == 0.5):
        A = np.array([
            [1,J*(J+1),(-1)**(J+1/2)*(J+1/2)]
            for J in band.J_list_energy
        ])
    else:
        A = np.array([
            [1,J*(J+1)]
            for J in band.J_list_energy
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

def band_fit_E2 (results,band):
    """Extracts band E2 moments for normalization.

    Args:
        results (MFDnRunData): results object
        band (BandDefinition): band definition

    Returns:
        ((np.array 1x2)): band intrinsic quadrupole moments (Q0p,Q0n)

    """

    K = band.K
    J = band.J_Q
    qn = band.members[J]

    factor = (3*K**2-J*(J+1))/((J+1)*(2*J+3))
    if (factor == 0):
        raise ValueError("attempt to extract Q0 in case where Q(J) factor is 0")
    Q_values = np.array(results.moments.get((qn,"E2"),np.nan))
    Q0_values = Q_values/factor

    return Q0_values

def band_fit_M1 (results,band,verbose=False):
    """Extracts band M1 fit parameters.

    Args:
        results (MFDnRunData): results object
        band (BandDefinition): band definition

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
    K = band.K

    # accumulate moment entries
    A_moment = []
    b_moment = []
    for J in band.J_list_M1_moment:
        A_moment.append(f_moment(K,J))
        b_moment.append(results.moments[(band.members[J],"M1")])

    # accumulate transition entries
    A_trans = []
    b_trans = []
    for J in band.J_list_M1_trans:
        A_trans.append(f_trans(K,J))
        Ji = J
        Jf = J - 1
        Mj = band.Mj[Ji]
        values = np.array(results.get_rme(band.members[Jf],band.members[Ji],"M1",Mj))
        values *= band.signs[(Ji,Mj)]*band.signs[(Jf,Mj)]
        b_trans.append(values)

    # combine moment and transition arrays
    A = np.array(A_moment+A_trans,float)
    b = np.array(b_moment+b_trans,float)
    if (verbose):
        print("J_list_M1_moment:",band.J_list_M1_moment)
        print("J_list_M1_trans:",band.J_list_M1_trans)
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

def write_band_fit_parameters(results,filename,band,fields=None,verbose=False):
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
        band (BandDefinition): band definition
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
        parameters = band_fit_energy(results,band,verbose)
        line += (3*value_format).format(*parameters)

    # E2 parameters
    if ("E2" in fields):
        parameters = band_fit_E2(results,band)
        line += (2*value_format).format(*parameters)

    # M1 parameters
    if ("M1" in fields):
        parameters = band_fit_M1(results,band,verbose)
        flat_parameters = parameters.T.flatten()
        line += (12*value_format).format(*flat_parameters)

    # finalize line
    line += "\n"

    # write to file
    with open(filename,"wt") as fout:
        fout.write(line)



if (__name__ == "__main__"):
    pass
