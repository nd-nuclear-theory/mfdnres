""" mfdn_v14b05.py -- provide res file parser for MFDn version 14 beta 05 with VXX input

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame
    7/26/15 (mac): Adapted from res_parser_v14b06:
        -- Change expected positioning and format for basis parameter lines.
        -- Set up for generic input file format.
        -- Make less strict on expected name of radius observables.
        -- Remove unneeded special case patch for mac early squared am operator runs.
        -- AD HOC: Removed TBO block, but in general will need.  Better to restucture input
          to block identification loop and block handlers.
    10/10/16 (mac): Update to return list containing single mesh point.
    04/27/18 (mac): Rename parameter Mj to M.

"""

import re

# intra-package references
from .. import (
    mfdn_results_data_v14,
    input,
    tools,
    )

def parser(fin,verbose):
    """ Read result file data into MFDnResultsData objects.

    Args:
        fin (stream): results file to read
        verbose (bool): verbose output for debugging

    """

    ################################################################
    # set up container
    ################################################################

    results = mfdn_results_data_v14.MFDnResultsDataV14()

    ################################################################
    # read header
    ################################################################

    tools.parse_line(fin.readline(),r"")
    tools.parse_line(fin.readline(),r"OUTPUT from MFD-nuclear-physics Version 14")
    tools.parse_line(fin.readline(),r"")
    match = tools.parse_line(fin.readline(),r"hbar-omega\s+(?P<hw>[\-0-9.]+)")
    results.params["hw"] = float(match.group("hw"))
    tools.parse_line(fin.readline(),r"nucleon mass\s+(?P<mN>[\-0-9.]+)")
    tools.parse_line(fin.readline(),r"")
    match = tools.parse_line(fin.readline(),r"2-body interaction data in (?P<format>\S+) format")
    interaction_format = match.group("format")
    if (interaction_format == "_vxx_v4"):
        # VXX format
        # eat Hamiltonian input lines
        while (fin.readline().strip() != ""):
            pass
    else:
        raise ValueError("Hamiltonian file format mode {} not recognized by MFDn results parser".format(interaction_format))

    tools.parse_line(fin.readline(),r"INPUT: read in TBME files")
    tools.parse_line(fin.readline(),r"TBME binary file  ([^\.]+)\.bin")  # radius observable

    # read names of additional two-body observables
    tbo_operators = []
    done = False
    while (not done):
        line=fin.readline().strip()
        if (line == ""):
            done = True
            continue

        match = tools.parse_line(line,r"TBME binary file\s+tbme-(\w+)")
        operator_name = match.group(1)
        tbo_operators.append(operator_name)

    tools.parse_line(fin.readline(),r"")

    # Z and N
    match = tools.parse_line(fin.readline(),r"Number of protons \(and neutrons\)\s+(?P<Z>\d+)\s+(?P<N>\d+)")
    results.params["nuclide"] = (int(match.group("Z")),int(match.group("N")))

    # Nmin and Nmax
    match = tools.parse_line(fin.readline(),r"Number of HO quanta above minimal conf\s+(?P<Nmin>\d+)\s+to\s+(?P<Nmax>\d+)")
    results.params["Nmin"] = int(match.group("Nmin"))
    results.params["Nmax"] = int(match.group("Nmax"))

    tools.parse_line(fin.readline(),"")

    ################################################################
    # read static properties
    ################################################################

    tools.parse_line(fin.readline(),r"Summary of static properties")
    tools.parse_line(fin.readline(),r"")

    # set up data structure to hold states on input
    state_list = []  # states within current run
    state_lookup = {}  # lookup from current run's sequence number to (J,g,n) indices

    # occupation group -- provides basic state data
    tools.parse_line(fin.readline(),r"Seq    J    NP  n    T        Eabs      Occupation Prob.")
    done = False
    while (not done):
        # read line (until blank line reached)
        line=fin.readline().strip()
        if (line == ""):
            done = True
            continue
        if (verbose):
            print(">>> Occupations:",line)


        # read state data
        match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<Eabs>\S+)\s+(?P<occupations>.*)")
        seq = int(match.group("seq"))  # do not save seq in state info, since seq is not unique when runs are combined
        J = float(match.group("J"))  # floating point J (might not be converged to half integer!)
        g = int(match.group("g"))    # relative parity "grade", i.e., as integer in {0,1}
        n = int(match.group("n"))    # MFDn-assigned ordering number within (J,g) sector
        T = float(match.group("T"))  #  floating point isospin T
        E = float(match.group("Eabs"))  #  floating point energy eigenvalue
        qn = (J,g,n)

        # set up state instance
        # initialize new state
        state = mfdn_results_data_v14.MFDnStateData(qn,T,E)
        # save state to working list for this res file
        state_list.append(state)
        # save lookup information (by sequence number) for processing transitions below
        state_lookup[seq] = qn

        # parse occupations
        state.occupations = list(map(float,match.group("occupations").split()))

        # save observables to dictionaries
        results.properties[state.qn] = {"J":J,"g":g,"n":n,"T":T}
        results.energies[state.qn] = E


    # extract poor man's M
    # needed below when we store transitions (it would be better if M were stated in MFDn res file)
    M = min([state.properties["J"] for state in state_list])

    # amplitude group
    tools.parse_line(fin.readline(),r"Seq    J    NP  n    T        Eabs  \(amp\(N\)\)\^2 for N=")
    for state in state_list:
        # read line
        line=fin.readline().strip()
        if (verbose):
            print(">>> Amplitudes:",line)

        # read state data
        match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<Eabs>\S+)\s+(?P<amplitudes>.*)")
        state.amplitudes = list(map(float,match.group("amplitudes").split()))

        # make sure TBO dictionary available
        results.tbo.setdefault(state.qn,{})

    tools.parse_line(fin.readline(),r"")

    # radii TBO group
    tools.parse_line(fin.readline(),r"Seq    J    NP  n    T        E-Egs     Eabs      Error    r\(p\)     r\(n\)    r\(mass\)")
    for state in state_list:
        # read line
        line=fin.readline().strip()
        if (verbose):
            print(">>> Radii:",line)

        # read state data
        match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<Ex>\S+)\s+(?P<Eabs>\S+)\s+(?P<error>\S+)\s+(?P<radii>.*)")
        Ex = float(match.group("Ex"))  # do not store since not uniquely defined when runs are combined
        error = float(match.group("error"))  # do not store since not uniquely defined when runs are combined
        radii = list(map(float,match.group("radii").split()))
        ##(state.tbo["rp"],state.tbo["rn"],state.tbo["r"]) = radii
        (results.tbo[state.qn]["rp"],results.tbo[state.qn]["rn"],results.tbo[state.qn]["r"]) = radii



    tools.parse_line(fin.readline(),r"")

    # additional TBO group
##    tools.parse_line(fin.readline(),r"Seq    J    NP  n    T      additional Two-Body observables",strict=False)
##    if (verbose):
##        print("TBO operators:",tbo_operators)
##    for state in state_list:
##        # read line
##        line=fin.readline().strip()
##        if (verbose):
##            print(">>> Additional TBO:",line)
##
##        # read state data
##        match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<TBO>.*)")
##        tbo_values = list(map(float,match.group("TBO").split()))
##        for i in range(len(tbo_values)):
##            ##state.tbo[tbo_operators[i]] = tbo_values[i]
##            results.tbo[state.qn][tbo_operators[i]] = tbo_values[i]
##    tools.parse_line(fin.readline(),r"")

    # magnetic moments group
    line=fin.readline().strip()
    if (line.split()[0] == "No"):
        # must allow for group being replaced by warning message when M_j = 0
        tools.parse_line(line,r"No magnetic moments because M_j = 0")
        ## state["moments"]["M1EM"]= None
        ## state["moments"]["M1"]= None
    else:
        tools.parse_line(line,r"Seq    J    NP  n    T   M1 moment   Dl\(pro,neu\)       Ds\(pro,neu\)     sum\(Dl\+Ds\)=J")
        for state in state_list:
            # read line
            line=fin.readline().strip()
            if (verbose):
                print(">>> Magnetic moments:",line)

            # read state data
            match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<moments>.*)")
            moment_values = list(map(float,match.group("moments").split()))
            results.moments[(state.qn,"M1EM")] = moment_values[0]
            results.moments[(state.qn,"M1")] = moment_values[1:5]
    tools.parse_line(fin.readline(),r"")

    # one-body radii and E2 group
    tools.parse_line(fin.readline(),r"Seq    J    NP  n    T     R2-OneBody\(pro,neu\)     E2\(pro,neu\)moment")
    for state in state_list:
        # read line
        line=fin.readline().strip()
        if (verbose):
            print(">>> One-body radii & E2 moments:",line)

        # read state data
        match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<moments>.*)")
        moment_values = list(map(float,match.group("moments").split()))
        state.obo["r2_ob"] = moment_values[:2]
        if (len(moment_values) == 4):
            # E2 moments defined
            results.moments[(state.qn,"E2")] = moment_values[2:]
    tools.parse_line(fin.readline(),r"")

    ################################################################
    # accumulate states to total dictionary by (J,g,n)
    ################################################################

    for state in state_list:
        qn = state.qn
        results.states[qn] = state

    ################################################################
    # access transition data
    ################################################################

    # check for transition section header
    line=fin.readline().strip()
    if (line != "Summary of transition properties"):
        # if no transition properties to read, we are done...
        return
    tools.parse_line(fin.readline(),r"")


    # read each transition group
    done_with_transitions = False
    while (not done_with_transitions):

        # check block header line
        line=fin.readline().strip()
        if (line != "Transitions to reference state  with J, NP, n, T, Eabs"):
            done_with_transitions = True
            continue

        # read final state labels
        match = tools.parse_line(fin.readline(),r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<Eabs>\S+)")
        seqf = int(match.group("seq"))

        # prime the following subblock reading seqence
        line=fin.readline().strip()

        # Gamow-Teller block
        tools.parse_line(line,r"Seq    J    NP  n    T      Eabs      Seq\(Z\+1,N\-1\)  Seq\(Z\-1,N\+1\)  Ref\(Z\+1,N\-1\)  Ref\(Z\-1,N\+1\)")
        line=fin.readline().strip()
        while (re.match(r"\d+ ",line) is not None):
            # process data line
            match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<me>.*)")
            seqi = int(match.group("seq"))
            me_values = list(map(float,match.group("me").split()))
            results.transitions[(state_lookup[seqf],state_lookup[seqi],"GT",M)] = me_values[0:2]

            # read ahead
            line=fin.readline().strip()

        # M1 block
        tools.parse_line(line,r"Seq    J    NP  n    T      Eabs       Dl\(pro,neu\)       Ds\(pro,neu\)      B\(M1: \-\>ref; ref\-\>\)")
        line=fin.readline().strip()
        while (re.match(r"\d+ ",line) is not None):
            # process data line
            match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<Eabs>\S+)\s+(?P<me>.*)")
            seqi = int(match.group("seq"))
            me_values = list(map(float,match.group("me").split()))
            results.transitions[(state_lookup[seqf],state_lookup[seqi],"M1",M)] = me_values[0:4]

            # read ahead
            line=fin.readline().strip()

        # E2 block
        tools.parse_line(line,r"Seq    J    NP  n    T      Eabs       E2\(pro,neu\)      B\(E2: \-\>ref; ref\-\>\)")
        line=fin.readline().strip()
        while (re.match(r"\d+ ",line) is not None):
            # process data line
            match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<Eabs>\S+)\s+(?P<me>.*)")
            seqi = int(match.group("seq"))
            me_values = list(map(float,match.group("me").split()))
            results.transitions[(state_lookup[seqf],state_lookup[seqi],"E2",M)] = me_values[0:2]

            # read ahead
            line=fin.readline().strip()

    # package results
    mesh_data = [results]
    return mesh_data

# register parser
input.register_data_format("mfdn_v14b05",parser)

################################################################
# test code
################################################################

if (__name__ == "__main__"):
    pass
