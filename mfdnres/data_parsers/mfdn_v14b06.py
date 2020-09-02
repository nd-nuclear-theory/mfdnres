""" mfdn_v14b06.py -- provide res file parser for MFDn version 14 beta 06

    Language: Python 3
    Mark A. Caprio
    University of Notre Dame
    5/31/15 (mac): Initiated as part of mfdn_res.py.
    6/5/15 (mac): Extract from mfdn_res.py.
    10/10/16 (mac): Update to return list containing single mesh point.
    04/27/18 (mac): Rename parameter Mj to M.
    03/01/19 (mac): Fix parser failure when transitions are missing.
    06/17/20 (pjf): Add code registration.

"""

import re

# intra-package references
from .. import (
    mfdn_results_data_v14,
    input,
    tools,
    )

def read_occupations(self,fin):
    """

      More details

         RMS radius:                      proton, neutron, matter

         Orb. Occ. :         n, l, 2j,     prot,  neut
         Shell Occ.:             1+2n+l,   prot,  neut

      Next state :    1    binding energy :  -31.17544
         J, T, Ex  :    4.00000    1.00010    0.00000
         RMS radius:                          2.10634    2.33128    2.24401

         Orb. Occ. :           0   0   1       1.847836   1.876969
         ...

         Shell Occ.:                   1       1.847836   1.876969
         ...

         Total # of Nucleons:              4.000000   6.000001



    Fields:

       Orbital occupation: n l 2j np nn
       Shell occupation: (N+1) np nn

    """


    # skip to occupations
    done = False
    while (not done):
        # read line (until blank line reached)
        line=fin.readline().strip()
        if (line[0:11] == "Shell Occ.:"):
            done = True
    tools.parse_line(fin.readline(),r"")

    # iterate through states
    for seq in range(1,999):
        # read state header
        tools.parse_line(fin.readline(),r"Next state :")
        tools.parse_line(fin.readline(),r"J, T, Ex  :")
        tools.parse_line(fin.readline(),r"RMS radius:")
        tools.parse_line(fin.readline(),r"")

        # read occupations
        orbital_occupations = {}
        done = False
        while (not done):
            # read line (until blank line reached)
            line=fin.readline().strip()
            if (line == ""):
                done = True
                continue
            # extract and store data from line
            match = tools.parse_line(line,r"Orb. Occ. :\s+(?P<n>\S+)\s+(?P<l>\S+)\s+(?P<jj>\S+)\s+(?P<np>\S+)\s+(?P<nn>\S+)")
            orbital_nlj = (int(match.group("n")),int(match.group("l")),int(match.group("j"))/2)
            orbital_occ = (float(match.group("np")),float(match.group("nn")))
            orbital_occupations[orbital_nlj] = orbital_occ

    # TODO: continue

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
    tools.parse_line(fin.readline(),r"2-body interaction data in H2full format")
    line = fin.readline().strip()
    if (line.split()[0] == "hbar-omega"):
        match = tools.parse_line(line,r"hbar-omega\s+(?P<hw>[\-0-9.]+)\s+nucleon mass\s+(?P<mN>[\-0-9.]+)")
        results.params["hw"] = float(match.group("hw"))
        line = fin.readline()
    else:
        results.params["hw"] = 0.
    tools.parse_line(line,r"TBME binary file")  # Hamiltonian
    tools.parse_line(fin.readline(),r"")
    tools.parse_line(fin.readline(),r"INPUT: read in TBME files")
    tools.parse_line(fin.readline(),r"TBME binary file  tbme-rrel2.bin")

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
        # ad hoc fix to mac's squared operator names through run0365
        if (operator_name in ["L", "Sp", "Sn", "S", "J"]):
            operator_name += "2"
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
    tools.parse_line(fin.readline(),r"Seq    J    NP  n    T      additional Two-Body observables")
    if (verbose):
        print("TBO operators:",tbo_operators)
    for state in state_list:
        # read line
        line=fin.readline().strip()
        if (verbose):
            print(">>> Additional TBO:",line)

        # read state data
        match = tools.parse_line(line,r"(?P<seq>\S+)\s+(?P<J>\S+)\s+(?P<g>\S+)\s+(?P<n>\S+)\s+(?P<T>\S+)\s+(?P<TBO>.*)")
        tbo_values = list(map(float,match.group("TBO").split()))
        for i in range(len(tbo_values)):
            ##state.tbo[tbo_operators[i]] = tbo_values[i]
            results.tbo[state.qn][tbo_operators[i]] = tbo_values[i]
    tools.parse_line(fin.readline(),r"")

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

    # one-body radii and E2 group
    line=fin.readline().strip()
    if (line == ""):
        line=fin.readline().strip()
    if (line.split()[0] == "No"):
        # must allow for group being replaced by warning message when not HO basis
        tools.parse_line(line,r"No R2OB nor E2 moments because not HO basis")
    else:
        tools.parse_line(line,r"Seq    J    NP  n    T     R2-OneBody\(pro,neu\)     E2\(pro,neu\)moment")
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
    have_transitions = (line == "Summary of transition properties")

    # read transition data
    if (have_transitions):
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
input.register_data_format("mfdn_v14b06",parser)
input.register_code_name('mfdn', 'mfdn_v14b06')

################################################################
# test code
################################################################

if (__name__ == "__main__"):
    pass
