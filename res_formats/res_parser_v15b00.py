"""
    MFDn Version 15 Parser Revision 13
    Julie Butler
    June 21, 2017
    Python3

    Parses the results file from MFDN Version 15.  This revision
    has been designed to fit into and work with the already
    existing mfdnres structure.  Instead of using regular expressions
    like the previous version of the parser, this parser stores
    all the data in a dictionary, with the key being the section label
    and the structure of the value depending on the type of data in the
    section. (See the documentation in make_dict for more information).

    This file is currently setup with test code at the bottom for debuggin
    purposes.
"""


# File needed from the mfdnres code (Will need to be redone with appropriate path at the end)
import mfdnres.res

# The location of the parser
from mfdnres.make_dict import make_dict_mfdn15


def res_parser_v15b00(self, fin, verbose):
    """
        Arguments:
            self: an instance of MFDnRunData.
            fin: a file pointer.  The results file to be read
            verbose: a boolean. for debugging purposes
        Returned:
            None.

        Parses a MFDn version 15 results file and extract data to be assigned to either
        members of self or as properties of instances of MFDnStateData.
    """
    # Takes all the lines from the file pointer and puts then into results as a strings
    # This makes results a list of strings
    results = []
    for row in fin:
        results.append(row)

    # Stores the values returned by the parser
    results_dict, order = make_dict_mfdn15(results)

    # Extract the data needed to fill the dictionary self.params
    self.params['hw'] = results_dict['Interaction']['hbomeg']
    self.params['nuclide'] = (results_dict['Basis']['Nprotons'], results_dict['Basis']['Nneutrons'])
    self.params['Nmax'] = results_dict['Basis']['Nmax']
    self.params['Nmin'] = results_dict['Basis']['Nmin']

    # Extract the value of the parity variable and calculates the value of
    # g from the value of parity  (parity = (-1)^g)
    parity = float(results_dict['Basis']['parity'])
    if parity == 1:
        g = 0
    elif parity == -1:
        g = 1
    elif parity == 0:
        g = 'None'
    else:
        raise ValueError('Parity is an undefined value.')

    # set up data structure to hold states on input
    state_list = []  # states within current run
    state_lookup = {}  # lookup from current run's sequence number to (J,g,n) indices

    # This section fills the energy list of MDFnRunData, creates instances of MFDnStateData
    # and assigned this instances properties from Energies and Occupation Probabilities
    # if also fills self.properties, self.energies, and self.orbital_occupations
    energies = results_dict['Energies']
    occupation_probabilities = results_dict['Occupation probabilities']
    occupation_offset = 4
    for x in energies:
        # Extracting values from the 'Energies' section and storing the
        # results in dictionaries and making instances of MFDnStateData
        # for each sequence nu,ber
        seq = x[0]
        qn = (float(x[1]), g, float(x[2]))
        T = float(x[3])
        E = float(x[4])
        state = res.MFDnStateData(qn, T, E)
        state_list.append(state)
        state_lookup[seq] = qn
        self.properties[state.qn] = {'J': x[1], 'g': g, 'n': x[2], 'T': T}
        self.energies[state.qn] = E
        # Extract the column of data from the 'Occupation probabilities'
        # section that corresponds to the current sequence number.  This
        # data is then stored in state.occupations.
        occupations = []
        index = occupation_offset + int(seq)
        for y in occupation_probabilities:
            occupations.append(y[index])
        state.occupations = occupations
        # This section deals with self.orbital_occupations.  The condition
        # should be true because there should be an equivalent number of
        # proton and neutron occupation probabilities.
        if len(occupation_probabilities) % 2 == 0:
            orb_occ = {}
            proton_index = 0
            neutron_index = int(len(occupation_probabilities)/2)
            while (
                    (proton_index < len(occupation_probabilities)/2 - 1)
                    and (neutron_index < len(occupation_probabilities))
            ):
                key = (float(occupation_probabilities[proton_index][2]),
                    float(occupation_probabilities[proton_index][3]), 
                    float(occupation_probabilities[proton_index][4])/2)
                value = (float(occupation_probabilities[proton_index][index]),
                    float(occupation_probabilities[neutron_index][index]))
                orb_occ[key] = value
                proton_index = proton_index + 1
                neutron_index = neutron_index + 1
            self.orbital_occupations[qn] = orb_occ

    # This section assigns the amplutudes property of the MFDnStateData instances
    # from the 'Oscillator quanta' section of the results file
    oscillator = results_dict['Oscillator quanta']
    for x in oscillator:
        seq = x[0]
        qn = (float(x[1]), g, float(x[2]))
        amp = x[4:]
        amp = [float(i) for i in amp]
        state_list[int(seq)-1].amplitudes = amp
        if state_lookup[seq] == qn:
            # makes sure TBO dictionary available
            self.tbo.setdefault(qn, {})
        else:
            print('Quantum numbers do not match in oscillator')

    # This section extract the radii from the results file
    radii = results_dict['Relative radii']
    for x in radii:
        # Extract the relevant data from each data entry in 'Relative radii'
        seq = x[0]
        qn = (float(x[1]), g, float(x[2]))
        radii = x[4:]
        radii = [float(i) for i in radii]
        # Check to make sure the quantum numbers of the row matches the quantum
        # number of the state that corresponds to the sequence number
        if state_lookup[seq] == qn:
            (self.tbo[qn]['rp'], self.tbo[qn]['rn'], self.tbo[qn]['r']) = radii
        else:
            print('Quantum numbers do not match in oscillator')

    # read names of additional two-body observables from the names of the
    # TBME file names
    observables = results_dict['Observables']
    tbo_operators = []
    tbme_filenames = []
    for key, value in observables.items():
        if key[0:5] == 'TBMEf':
            tbme_filenames.append(value)
    for operator_name in tbme_filenames:
        operator_name = operator_name[5:]
        operator_name = operator_name[:len(operator_name)-4]
        if not (operator_name == 'H' or operator_name == 'rrel2'):
            # ad hoc fix to mac's squared operator names through run0365
            if (operator_name in ["L", "Sp", "Sn", "S", "J"]):
                operator_name += "2"
            tbo_operators.append(operator_name)

    other_2_body_observables = results_dict['Other 2-body observables']
    for x in other_2_body_observables:
        # Extracts the relevant data from each line in 'Other 2-body
        # observables.
        seq = x[0]
        qn = (float(x[1]), g, float(x[2]))
        other = x[4:]
        other = [float(i) for i in other]
        # Check to make sure the quantum numbers of the row matches the quantum
        # number of the state that corresponds to the sequence number
        if state_lookup[seq] == qn:
            # matches the extracted observable to the variable names extracted from
            # the TBME file names (variable names are stored in tbo_operators)
            for j in range(len(other)):
                self.tbo[qn][tbo_operators[j]] = other[j]
        else:
            print('Quantum numbers do not match in oscillator')

    # This section handles the label '[M1 moments]' and puts the appropriate data into
    # self.moments
    # Check to make sure the section 'M1 moments' exist and contains data
    if 'M1 moments' in results_dict and (not results_dict['M1 moments'] == 'None'):
        m1_moments = results_dict['M1 moments']
        for x in m1_moments:
            # Extracts the relevant data, which is characterized as either 'M1EM' or 'M1'
            seq = x[0]
            qn = (float(x[1]), g, float(x[2]))
            m1EM = float(x[4])
            m1 = x[5:]
            m1 = [float(i) for i in m1]
        # Check to make sure the quantum numbers of the row matches the quantum
        # number of the state that corresponds to the sequence number
            if state_lookup[seq] == qn:
                self.moments[(qn, 'M1EM')] = m1EM
                self.moments[(qn, 'M1')] = m1
            else:
                print('Quantum numbers do not match in oscillator')
    else:
        print('No M1 moments data section in the file, skipping to the next section.')

    # This section handles the label '[E2 Moments]' and puts the appropriate data into
    # self.moments.  It also assigns to each state an obo.
    # Check tp make sure the 'E2 moments' data section exist and contains data
    if 'E2 moments' in results_dict and (not results_dict['E2 moments'] == 'None'):
        e2_moments = results_dict['E2 moments']
        for x in e2_moments:
            # Extracts the relevant data
            seq = x[0]
            qn = (float(x[1]), g, float(x[2]))
            moments = x[4:]
            # Assigns to the state corresponding to the sequence number values for obo
            state_list[int(seq)-1].obo['r2_ob'] = moments[:2]
            if len(moments) == 4:
                # Check to make sure the quantum numbers of the row matches the quantum
                # number of the state that corresponds to the sequence number
                if state_lookup[seq] == qn:
                    self.moments[(qn, 'E2')] = moments[2:]
                else:
                    print('Quantum numbers do not match in oscillator')
    else:
        print ('No E2 moments data section in the file, skipping to the next section.')

    # Extract Mj which is stored in the section 'Basis'
    Mj = results_dict['Basis']['TwoMj']

    # This sections handles the Transitions data label and its data from the results file
    # Check to make sure the 'Transitions' section exist and contains data
    if 'Transitions' in results_dict and (not results_dict['Transitions'] == 'None'):
        transitions = results_dict['Transitions']
        for x in transitions:
            # Extracts out the relevant information from both the initial and final states
            seq_f = x[0]
            qn_f = (float(x[1]), g, float(x[2]))
            category = x[3]
            seq_i = x[4]
            qn_i = (float(x[5]), g, float(x[6]))
            rest = x[7:]
            rest = [float(i) for i in rest]
        # Check to make sure the quantum numbers of the row matches the quantum
        # number of the state that corresponds to the sequence number, for both
        # the initial and final states.  Assigns values to self.transitions based
        # on the type of transition.
            if state_lookup[seq_f] == qn_f and state_lookup[seq_i] == qn_i:
                if category == 'GT':
                    self.transitions[(qn_f, qn_i, category, Mj)] = rest[0:2]
                elif category == 'M1':
                    self.transitions[(qn_f, qn_i, category, Mj)] = rest[0:4]
                elif category == 'E2':
                    self.transitions[(qn_f, qn_i, category, Mj)] = rest[0:2]
                else:
                    raise ValueError('Unknown transition type')
            else:
                print('Quantum numbers do not match in oscillator')
    else:
        print ('No Tranistions data section in the file, skipping to the next section.')

    # Fills the self.states with state_list made in the energy section
    for state in state_list:
        qn = state.qn
        self.states[qn] = state


# Registers the parse
mfdnres.res.register_res_format('v15b00', res_parser_v15b00)

################################################################
# test code
################################################################
if __name__ == "__main__":
    # unit test is currently somewhat inaccessible, since running this
    # module direcly in its cwd means mfdnres packages cannot be found
    # unless added to PYTHONPATH
    filename = 'HeReal15/runaem0007-mfdn15-Z2-N1-JISP16-coul0-hw07.711-a_cm20-Nmax02-Mj0.5-lan2000-tol1.0e-06.res'
    data = res.MFDnRunData()
    data.read_file(filename, res_format="v15b00", verbose=True)
    print("states {}, moments {}, transitions {}".format(len(data.states), len(data.moments), len(data.transitions)))
