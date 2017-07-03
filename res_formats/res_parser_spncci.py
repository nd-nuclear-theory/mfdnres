"""spncii Parser Revision 1
    Julie Butler
    June 21, 2017
    Python3

    The controlling file for the spncci results file
    parser.  The parser used here is the method
    make_dict_spncci, stored in the file make_dict.py
    (it is imported below).  Currently this file is set up
    to parse the file whos name is provided below and make
    a graph and table of energy versus hw.  This can be
    expanded later to include more analysis.

    Once the format of the results file is on a more
    finalized form, this should be rearranged into a
    class structure for ease of access to variables.

TODO:

   Needed fields:

     need to count (J,g) ->  num e-states
       can extract by tallying from [Energies]
       but maybe SpNCCI should just output table of num actual eigenstates
         by (J,g) as min::(dimension,num_eigenstates)



     params: key-value pairs from -- Space, Interaction, Calculation
     (i.e., hw) (if not counting run)

     basis: SpJ listing, BabySpNCCI listing

     For actual calculation runs:

         energies

         state.amplitudes [from Decompositions: Nex]
           => state.decomposition_Nex

         state.decompositin_baby_spncci  [from Decompositions: BabySpNCCI]

         [rms radii] -- from diagonals of the r^2 observable
            can store in tbo until split out to more appropriate name?
            no, maybe treat as a transition observable from the get-go

         [quadrupole transitions] -- actually, don't use the MFDn
         transitions structure here, either...

observables matrices -- store as matrices between
           (observable_filename,(J,g)_final,(J,g)_initial)  [from Observables::filenames]

         Modify accessors accordingly for spncci!

"""


import mfdnres.res
from numpy import array_split as array_split

# Imports the parser from the make_dict file
from mfdnres.make_dict import make_dict_spncii


def res_parser_spncci(self, fin, verbose):
    """
        Arguments:
            self: an instance of MFDnRunData.
            file_name: a string.  The filename of the results
            file that needs to be analyzed.
        Returned:
            results_dict: a dictionary.  What is returned
            from make_dic_spncci.

        This method currently takes in a file names and
        converts the file contents into a list of strings,
        where each string is a line from the file.  It then
        takes the list and passes it to the parser
        make_dict_spncci.  The dictionary returned from the
        parser is then printed and returned by the method.
        This method will be expanded and possibly incoporated
        into a class structure when the format of the results
        files is finalized.
    """
    # Stores the results from inputting the file
    results = []

    # Inputs the file passed as an argument and stores the
    # parsed file contents in results
   
    for row in fin:
        results.append(row)

    # Passes results to the parser make_dict_spncci and
    # stores the returned dictionary and list
    results_dict, order = make_dict_spncii(results)

    hw_listing = results_dict['Mesh']['hw']
    j_listing = results_dict['Branching']['J']
    print(j_listing)
    num_j = len(j_listing)

    space = results_dict['Space']
    interaction = results_dict['Interaction']
    mesh = results_dict['Mesh']
    for key,value in space.items():
        self.params[key] = value
    for key,value in interaction.items():
        self.params [key] = value
    for key,value in mesh.items():
        self.params[key] = value

    spj_listing = results_dict['SpJ (listing)']
    for x in spj_listing:
        self.spj_listing.append((float(x[1]), int(x[2])))

    baby_spncci_listing = results_dict['BabySpNCCI (listing)']
    for x in baby_spncci_listing:
        subspace_index = int(x[0])
        irrep_family_index = int(x[1])
        Nsigmaex = int(x[2])
        sigma_N = float(x[3])
        sigma_lambda = int(x[4])
        sigma_mu = int(x[5])
        sp = float(x[6])
        sn = float(x[7])
        s = float(x[8])
        nex = int(x[9])
        omega_N = float(x[10])
        omega_lambda = int(x[11])
        omega_mu = int(x[12])
        gamma_max = int(x[13])
        upsilon_max = int(x[14])
        dim = int(x[15])
        sigma = (sigma_N, sigma_lambda, sigma_mu)
        omega = (omega_N, omega_lambda, omega_mu)
        spin = (sp, sn, s)
        self.dimensions_by_omega[omega] = self.dimensions_by_omega.setdefault(omega,0) + dim
        self.baby_spncci_listing.append([sigma, omega, spin])

    state_list = []
    state_lookup = {}
    print(hw_listing)
    if len(hw_listing) > 0:
        for hw in hw_listing:
            energy = results_dict[float(hw)]['Energies'] 
            decomp_nex = results_dict[float(hw)]['Decompositions: Nex']
            decomp_baby_spncci = results_dict[float(hw)]['Decompositions: BabySpNCCI']
            observables = results_dict[float(hw)]['Observables']
            print(hw)
            state = mfdnres.res.SpNCCIStateData(hw)
            state_list.append (state)
            state_lookup[hw] = state
            for x in energy:
                qn = (float(x[0]), float(x[1]), float(x[2]))
                key = (float(hw), qn)
                E = float(x[3])
                self.properties[key] = {'J': float(x[0]), 'gex': float(x[1]), 'i': float(x[2])}
                self.energies[key] = E
            decomp_nex_split = array_split(decomp_nex, num_j)
            for i in range(0, num_j):
                state.decomposition_nex[j_listing[i]] = decomp_nex_split[i].tolist()
            decomp_baby_spncci_split = array_split(decomp_baby_spncci, num_j)
            for i in range(0, num_j):
                 state.decomposition_baby_spncci[j_listing[i]] = decomp_baby_spncci_split[i].tolist()
            observables_split = array_split(observables, num_j)
            op = str(results_dict['Observables']['filenames'])
            for i in range(0, num_j):
                 header_info = observables_split[i][0]
                 matrix = observables_split[i][1:].tolist()
                 tup = (op, (float(header_info[4]), float(header_info[5])),
                        (float(header_info[2]), float(header_info[3])))
                 state.observables[tup] = matrix


    for state in state_list:
        self.states[state.hw] = state
    
    
 
# Register the parser
mfdnres.res.register_res_format('spncci', res_parser_spncci)
#res = res_parser_spncci('type_specimens/runmac0415-Z3-N3-Nsigmamax02-Nmax02.res')
