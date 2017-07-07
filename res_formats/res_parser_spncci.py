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

    Methods:
        res_parser_spncci: Sorts the parsed SpnCCI results file
            into an instance of SpNCCIMeshPointData.  Takes as
            arguments an instance of SpNCCIMeshPointData, a 
            dictionary with comes from the method make_dict_spncci,
            from make_dict.py, and a boolean for debugging.  Returns
            nothing.  Extracts the relevant information from results_dict
            and stores it in the appropriate attributes of the 
            instance of SpNCCIMeshPointData.


"""


import mfdnres.res
from numpy import array_split as array_split

# Imports the parser from the make_dict file
from mfdnres.make_dict import make_dict_spncci


def res_parser_spncci(self, results_dict, verbose):
    """
        Arguments:
            self (instance of SpNCCIMeshPointData)
            results_dict (dictionary): Created by make_dict.make_dict_spncci. 
        Returned:
            None.

        Takes in the dictionary made by the make_dict_spncci method
        of make_dict.py as well as an instance of SpNCCIMeshPointData.
        Extract the relevant information from results_dict and stores
        it in the appropriate attribute of self.
    """
    # Determines what J-values are present in the run and how many.
    j_listing = results_dict['Branching']['J']
    num_j = len(j_listing)

    ################################################################
    # populate params
    ################################################################

    # Set values based on entries in Space and Interaction sections, plus
    # hw from Calculation section.

    # nuclide = 3 3 -- list of int
    # A = 6  -- int
    # Nsigma0 = 9.5  -- float
    # Nsigmamax = 2  -- int
    # N1v = 1  -- int
    # Nmax = 2  -- int
    # 
    # interaction = RESERVED  -- str
    # use_coulomb = 0  -- int -> bool

    space = results_dict['Space']
    interaction = results_dict['Interaction']
    for key,value in space.items():
        self.params[key] = value
    for key,value in interaction.items():
        self.params[key] = value

    # convert fields
    conversions = {
        "Nsigmamax" : int,
        "Nmax" : int
    }
    for key in conversions:
        conversion = conversions[key]
        self.params[key] = conversion(self.params[key])

    self.params["hw"] = self.hw

    # populate spj_listing
    
    # Stores the information under the heading 'SpJ (listing)'
    # in self.spj_listing
    spj_listing = results_dict['SpJ (listing)']
    for x in spj_listing:
        self.spj_listing.append((float(x[1]), int(x[2])))

    # Stores the information under the heading 'BabySpNCCI (listing)'
    # in self.baby_spncci_listing and self.dimensions_by_omega
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

    # Retrieves the information from the headings 'Energies', 'Decompositions: Nex', 
    # 'Decomposition: BabySpNCCI' and 'Observables' by the hw value supplied when
    # self was created.
    energy = results_dict[self.hw]['Energies'] 
    decomp_nex = results_dict[self.hw]['Decompositions: Nex']
    decomp_baby_spncci = results_dict[self.hw]['Decompositions: BabySpNCCI']
    observables = results_dict[self.hw]['Observables']
   
    # Fills self.energies with the data stored under the heading 'Energies'
    for x in energy:
        qn = (float(x[0]), float(x[1]), float(x[2]))
        E = float(x[3])
        self.energies[qn] = E

    # Stores the information from 'Decompositions: Nex' and 'Decompositions:
    # BabySpNCCI' in self.decompositions
    decomposition_nex = {}
    decomposition_baby_spncci = {}
    decomp_nex_split = array_split(decomp_nex, num_j)
    for i in range(0, num_j):
        decomposition_nex[j_listing[i]] = decomp_nex_split[i].tolist()
    self.decomposition['Nex'] = decomposition_nex
    decomp_baby_spncci_split = array_split(decomp_baby_spncci, num_j)
    for i in range(0, num_j):
         decomposition_baby_spncci[j_listing[i]] = decomp_baby_spncci_split[i].tolist()
    self.decomposition['BabySpNCCI'] = decomposition_baby_spncci

    # Stores the information under the heading 'Observables' in self.observables
    observables_split = array_split(observables, num_j)
    op = str(results_dict['Observables']['filenames'])
    for i in range(0, num_j):
        header_info = observables_split[i][0]
        if len(header_info) > 6:
            matrix = observables_split[i][1:].tolist()
            tup = (op, (float(header_info[4]), float(header_info[5])),
                (float(header_info[2]), float(header_info[3])))
            self.observables[tup] = matrix


# Register the parser
mfdnres.res.register_res_format('spncci', res_parser_spncci)
#res = res_parser_spncci('type_specimens/runmac0415-Z3-N3-Nsigmamax02-Nmax02.res')
