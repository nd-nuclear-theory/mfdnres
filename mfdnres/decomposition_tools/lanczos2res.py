"""lanczos2res.py

    Generate decomposition results files by combining (legacy) Lanczos
    alpha-beta files (.lanczos) and decomposition data files (.decomp).

    Usage:

        lanczos2res run_dir decomp_dir

    Processes all .lanczos files under run_dir, taking decomposition data from
    decomp_dir.

    Ex:
   
        % lanczos2res ${MFDNRES_RESULTS_DIR}/mcaprio/mfdn/runmac0967 ${GROUP_HOME}/data/decomposition/201112-mcaprio

        Reading decomp file /home/mcaprio/data/decomposition/201112-mcaprio/Z03-N08-Nmax02-U3SpSnS.decomp...
        Reading lanczos file /home/mcaprio/results/mcaprio/mfdn/runmac0967/results/lanczos/runmac0967-mfdn15-Z3-N8-Daejeon16-coul1-hw12.500-Nmax04-J01.5-g1-n02-U3SpSnS-dNmax02-dlan0800.lanczos...
        Writing res file /home/mcaprio/results/mcaprio/mfdn/runmac0967/results/res/mac0967-decomp-Z3-N8-Daejeon16-coul1-hw12.500-Nmax04-J01.5-g1-n02-U3SpSnS-dNmax02-dlan0800.res...
        Reading decomp file /home/mcaprio/data/decomposition/201112-mcaprio/Z03-N08-Nmax02-U3SpSnS.decomp...
        Reading lanczos file /home/mcaprio/results/mcaprio/mfdn/runmac0967/results/lanczos/runmac0967-mfdn15-Z3-N8-Daejeon16-coul1-hw10.000-Nmax08-J01.5-g1-n02-U3SpSnS-dNmax02-dlan0800.lanczos...
        Writing res file /home/mcaprio/results/mcaprio/mfdn/runmac0967/results/res/mac0967-decomp-Z3-N8-Daejeon16-coul1-hw10.000-Nmax08-J01.5-g1-n02-U3SpSnS-dNmax02-dlan0800.res...
        ...

    Mark A. Caprio
    University of Notre Dame

    + 10/05/26 (mac): Created.

"""

import argparse
import glob
import os
import os.path

import numpy as np

import mfdnres.input
import mfdnres.decomposition
import mfdnres.decomposition_io


################################################################
# parse arguments
################################################################

def parse_args():
    """Parse arguments.

    Returns:
        (argparse.Namespace) parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Combine legacy lanczos files with decomposition data to yield decomposition results files.",
        usage="%(prog)s run_dir decomp_dir\n",
        epilog="The given run directory should contain a subdirectory results/lanczos.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

    # positional arguments
    parser.add_argument("run_dir", help="Run directory")
    parser.add_argument("decomp_dir", help="Decomposition data file directory")

    # options
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress output text")
    
    return parser.parse_args()


################################################################
# main
################################################################

def main():
    
    args = parse_args()
    run_dir = args.run_dir
    decomp_dir = args.decomp_dir

    # create target directory for res files
    os.makedirs(os.path.join(run_dir, "results", "res"), exist_ok=True)

    # process lanczos files
    lanczos_filenames = glob.glob(os.path.join(run_dir, "results", "lanczos", "*.lanczos"))
    for lanczos_filename in lanczos_filenames:

        # extract info from lanczos filename
        ## print(lanczos_filename)
        lanczos_basename = os.path.basename(lanczos_filename)
        info_from_filename = mfdnres.input.parse_filename(lanczos_basename)
        ## print(info_from_filename)

        # read decomposition data
        decomp_filename_template = "Z{nuclide[0]:02d}-N{nuclide[1]:02d}-Nmax{Nmax:02d}-{decomposition_type}.decomp"
        nuclide = info_from_filename["nuclide"]
        decomposition_Nmax = info_from_filename["decomposition_Nmax"]
        decomposition_type = info_from_filename["decomposition_type"]
        decomp_filename = os.path.join(
            decomp_dir,
            decomp_filename_template.format(
                nuclide=nuclide, Nmax=decomposition_Nmax, decomposition_type=decomposition_type,
            ),
        )
        if not args.quiet:
            print("Reading decomp file {}...".format(decomp_filename))
        decomposition_data = mfdnres.decomposition_io.parse_decomp_file(decomp_filename)
        
        # read lanczos resuilts
        if not args.quiet:
            print("Reading lanczos file {}...".format(lanczos_filename))
        alpha_beta = np.loadtxt(lanczos_filename, usecols=(1, 2))
        decomposition_data["lanczos"] = alpha_beta

        # write decomposition results
        res_filename_template = "{run}-{code_name}-{descriptor}.res"
        run = info_from_filename["run"]
        code_name = "decomp"
        descriptor = info_from_filename["descriptor"]
        res_filename = os.path.join(
            run_dir, "results", "res",
            res_filename_template.format(run=run, code_name=code_name, descriptor=descriptor),
        )
        header_comment_lines = [
            "lanczos2res",
            "Run: {}".format(run),
            "Descriptor: {}".format(descriptor),
            "Decomposition data: {}".format(decomp_filename),
        ]
        lines = mfdnres.decomposition_io.generate_decomp_file(decomposition_data, header_comment_lines)
        output_str = "\n".join(lines) + "\n"
        if not args.quiet:
            print("Writing res file {}...".format(res_filename))
        res_file = open(res_filename, "w")
        res_file.write(output_str)
        res_file.close()

        
if (__name__ == "__main__"):
    
    main()
