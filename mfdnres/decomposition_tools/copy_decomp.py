"""copy_decomp.py -- Parse and rewrite decomp data file (testbed code).

    Usage: copy_decomp.py INFILE OUTFILE

    Mark A. Caprio
    University of Notre Dame

    + 08/24/26 (mac): Created.
"""

import sys

import mfdnres.decomposition_io

################################################################
# main
################################################################

if (__name__ == "__main__"):
    
    in_filename = sys.argv[1]
    out_filename = sys.argv[2]

    # read decomp file
    print("Reading {}...".format(in_filename))
    decomp_data = mfdnres.decomposition_io.parse_decomp_file(in_filename)
    
    # write decomp file
    print("Writing {}...".format(out_filename))
    lines = mfdnres.decomposition_io.generate_decomp_file(decomp_data, header_comment_lines=["copy_decomp"])
    output_str = "\n".join(lines) + "\n"
    data_file = open(out_filename, "w")
    data_file.write(output_str)
    data_file.close()
