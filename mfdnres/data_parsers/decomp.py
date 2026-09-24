""" decomp

    Provides parser for decomposition results files.

    Language: Python 3

    Mark A. Caprio
    University of Notre Dame

    09/18/26 (mac): Created.
"""


from .. import (
    decomposition_io,
    input,
    mfdn_results_data,
    )

from ..mfdn_results_data import (
    MFDnResultsData,  # for typing
)


################################################################
# parser
################################################################

def parser(in_file, verbose):
    """ Parse full results file.

    Arguments:
        in_file (stream): input file stream (already opened by caller)
        verbose (bool,optional): enable verbose output
    """

    # read decomposition data
    decomp_data = decomposition_io.parse_decomp_file(in_file)

    # construct results data object
    results = mfdn_results_data.MFDnResultsData()
    decomposition_type = results.params["decomposition_type"]
    qn = results.params["decomposition_state"]
    results.mfdn_level_lanczos_decomposition_data = {
        decomposition_type: {
            qn: decomp_data
        }
    }

    
    # package results
    mesh_data = [results]
    return mesh_data


# register the parser
input.register_data_format('decomp', parser)
input.register_code_name('decomp', 'decomp')


if (__name__=="__main__"):
    pass
