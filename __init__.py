
# res import and retrieval
import mfdnres.res
from . import res
import mfdnres.res_formats  # register parsers

# descriptor/filename parsing
import mfdnres.descriptor
from . import descriptor
import mfdnres.filename_formats  # register parsers

# analysis control
import mfdnres.analysis
from . import analysis

# optional: user must load
## import mfdnres.band
## from . import band
