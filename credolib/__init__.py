"""credolib

Packagified variant of the Python v3 version of the CREDO data processing library.

This is the package-ified version of credo_processing_library3.py.
It contains utilities (``macros``) to be used with the
IPython notebook interface of the CREDO data reduction procedure. For
the sake of usability, functions accept only a limited range of
arguments, most of the data transfer is done using the IPython user
namespace, i.e. by variables global across the cells of the notebook.
"""
from . import atsas
from . import calculation
from . import initialization
from . import io
from . import persistence
from . import plotting
from . import procedures
from . import qualitycontrol
from . import utils

from .atsas import *
from .calculation import *
from .initialization import *
from .io import *
from .persistence import *
from .plotting import *
from .procedures import *
from .qualitycontrol import *
from .utils import *

__all__=sorted(atsas.__all__+ calculation.__all__ + initialization.__all__ + io.__all__ + persistence.__all__ + plotting.__all__ +
               procedures.__all__ + qualitycontrol.__all__ + utils.__all__)
