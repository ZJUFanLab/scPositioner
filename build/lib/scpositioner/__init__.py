import sys
from importlib.metadata import version

from . import plot as pl
from . import tools as tl
from . import analysis as al

# has to be done at the end, after everything has been imported
sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pl", "al"]})

__version__ = version("scpositioner")

__all__ = ["__version__", "tl", "pl", "al"]
