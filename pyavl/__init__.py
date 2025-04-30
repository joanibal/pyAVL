import importlib.metadata

try:
    # __package__ allows for the case where __name__ is "__main__"
    __version__ = importlib.metadata.version(__package__ or __name__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"

print('Warning: pyAVL is deprecated and will be removed soon, use OptVL instead (pip install optvl). '\
      'See https://joanibal.github.io/OptVL/pyavl_changes/ for an overview of the changes.')
from .pyAVL import AVLSolver

try:
    from .om_wrapper import AVLGroup, AVLMeshReader
except ImportError:
    # if openmdao is not installed, then we can't use the wrapper
    pass
