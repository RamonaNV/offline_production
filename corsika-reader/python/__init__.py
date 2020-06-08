from icecube.load_pybindings import load_pybindings
try:
    load_pybindings(__name__,__path__)
    from . import converters
except ImportError:
    from icecube.icetray import load
    load('corsika-reader', False)
    del load
del load_pybindings

from .ReadCorsika import ReadCorsika

