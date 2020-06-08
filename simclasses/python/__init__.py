# this file is based on dataclasses/python/__init__.py
from icecube.load_pybindings import load_pybindings
import icecube.icetray  # pull in our dependencies

icecube.icetray.load('simclasses', False)
load_pybindings(__name__, __path__)
