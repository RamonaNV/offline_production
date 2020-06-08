#!/usr/bin/env python

from optparse import OptionParser
from icecube.simprod.jcorsika import Corsika

parser = OptionParser()
corsika = Corsika(parser)
corsika.Execute()

