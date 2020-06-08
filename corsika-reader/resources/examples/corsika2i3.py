#!/usr/bin/env python

"""
Read showers from CORSIKA files and write to I3 files
"""

from os.path import expandvars
from argparse import ArgumentParser
parser = ArgumentParser(description=__doc__)
parser.add_argument('-i', '--infile', nargs='*', help="Input CORSIKA files",
    default=[expandvars('$I3_TESTDATA/DAT010000')])
parser.add_argument('-o', '--outfile', help="Output I3 file",
    default="DAT010000.i3.bz2")
parser.add_argument('-g', '--gcdfile', default="", help="GCD file to prepend")
parser.add_argument('-n', '--nevents', default=1, type=int, help="Number of simulated events per file")
opts = parser.parse_args()

from icecube import icetray, dataclasses, dataio, phys_services, corsika_reader
from I3Tray import I3Tray

tray = I3Tray()

tray.context['I3FileStager'] = dataio.get_stagers()
tray.context['I3RandomService'] = phys_services.I3GSLRandomService(42)
tray.AddModule('I3CORSIKAReader',
               FilenameList = opts.infile,
               Prefix = opts.gcdfile,
               NEvents = opts.nevents)

tray.Add("I3Writer", Streams=map(icetray.I3Frame.Stream, 'QP'),
    filename=opts.outfile)

tray.Execute()
