#!/usr/bin/env python

import os
import sys

from I3Tray import *
from os.path import expandvars

import argparse
parser = argparse.ArgumentParser(description = "Propagates photons with ppc.")
parser.add_argument('gpu', type = int) 
parser.add_argument('nevents', type = int) 
parser.add_argument('seed', type = int, default = 42) 
parser.add_argument('infile') 
parser.add_argument('outfile') 
args = parser.parse_args()

gcdfile = expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55380_corrected.i3.gz")
os.putenv("PPCTABLESDIR", expandvars("$I3_BUILD/ppc/resources/ice"))
os.putenv("OGPU", "1")  # makes sure only GPUs are used (with OpenCL version)

from icecube import icetray, dataclasses, dataio, phys_services, sim_services, ppc

tray = I3Tray()

tray.context["I3RandomService"] = phys_services.I3GSLRandomService(args.seed)
tray.AddModule("I3Reader", FileNameList = [gcdfile, args.infile])

tray.AddModule("i3ppc", gpu = args.gpu)

tray.AddModule("I3Writer", 
               streams = [icetray.I3Frame.DAQ],
               filename = args.outfile)

if(args.nevents):
    tray.Execute(args.nevents + 3)
else:
    tray.Execute()

