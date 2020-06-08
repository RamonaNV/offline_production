#!/usr/bin/env python

# import required icecube-related stuff
from icecube import icetray, dataclasses, simclasses, dataio
from icecube.icetray import I3Units
from I3Tray import I3Tray

# command line options required to configure the simulation
from argparse import ArgumentParser
from os.path import expandvars

#./photon_syst.py -s 1 -r 1 -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_IC86_Merged.i3.gz -o /data/user/nwandkowsky/data/test/NuE.i3.zst --flavor NuE -n 10 --gamma 1.5 --energy_min 5 --energy_max 100000 --usegpu False

icetray.logging.set_level("INFO")

usage = "usage: %prog [options] "
#outputfile"
parser = ArgumentParser()
parser.add_argument("-g", "--gcd",type=str,
        default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_IC86_Merged.i3.gz",
        dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_argument("-s", "--seed",type=int,default=12345, dest="SEED", help="Initial seed for the random number generator")
parser.add_argument("-o","--output", default="generator_output.i3", dest="OUTPUT", help="output file name")
parser.add_argument("--domeff", type=float, default=1.17, dest="DOMEFF", help="dom efficiency")
parser.add_argument("--domos", type=int, default=1, dest="DOMOS", help="dom oversizing")
parser.add_argument("--mcpeseriesname", default="I3MCPESeriesMap", dest="MCPES", help="Name of I3MCPESeriesMap in frame")
parser.add_argument("--disable-rel-dom-pos", default=True, action="store_false",dest="relativedompos", help="PhotonPositionsAreRelative")
parser.add_argument("--runnumber", type=int, default=0, dest="RUNNUMBER", help="Run number")
parser.add_argument("--holeice",
        default="/cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/metaprojects/simulation/V06-01-02/ice-models/resources/models/angsens/as.h2-50cm", 
        dest="HOLEICE", help="Holeice file")
parser.add_argument("--icemodellocation",
        default="/cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/metaprojects/simulation/V06-01-02/ice-models/resources/models/spice_3.2.1", 
        dest="icemodellocation", help="Holeice file")
parser.add_argument("-i","--input", default="generator_output.i3", dest="INPUT", help="output file name")
# parse cmd line args, bail out if anything is not understood
args = parser.parse_args()

outfile = args.OUTPUT

print("Command line options parsed...")

import os, sys

tray = I3Tray()

# import phys_services which includes rng
from icecube import phys_services

print("RNG being initiated...")
# set up a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = args.SEED,
    nstreams = 100000000,
    streamnum = args.RUNNUMBER)
    
tray.context['I3RandomService'] = randomService

tray.AddModule('I3Reader', 'reader',  filenamelist = [args.GCDFILE, args.INPUT])

tray.Add('Delete','initial_clean',Keys=['I3MCPESeriesMap'])

holeiceparameterization=expandvars(args.HOLEICE)

print("mcpes:",args.MCPES)

tray.AddModule("I3GeometryDecomposer", "decomposeGeometry")

from icecube import clsim
tray.AddSegment(clsim.I3CLSimMakeHitsFromPhotons,"makePhotons",
            PhotonSeriesName="PhotonSeriesMap",
            MCPESeriesName=args.MCPES,
            RandomService='I3RandomService',
            IceModelLocation=args.icemodellocation,
            DOMOversizeFactor=args.DOMOS,
            GCDFile=args.GCDFILE,
            UnshadowedFraction=args.DOMEFF,
            PhotonPositionsAreRelative=args.relativedompos,
            HoleIceParameterization=args.HOLEICE
        )

# Delete all MCPEs we're not operating on
def empty_mcpe(frame):
    entries = 0
    for k in frame.keys():
        if isinstance(frame[k], simclasses.I3MCPESeriesMap):
            entries = entries + len(frame[k])
    return entries>0
tray.AddModule(empty_mcpe, Streams=[icetray.I3Frame.DAQ])

# clean up unnecessary keys
#tray.Add('Delete', Keys=['PhotonSeriesMap', 'I3MCTree', 'I3MCTree_sliced', 'MMCTrackList'])
tray.Add('Delete', Keys=['I3MCTree_sliced'])

# dump content of frame to see if relevant frame objects are there
#tray.AddModule('Dump', "dumpy")

# 5. write some output data
tray.AddModule('I3Writer', 'writer',
    Streams=[icetray.I3Frame.Stream('S'), icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    filename=args.OUTPUT)

# final clean-up and execution
tray.AddModule('TrashCan', 'YesWeCan')
print("Executing...")
tray.Execute()
print("Finish!")
tray.Finish()

os.system("mv %s %s" % (args.OUTPUT,args.INPUT))

del tray
