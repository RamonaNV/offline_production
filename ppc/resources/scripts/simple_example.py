#!/usr/bin/env python

import os, sys

if(os.getenv("I3_BUILD") == None):
    print("I3_BUILD not set.")
    sys.exit()

from os.path import expandvars

from I3Tray import *
from icecube import dataclasses, dataio, ppc

def particle(f):
    p = dataclasses.I3Particle(dataclasses.I3Position(0,0,0),
                               dataclasses.I3Direction(0,0,1), 0,
                               dataclasses.I3Particle.Cascade, 0)
    p.type = dataclasses.I3Particle.ParticleType.EMinus
    p.location_type = dataclasses.I3Particle.LocationType.InIce
    p.energy = 1.e5*I3Units.GeV
    t = dataclasses.I3MCTree()
    t.add_primary(p)
    f["particle"] = t

tray = I3Tray()

gcdfile=expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz")
tray.AddModule("I3InfiniteSource", "muxer")(
#   ("Prefix", gcdfile),
    ("Stream", icetray.I3Frame.DAQ)
    )

if(True):  # ppc part
    os.putenv("OGPU", "1")
    os.putenv("PPCTABLESDIR", expandvars("$I3_BUILD/ppc/resources/ice"))

    if(True):  # flasher simulation
        tray.AddModule("i3ppc", "ppc")(
            ("nph", 1.e9),
            ("fla", OMKey(63, 20))
            )
    else:      # particle simulation
        tray.AddModule(particle, "particle", Streams=[icetray.I3Frame.DAQ])
        tray.AddModule("i3ppc", "ppc")

tray.AddModule("I3Writer", "writer")(
    ("streams", [icetray.I3Frame.DAQ]),
    ("filename", "out.i3")
    )

tray.Execute(10)

del tray
