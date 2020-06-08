#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import os
import string
from os.path import expandvars

from I3Tray import *
from os.path import expandvars
import os
import sys

NUMEVENTS=10
DOMoversizing=1
SEED=12332
ICEMODEL= expandvars("$I3_SRC/ice-models/resources/models/spice_3.2.1")
gcdfile = expandvars("$I3_TESTDATA/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz")


from icecube import icetray, dataclasses, dataio, phys_services
from icecube import clsim

class generateEvent(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("Type", "", dataclasses.I3Particle.ParticleType.EMinus)
        self.AddParameter("Energy", "", 10.*I3Units.TeV)
        self.AddParameter("NEvents", "", 1)
        self.AddParameter("XCoord", "", 0.)
        self.AddParameter("YCoord", "", 0.)
        self.AddParameter("ZCoord", "", 0.)

        self.AddOutBox("OutBox")        

    def Configure(self):
        self.rs = self.GetParameter("I3RandomService")
        self.particleType = self.GetParameter("Type")
        self.energy = self.GetParameter("Energy")
        self.nEvents = self.GetParameter("NEvents")
        self.xCoord = self.GetParameter("XCoord")
        self.yCoord = self.GetParameter("YCoord")
        self.zCoord = self.GetParameter("ZCoord")
        
        self.eventCounter = 0

    def DAQ(self, frame):
        daughter = dataclasses.I3Particle()
        daughter.type = self.particleType
        daughter.energy = self.energy
        daughter.pos = dataclasses.I3Position(self.xCoord,self.yCoord,self.zCoord)
        daughter.dir = dataclasses.I3Direction(0.,0.,-1.)
        daughter.time = 0.
        daughter.location_type = dataclasses.I3Particle.LocationType.InIce

        primary = dataclasses.I3Particle()
        primary.type = dataclasses.I3Particle.ParticleType.NuE
        primary.energy = self.energy
        primary.pos = dataclasses.I3Position(self.xCoord,self.yCoord,self.zCoord)
        primary.dir = dataclasses.I3Direction(0.,0.,-1.)
        primary.time = 0.
        primary.location_type = dataclasses.I3Particle.LocationType.Anywhere

        mctree = dataclasses.I3MCTree()
        mctree.add_primary(primary)
        mctree.append_child(primary,daughter)

        frame["I3MCTree"] = mctree

        self.PushFrame(frame)

        self.eventCounter += 1
        if self.eventCounter==self.nEvents:
            self.RequestSuspension()

# a random number generator
randomService = phys_services.I3GSLRandomService(SEED)

tray = I3Tray()

# this is how you can dump some of the simulation timings&statistics to an XML file:
summary = dataclasses.I3MapStringDouble()
tray.context['I3SummaryService'] = summary

tray.AddModule("I3InfiniteSource","streams",
    Prefix = gcdfile,
    Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
    Year=2009,
    DAQTime=158100000000000000,
    RunNumber=1,
    EventID=1,
    IncrementEventID=True)

tray.AddModule(generateEvent, "generateEvent",
    I3RandomService = randomService,
    NEvents = NUMEVENTS,
    Energy = 10*I3Units.GeV,
    )


photonSeriesName = "PropagatedPhotons"
MCTreeName="I3MCTree"

usegpus = any([device.gpu for device in clsim.I3CLSimOpenCLDevice.GetAllDevices()])
tray.AddSegment(clsim.I3CLSimMakePhotons, "makeCLSimPhotons",
    UseGPUs = usegpus,
    UseOnlyDeviceNumber=0,
    UseCPUs = not usegpus,    
    GCDFile=gcdfile,
    PhotonSeriesName = None,
    MCTreeName = MCTreeName,
    RandomService = randomService,
    MCPESeriesName = "MCPESeriesMap",
    UnshadowedFraction = 1.0,
    DOMOversizeFactor=DOMoversizing,
    StopDetectedPhotons=False,
    IceModelLocation=ICEMODEL)


tray.Execute()
