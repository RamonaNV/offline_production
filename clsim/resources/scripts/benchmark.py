#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from argparse import ArgumentParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = ArgumentParser(usage)
parser.add_argument("-n", "--numevents", type=int, default=1,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_argument("-s", "--seed",type=int,default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_argument("-r", "--runnumber", type=int, default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_argument("-x", "--xmlfile", default=None,
                  dest="JSONFILE", help="Write statistics to JSONFILE")
parser.add_argument("--oversize", default=1,
                  dest="OVERSIZE", help="DOM oversize factor")
parser.add_argument("--energy", default=1e3, type=float,
                  dest="ENERGY", help="Particle energy in GeV")
parser.add_argument("--type", default="EMinus",
                  dest="PARTICLE_TYPE", help="Particle type")
parser.add_argument("--icemodel", default=expandvars("$I3_BUILD/ice-models/resources/models/spice_lea"),
                  dest="ICEMODEL", help="A clsim ice model file/directory (ice models *will* affect performance metrics, always compare using the same model!)")
parser.add_argument("--unweighted-photons", action="store_true",
                  help="Propagate all Cherenkov photons. This is ~13x slower than downsampling first.")
parser.add_argument("--cable-position", help='explicitly simulate cable shadow in given position',
                  choices=('cable_shadow','led7'))

group = parser.add_mutually_exclusive_group()
group.add_argument("--minimal-gcd",  action="store_true", default=False,
                  dest="MINIMALGCD", help="generate a trivial GCD from scratch with only 24 DOMs. There are fewer collision checks, so usually things are faster, but unrealistic.")
group.add_argument("-g", "--gcd-file",
                  default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withStdNoise.i3.gz", dest="GCDFILE")

parser.add_argument("-d", "--device", type=int, default=None,
                  dest="DEVICE", help="device number")
parser.add_argument("--use-cpu",  action="store_true", default=False,
                  dest="USECPU", help="simulate using CPU instead of GPU")

# parse cmd line args, bail out if anything is not understood
options = parser.parse_args()

if options.DEVICE is not None:
    print(" ")
    print(" ** DEVICE selected using the \"-d/--device\" command line option. Only do this if you know what you are doing!")
    print(" ** You should be using the CUDA_VISIBLE_DEVICES and/or GPU_DEVICE_ORDINAL environment variables instead.")

if options.MINIMALGCD:
    parser.error("--minimal-gcd does not work with I3CLSimClientModule; it needs an external GCD file")
    print(" ")
    print(" ** You chose to not use a standard IceCube GCD file but instead to create a trivial geometry from scratch.")
    print(" ** This geometry only has 24 DOMs, so there are fewer collision checks.")
    print(" ** This usually means propagation is faster, but unrealistic. Might differ from GPU type to GPU type.")


from I3Tray import *
import os
import sys
import math
import numpy

from icecube import icetray, dataclasses, dataio, phys_services, sim_services, simclasses, clsim

icetray.logging.rotating_files('logtrace')


#icetray.I3Logger.global_logger.set_levelicetray.I3LogLevel.LOG_INFO)
#icetray.I3Logger.global_logger.set_level(icetray.I3LogLevel.LOG_WARN)
icetray.I3Logger.global_logger.set_level(icetray.I3LogLevel.LOG_TRACE)

radius = 120.*I3Units.m

omPos = numpy.array(
    [[ 0.,  1., 0.],
     [ 1.,  1., 0.],
     [ 1.,  0., 0.],
     [ 1., -1., 0.],
     [ 0., -1., 0.],
     [-1., -1., 0.],
     [-1.,  0., 0.],
     [-1.,  1., 0.]]
    )
# normalize and scale
omPos = (omPos.T/numpy.sqrt(numpy.sum(omPos**2, 1))).T * radius

omPosLower = numpy.array(omPos)
omPosLower.T[2] = omPosLower.T[2] - radius
omPosUpper = numpy.array(omPos)
omPosUpper.T[2] = omPosUpper.T[2] + radius

omPositions = numpy.concatenate((omPosUpper, omPos, omPosLower), axis=0)


omKeys = [
    icetray.OMKey(1,1),
    icetray.OMKey(2,1),
    icetray.OMKey(3,1),
    icetray.OMKey(4,1),
    icetray.OMKey(5,1),
    icetray.OMKey(6,1),
    icetray.OMKey(7,1),
    icetray.OMKey(8,1),

    icetray.OMKey(1,2),
    icetray.OMKey(2,2),
    icetray.OMKey(3,2),
    icetray.OMKey(4,2),
    icetray.OMKey(5,2),
    icetray.OMKey(6,2),
    icetray.OMKey(7,2),
    icetray.OMKey(8,2),

    icetray.OMKey(1,3),
    icetray.OMKey(2,3),
    icetray.OMKey(3,3),
    icetray.OMKey(4,3),
    icetray.OMKey(5,3),
    icetray.OMKey(6,3),
    icetray.OMKey(7,3),
    icetray.OMKey(8,3),

]


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


class injectFakeGCD(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("OMKeys",      "", [])
        self.AddParameter("OMPositions", "", [])
        self.AddParameter("XCoord", "", 0.)
        self.AddParameter("YCoord", "", 0.)
        self.AddParameter("ZCoord", "", 0.)

        self.AddOutBox("OutBox")        

    def Configure(self):
        self.omkeys = self.GetParameter("OMKeys")
        self.ompositions = self.GetParameter("OMPositions")
        self.xCoord = self.GetParameter("XCoord")
        self.yCoord = self.GetParameter("YCoord")
        self.zCoord = self.GetParameter("ZCoord")

        self.has_been_injected = False

    def DAQ(self, frame):
        # only inject it once
        if self.has_been_injected:
            self.PushFrame(frame)
            return
        self.has_been_injected = True

        geometry = dataclasses.I3Geometry()
        calibration = dataclasses.I3Calibration()
        detectorStatus = dataclasses.I3DetectorStatus()

        # fill the geometry map
        omgeomap = geometry.omgeo
        domcalmap = calibration.dom_cal
        domstatusmap = detectorStatus.dom_status

        for i, pos in enumerate(omPositions):
            shiftedPos = pos
            shiftedPos[0] += self.xCoord*I3Units.m
            shiftedPos[1] += self.yCoord*I3Units.m
            shiftedPos[2] += self.zCoord*I3Units.m

            omkey = omKeys[i]

            newomgeo = dataclasses.I3OMGeo()
            newomgeo.omtype = dataclasses.I3OMGeo.OMType.IceCube
            newomgeo.orientation = dataclasses.I3Orientation(dataclasses.I3Direction(0.,0.,-1.))
            newomgeo.position = dataclasses.I3Position(shiftedPos[0], shiftedPos[1], shiftedPos[2])
            omgeomap[omkey] = newomgeo


            newdomcal = dataclasses.I3DOMCalibration()
            newdomcal.relative_dom_eff = 1.0
            domcalmap[omkey] = newdomcal


            newdomstatus = dataclasses.I3DOMStatus()
            newdomstatus.pmt_hv = 1345.*I3Units.V # some arbitrary setting: >0 and not NaN
            domstatusmap[omkey] = newdomstatus


        # make GCD frames and fill them with objects
        Gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
        Cframe = icetray.I3Frame(icetray.I3Frame.Calibration)
        Dframe = icetray.I3Frame(icetray.I3Frame.DetectorStatus)

        Gframe["I3Geometry"] = geometry
        Cframe["I3Calibration"] = calibration
        Dframe["I3DetectorStatus"] = detectorStatus

        # push the new GCD frames
        self.PushFrame(Gframe)
        self.PushFrame(Cframe)
        self.PushFrame(Dframe)

        # push the original Q-frame
        self.PushFrame(frame)


tray = I3Tray()

summary = dataclasses.I3MapStringDouble()
tray.context['I3SummaryService'] = summary

# a random number generator
try:
    randomService = phys_services.I3SPRNGRandomService(
        seed = options.SEED,
        nstreams = 10000,
        streamnum = options.RUNNUMBER)
except AttributeError:
    randomService = phys_services.I3GSLRandomService(
        seed = options.SEED*1000000 + options.RUNNUMBER)

if options.MINIMALGCD:
    tray.AddModule("I3InfiniteSource","streams",
        Stream=icetray.I3Frame.DAQ)

    tray.AddModule(injectFakeGCD,"gcd",
        OMKeys = omKeys,
        OMPositions = omPositions,
        # XCoord = xCoord,
        # YCoord = yCoord,
        # ZCoord = zCoord,
        )
else:
    tray.AddModule("I3InfiniteSource","streams",
        Prefix = options.GCDFILE,
        Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
    Year=2009,
    DAQTime=158100000000000000,
    RunNumber=1,
    EventID=1,
    IncrementEventID=True)

tray.AddModule(generateEvent, "generateEvent",
    I3RandomService = randomService,
    NEvents = options.NUMEVENTS,
    Energy = options.ENERGY,
    Type = getattr(dataclasses.I3Particle.ParticleType, options.PARTICLE_TYPE),
    # Energy = 1000.*I3Units.TeV,
    # XCoord = xCoord,
    # YCoord = yCoord,
    # ZCoord = zCoord,
    )

MCTreeName="I3MCTree"
photonSeriesName = None

kwargs = {}
if options.cable_position:
    from icecube.clsim import GetIceCubeCableShadow
    kwargs['CableOrientation'] = GetIceCubeCableShadow.GetIceCubeCableShadow(getattr(GetIceCubeCableShadow, 'from_{}'.format(options.cable_position)))
tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
    GCDFile = options.GCDFILE,
    PhotonSeriesName = photonSeriesName,
    MCTreeName = MCTreeName,
    RandomService = randomService,
    MCPESeriesName = "MCPESeriesMap",
    UnshadowedFraction = 0.95,
    UseGPUs=not options.USECPU,
    UseCPUs=options.USECPU,
    UseOnlyDeviceNumber=options.DEVICE,
    IceModelLocation=options.ICEMODEL,
    DOMOversizeFactor=options.OVERSIZE,
    UnWeightedPhotons=options.unweighted_photons,
    **kwargs
    )

icetray.logging.set_level_for_unit('I3CLSimServer', 'INFO')
icetray.logging.set_level_for_unit('I3CLSimStepToPhotonConverterOpenCL', 'INFO')

from datetime import datetime
t0 = datetime.now()
tray.Execute()
walltime_in_execute = ((datetime.now() - t0).total_seconds())*1e9

del tray

if options.JSONFILE:
    import json
    with open(options.JSONFILE, 'w') as f:
        json.dump(dict(summary), f, indent=1)

########### this is optional and just parses the generated summary
import numpy as np
def get(key, aggregate=np.mean, default=0):
    pkey = 'I3CLSimModule_makeCLSimHits_makePhotons_clsim_'+key
    return aggregate([summary.get(k,default) for k in summary.keys() if k.startswith(pkey)])
ns_per_photon = get('AverageDeviceTimePerPhoton')
ns_per_photon_with_util = get('AverageHostTimePerPhoton')
device_util = get('DeviceUtilization')
ncalls = get('NumKernelCalls', sum)

if ncalls == 0:
    sys.stderr.write("Not enough kernel calls to estimate performance! Trying increasing the number of events.\n")
    sys.exit(1)

total_host_time = get('TotalHostTime', sum)
total_queue_time = get('TotalQueueTime', sum)

class duration(float):
    def __format__(self, format_spec):
        if self > 2e9:
            return format(self/1e9, format_spec) + ' s'
        elif self > 2e6:
            return format(self/1e6, format_spec) + ' ms'
        elif self > 2e3:
            return format(self/1e3, format_spec) + ' µs'
        else:
            return format(self/1e6, format_spec) + ' ns'

print(" ")
print("# these numbers are performance figures for the GPU:")
print("time per photon (GPU):", ns_per_photon, "ns")
print("photons per second (GPU):", 1e9/ns_per_photon, "photons per second")

print(" ")
print("# these numbers include the host utilization and are probably not meaningful for --numevents=1 (the default). You need more events to even out the startup/setup time.")
print("(avg)   time per photon (actual, including under-utilization):", ns_per_photon_with_util, "ns")
print("(avg)   photons per second (actual, including under-utilization):", 1e9/ns_per_photon_with_util, "photons per second")
print("(total) host time: {:.1f}".format(duration(total_host_time)))
print("(total) waiting time: {:.1f} ({:.3f}%)".format(duration(total_queue_time), 100.*total_queue_time/total_host_time))
print("(total) number of kernel calls: {:.0f}".format(ncalls))
print("        wallclock time: {:.1f}".format(duration(walltime_in_execute)))

print("(avg)   device utilization:", device_util*100., "%")

