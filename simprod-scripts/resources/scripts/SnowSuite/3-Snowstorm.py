# -*- coding: utf-8 -*-
# Copyright (c) 2019
# Jakob van Santen <jakob.van.santen@desy.de>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
# $Id$
#
# @file hobo-multisim.py
# @version $Revision$
# @date $Date$
# @author Jakob van Santen

# Benjamin Smithers  benjamin.smithers@mavs.uta.edu
# Erik Ganster erik.ganster@rwth-aachen.de

print("        ===== Welcome to SnowSuite =====")
print("")
print("                       /\ ")
print("                  __   \/   __")
print("                  \_\_\/\/_/_/ ")
print("                    _\_\/_/_ ")
print("                   __/_/\_\__ ")
print("                  /_/ /\/\ \_\ ")
print("                       /\ ")
print("                       \/ ")
print("")
print("        ====   Running Snowstorm   ====")
print("  Some Novel Operations Build Onto Snowstorm")
print("")

# just load in argparse so the help function can be quickly accessed
from argparse import ArgumentParser
from utils import str2bool

parser = ArgumentParser()
parser.add_argument("-i","--infile",
                    type=str, required=True,
                    help="File to process")

parser.add_argument("-o","--outfile",
                    type=str, required=True,
                    help="Output file")

parser.add_argument("-g", "--gcdfile",
                    type=str, required=True,
                    help="Calibration file to use")

parser.add_argument("-c", "--config_file",
                    type=str, required=True,
                    help="YAML/JSON config file for Snwostorm")

parser.add_argument("-s", "--seed",
                    type=int, required=True,
                    help="Seed to pass to the random service")

parser.add_argument("--domoversizefactor",
                    type=float, required=True,
                    help="DOM oversize factor")

parser.add_argument("--summaryfile",
                    type=str, default="summary_snowstorm.json",
                    help="(optional) Name and path of a summary file. By default a 'summary_snowstorm.json' will be created in the current directory.")

parser.add_argument("--events-per-model",
                    type=int, required=True,
                    help="Number of frames to process before scrambling systematics")

parser.add_argument("--cpu",
                    type=str2bool, default=False,
                    nargs='?', const=True,
                    help="Call this flag when NOT using GPUs.")

parser.add_argument("--UseI3PropagatorService", dest="UseI3PropagatorService",
                    type=str2bool, default=True,
                    nargs='?', const=False,
                    help="Call this flag for NOT using I3PropagatorService.")

parser.add_argument("--log-level", dest="log_level",
                    type=str, default="WARN",
                    help="Sets the icetray logging level (ERROR, WARN, INFO, DEBUG, TRACE)")

args = parser.parse_args()

print("Called with:")
for key, v in vars(args).items():
    print("{:32s}: {}".format(key, v))
print("")

# Check for CPU processing
if args.cpu:
    import warnings
    warnings.warn("WARNING: cpu processing is really slow!")

# import time to measure runtime
import time
start_time = time.time()

print("Importing packages/modules... ", end='')
import os
#import sys
import copy
import itertools # scan through frames
import tempfile # temporary output files
import json # save summary file
import numpy as np
import yaml # read Snowstorm config

from icecube import icetray, dataclasses, dataio, clsim
from I3Tray import I3Tray
from icecube.clsim.traysegments.common import setupPropagators, setupDetector, configureOpenCLDevices
from icecube.clsim.traysegments.I3CLSimMakePhotons import I3CLSimMakePhotonsWithServer

from icecube.ice_models import icewave
from icecube.ice_models import angsens_unified
from icecube.snowstorm import Perturber, MultivariateNormal, DeltaDistribution, UniformDistribution
from icecube.snowstorm import all_parametrizations # Snowstorm parametrizations
from icecube import phys_services
print("done")

# set argparse arguments
GCDFile                = args.gcdfile
InputFiles             = args.infile
OutputFile             = args.outfile
RandomService          = phys_services.I3GSLRandomService(args.seed)
NumEventsPerModel      = args.events_per_model
SummaryFile            = args.summaryfile
UseCPUs                = args.cpu
UseGPUs                = not args.cpu
DOMOversizeFactor      = args.domoversizefactor
UseI3PropagatorService = args.UseI3PropagatorService
log_level              = args.log_level

# set icetray logging level
log_levels = {"error" : icetray.I3LogLevel.LOG_ERROR,
              "warn" : icetray.I3LogLevel.LOG_WARN,
              "info" : icetray.I3LogLevel.LOG_INFO,
              "debug" : icetray.I3LogLevel.LOG_DEBUG,
              "trace" : icetray.I3LogLevel.LOG_TRACE}
if log_level.lower() in log_levels.keys():
    icetray.set_log_level(log_levels[log_level.lower()])

# read either yaml or json config file for Snowstorm
print("Reading Snowstorm config... ", end="")
_, ext = os.path.splitext(args.config_file)
if ext == ".json":
    loader = json.load
elif ext in [".yaml", ".yml"]:
    loader = yaml.load
else:
    raise ValueError("Cannot load snowstorm config with extension {}".format(ext))

Snowstorm_config = loader(open(args.config_file))
print("done")

#
#
# Helper Functions
#
#

class FrameSequenceReader(icetray.I3Module):
    """
    Emit frames from an externally supplied dataio.I3FrameSequence, effectively
    making a persistent I3Reader.
    """
    def __init__(self, ctx):
        super(FrameSequenceReader,self).__init__(ctx)
        self.AddParameter("Sequence", "Iterable of frames to emit", None)
    def Configure(self):
        self._frames = self.GetParameter("Sequence")
    def Process(self):
        # this can run into issues if it's the last one
        try:
            frame = next(self._frames)
            if frame is not None:
                self.PushFrame(frame)
            else:
                self.RequestSuspension()
        except StopIteration:
            self.RequestSuspension()


class Bumper(icetray.I3Module):
    """
    Stop the tray after N Q-frames
    """
    def __init__(self, ctx):
        super(Bumper,self).__init__(ctx)
        self.AddParameter("NumFrames", "", 100)
    def Configure(self):
        self._numframes = self.GetParameter("NumFrames")
        self._count = 0
    def DAQ(self, frame):
        self._count += 1
        if self._count >= self._numframes:
            self.PushFrame(frame)
            self.RequestSuspension()
        else:
            self.PushFrame(frame)


class EnsureSFrame(icetray.I3Module):
    """
    Inject an S frame if none present, and ensure that M frames come after S
    """
    def __init__(self, ctx):
        super(EnsureSFrame,self).__init__(ctx)
        self.AddParameter("Enable", "", True)
    def Configure(self):
        self._disabled = not self.GetParameter("Enable")
        self._mframes = []
    def Process(self):
        frame = self.PopFrame()
        if self._disabled:
            self.PushFrame(frame)
            return
        elif frame.Stop.id == 'S':
            # got an existing S frame, emit buffered M frames
            self._disabled = True
            self.PushFrame(frame)
            for m in self._mframes:
                self.PushFrame(m)
            del self._mframes[:]
        elif frame.Stop.id == 'M':
            self._mframes.append(frame)
        elif frame.Stop.id == 'Q':
            # no S frame seen, emit SMQ
            self._disabled = True
            self.PushFrame(icetray.I3Frame('S'))
            for m in self._mframes:
                self.PushFrame(m)
            del self._mframes[:]
            self.PushFrame(frame)
        else:
            self.PushFrame(frame)


class GatherStatistics(icetray.I3Module):
    """Mimick the summary stage of I3CLSimModule::Finish()"""
    def Finish(self):
        if not 'I3SummaryService' in self.context:
            return
        summary = self.context['I3SummaryService']
        server = self.context['CLSimServer']

        if not "TotalNumGeneratedHits" in summary.keys():
            summary["TotalNumGeneratedHits"] = 0
        for k, v in summary.items():
            if k.startswith("I3PhotonToMCPEConverter") and k.endswith("NumGeneratedHits"):
                summary["TotalNumGeneratedHits"] += v
                summary.pop(k)

        for k, v in server.GetStatistics().items():
            if k in summary and (k.startswith('Total') or k.startswith('Num')):
                summary[k] += v
            else:
                summary[k] = v

#
#
# Define simulation parameters
#
#

UseOnlyDeviceNumber     = None
MCTreeName              = "I3MCTree"
OutputMCTreeName        = None
FlasherInfoVectName     = None
FlasherPulseSeriesName  = None
PhotonSeriesName        = None # "PhotonSeriesMap" # set a name to keep all photons (None means keeping MCPEs only)
MCPESeriesName          = "I3MCPESeriesMap"
IceModelLocation        = os.path.expandvars(Snowstorm_config["IceModelLocation"])
DisableTilt             = False
UnWeightedPhotons       = False
UnWeightedPhotonsScalingFactor= None
UseGeant4               = False
ParticleHistory         = True
ParticleHistoryGranularity= 20*icetray.I3Units.m
CrossoverEnergyEM       = None
CrossoverEnergyHadron   = None
UseCascadeExtension     = True
StopDetectedPhotons     = True
PhotonHistoryEntries    = 0
DoNotParallelize        = False
UnshadowedFraction      = 1.0
HoleIceParameterization = os.path.expandvars(Snowstorm_config["HoleIceParameterization"])
WavelengthAcceptance    = None
DOMRadius               = 0.16510*icetray.I3Units.m # 13" diameter
CableOrientation        = None
OverrideApproximateNumberOfWorkItems= None
IgnoreSubdetectors      = ["IceTop"]
ExtraArgumentsToI3CLSimClientModule= dict()

"""
Setup and run Snowstorm (aka MultiSim) by running a series of short
trays, each with a different ice model. This works by front-loading as much
of the expensive initialization (reading the GCD file, setting up
PROPOSAL/Geant4, etc) as possible, so that only the propagation kernel
needs to be recompiled for every tray.
"""

# instantiate baseline detector setup.
# this will help construct the baseline characteristics before applying the perturbers

print("Setting up detector... ", end="")
clsimParams = setupDetector(
    GCDFile=GCDFile,
    SimulateFlashers=bool(FlasherInfoVectName or FlasherPulseSeriesName),
    IceModelLocation=IceModelLocation,
    DisableTilt=DisableTilt,
    UnWeightedPhotons=UnWeightedPhotons,
    UnWeightedPhotonsScalingFactor=UnWeightedPhotonsScalingFactor,
    UseI3PropagatorService=UseI3PropagatorService,
    UseGeant4=UseGeant4,
    CrossoverEnergyEM=CrossoverEnergyEM,
    CrossoverEnergyHadron=CrossoverEnergyHadron,
    UseCascadeExtension=UseCascadeExtension,
    StopDetectedPhotons=StopDetectedPhotons,
    DOMOversizeFactor=DOMOversizeFactor,
    UnshadowedFraction=UnshadowedFraction,
    HoleIceParameterization=HoleIceParameterization,
    WavelengthAcceptance=WavelengthAcceptance,
    DOMRadius=DOMRadius,
    CableOrientation=CableOrientation,
    IgnoreSubdetectors=IgnoreSubdetectors,
)
print("done")

print("Setting up OpenCLDevices... ", end="")
openCLDevices = configureOpenCLDevices(
    UseGPUs=UseGPUs,
    UseCPUs=UseCPUs,
    OverrideApproximateNumberOfWorkItems=OverrideApproximateNumberOfWorkItems,
    DoNotParallelize=DoNotParallelize,
    UseOnlyDeviceNumber=UseOnlyDeviceNumber)
print("done")

#
#
# Setup perturbations
#
#

# create empty "perturber" object
perturber = Perturber()

# get perturbation_cfg dict to simplify calls
perturbation_cfg = Snowstorm_config["Perturbations"]

# loop over all perturbations in the perturbation_cfg
print("Setting up perturbers... ")
for name, params in perturbation_cfg.items():
    # catch special case of IceWavePlusModes
    if name == "IceWavePlusModes":
        if not params["apply"]:
            continue
        if params["type"] == "default":
            print("-> adding {} of type {}".format(name, params["type"]))
            perturber.add('IceWavePlusModes', *icewave.get_default_perturbation())
            continue
        else:
            raise NotImplementedError("IceWavePlusModes of type '{}' are not implemented (yet).".format(params["type"]))

    # all other cases
    if params["type"] == "delta":
        print("-> adding {} of type {}".format(name, params["type"]))
        params = params["delta"]
        perturber.add(name, all_parametrizations[name],
                      DeltaDistribution(params["x0"]))
    elif params["type"] == "gauss":
        print("-> adding {} of type {}".format(name, params["type"]))
        params = params["gauss"]
        # Caution: MultivariateNormal expect the covariance matrix as first argument, so we need to use sigma**2
        perturber.add(name, all_parametrizations[name],
                      MultivariateNormal(dataclasses.I3Matrix(np.diag(params["sigma"])**2), params["mu"]))
    elif params["type"] == "uniform":
        print("-> adding {} of type {}".format(name, params["type"]))
        params = params["uniform"]
        perturber.add(name, all_parametrizations[name],
                      UniformDistribution([dataclasses.make_pair(*limits) for limits in params["limits"]]))
    else:
        raise NotImplementedError("Perturbation '{}' of type '{}' not implemented.".format(name, params["type"]))
print("done")


# Setting up some other things
gcdFrames = list(dataio.I3File(GCDFile))
inputStream = dataio.I3FrameSequence([InputFiles])
summary = dataclasses.I3MapStringDouble()
intermediateOutputFiles = []

#
#
# Run PhotonProp
#
#

# start a model counter
model_counter = 0

# Execute photon propagation
print("Executing photon propagation...", end="")
while inputStream.more():
    # measure CLSimInit time
    time_CLSimInit_start = time.time()
    
    tray = I3Tray()
    tray.context['I3RandomService'] = RandomService
    tray.context['I3SummaryService'] = summary

    # make a mutable copy of the config dict
    config = dict(clsimParams)

    # populate the M frame with I3FrameObjects from clsimParams
    model = icetray.I3Frame('M')
    for k,v in config.items():
        if isinstance(v, icetray.I3FrameObject):
            model[k] = v

    # add EventsPerModel to M-frame
    model["SnowstormEventsPerModel"] = dataclasses.I3UInt64(NumEventsPerModel)

    # apply perturbations in the order they were configured
    perturber.perturb(RandomService, model)

    # check for items in the M-frame that were changed/added by the perturbers
    for k in model.keys():
        if k.startswith('Snowstorm'):
            # keep all Snowstorm keys
            continue
        if not k in config:
            raise KeyError("\n {} was put in the M frame, but does not appear in the CLSim configuration dict".format(k))
        if config[k] != model[k]:
            # if an items was changed, copy it back to clsimParams
            config[k] = model[k]
        else:
            # remove unmodified items from the M frame
            del model[k]

    # add "persistent" I3Reader
    tray.Add(FrameSequenceReader, Sequence=itertools.chain(gcdFrames, [model], inputStream))

    # inject an S frame if it doesn't exist
    tray.Add(EnsureSFrame, Enable=len(intermediateOutputFiles)==0)

    # write pertubations to frame
    def populate_s_frame(frame):
        perturber.to_frame(frame)
    tray.Add(populate_s_frame, Streams=[icetray.I3Frame.Stream('S')])

    # Add Bumper to stop the tray after NumEventsPerModel Q-frames
    tray.Add(Bumper, NumFrames=NumEventsPerModel)

    # initialize CLSim server and setup the propagators
    converters = setupPropagators(RandomService, config,
        UseGPUs=UseGPUs,
        UseCPUs=UseCPUs,
        OverrideApproximateNumberOfWorkItems=OverrideApproximateNumberOfWorkItems,
        DoNotParallelize=DoNotParallelize,
        UseOnlyDeviceNumber=UseOnlyDeviceNumber
    )
    server = clsim.I3CLSimServer("tcp://127.0.0.1:*", clsim.I3CLSimStepToPhotonConverterSeries(converters))
    server_location = server.GetAddress()

    # stash server instance in the context to keep it alive
    tray.context['CLSimServer'] = server

    # recycle StepGenerator to prevent repeated, expensive initialization
    if 'StepGenerator' in ExtraArgumentsToI3CLSimClientModule:
        stepGenerator = ExtraArgumentsToI3CLSimClientModule['StepGenerator']
        stepGenerator.SetMediumProperties(config['MediumProperties'])
        stepGenerator.SetWlenBias(config['WavelengthGenerationBias'])

    # add CLSim server to tray
    module_config = \
        tray.Add(I3CLSimMakePhotonsWithServer,
            ServerAddress=server_location,
            DetectorSettings=config,
            MCTreeName=MCTreeName,
            OutputMCTreeName=OutputMCTreeName,
            FlasherInfoVectName=FlasherInfoVectName,
            FlasherPulseSeriesName=FlasherPulseSeriesName,
            PhotonSeriesName=PhotonSeriesName,
            MCPESeriesName=MCPESeriesName,
            RandomService=RandomService,
            ParticleHistory=ParticleHistory,
            ParticleHistoryGranularity=ParticleHistoryGranularity,
            ExtraArgumentsToI3CLSimClientModule=ExtraArgumentsToI3CLSimClientModule
        )

    # recycle StepGenerator to prevent repeated, expensive initialization
    ExtraArgumentsToI3CLSimClientModule['StepGenerator'] = module_config['StepGenerator']

    # write to temporary output file
    intermediateOutputFiles.append(tempfile.mkstemp(suffix=(OutputFile.split("/"))[-1])[1])
    tray.Add("I3Writer", Filename=intermediateOutputFiles[-1],
            DropOrphanStreams=[icetray.I3Frame.TrayInfo],
            Streams=[icetray.I3Frame.TrayInfo,
                     icetray.I3Frame.Simulation,
                     icetray.I3Frame.Stream('M'),
                     icetray.I3Frame.DAQ,
                     icetray.I3Frame.Physics])

    # gather statistics in the "I3SummaryService"
    tray.Add(GatherStatistics)

    # measure CLSimInit time
    time_CLSimInit = time.time() - time_CLSimInit_start
    summary["CLSimInitTime_{:03d}".format(model_counter)] = time_CLSimInit
    if not "TotalCLSimInitTime" in summary:
        summary["TotalCLSimInitTime"] = time_CLSimInit
    else:
        summary["TotalCLSimInitTime"] += time_CLSimInit

    # measure CLSimTray time
    time_CLSimTray_start = time.time()
    
    # Execute Tray
    tray.Execute()
    
    # measure CLSimTray time
    time_CLSimTray = time.time() - time_CLSimTray_start
    summary["CLSimTrayTime_{:03d}".format(model_counter)] = time_CLSimTray
    if not "TotalCLSimTrayTime" in summary:
        summary["TotalCLSimTrayTime"] = time_CLSimTray
    else:
        summary["TotalCLSimTrayTime"] += time_CLSimTray
    
    # increase model counter
    model_counter += 1
print("done")

# Add number of models to summary
summary["TotalNumberOfModels"] = model_counter

# Concatenate intermediate files
print("Concatenating temporary files... ", end='')
tray = I3Tray()
tray.Add(dataio.I3Reader, "I3Reader", FilenameList=intermediateOutputFiles)
tray.Add("I3Writer", Filename=OutputFile,
         DropOrphanStreams=[icetray.I3Frame.TrayInfo],
         Streams=[icetray.I3Frame.TrayInfo,
                  icetray.I3Frame.Simulation,
                  icetray.I3Frame.Stream('M'),
                  icetray.I3Frame.DAQ,
                  icetray.I3Frame.Physics])
tray.Execute()
tray.Finish()
print("done")

print("Cleaning up Temporary files... ")
for fname in intermediateOutputFiles:
    os.unlink(fname)
print("done")

# Recalculate averages
print("Writing summary file... ", end='')
if not args.cpu:
    if summary['TotalHostTime'] > 0.0:
        summary['DeviceUtilization'] = summary['TotalDeviceTime']/summary['TotalHostTime']
    if summary['TotalNumPhotonsGenerated'] > 0.0:
        summary['AverageDeviceTimePerPhoton'] = summary['TotalDeviceTime']/summary['TotalNumPhotonsGenerated']
    if summary['TotalNumPhotonsGenerated'] > 0.0:
        summary['AverageHostTimePerPhoton'] = summary['TotalHostTime']/summary['TotalNumPhotonsGenerated']

if SummaryFile:
    with open(SummaryFile, 'w') as f:
        json.dump(dict(summary), f)
print("done")

# Hurray!
print("All finished!")

# say something about the runtime
end_time = time.time()
print("That took "+str(end_time - start_time)+" seconds.")
