#!/usr/bin/env python
"""
An example of a full simulation chain.
Note that this is only an example.
"""

import os

#### local settings:
if "_CONDOR_SCRATCH_DIR" in os.environ:
    scratchDirBase = os.environ["_CONDOR_SCRATCH_DIR"]
elif "TMPDIR" in os.environ:
    scratchDirBase = os.environ["TMPDIR"]
else:
    scratchDirBase = "/scratch"

cascadeTablePath="/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines"

if not os.path.isdir(scratchDirBase):
    print " *** scratch directory %s does not exist." % scratchDirBase
    print " *** please change the script code to point \"scratchDirBase\" it to scratch space or set $TMPDIR"
    exit(-1)



from I3Tray import *
from icecube import icetray

from icecube.simprod import segments


from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_full_simulation.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format) [no GCD]")
parser.add_option("--include-gcd-in-outfile",  action="store_true", default=False,
                  dest="INCLUDEGCDINOUTFILE", help="include GCD frames in the output file")

parser.add_option("-s", "--seed",type="int",default=437289,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",default="auto",
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format). Setting this to \"auto\" (the default) will select a GCD file according to the --detector parameter from $I3_TESTDATA")
parser.add_option("--datasetnumber", type="int", default=1,
                  dest="DATASETNUMBER", help="The dataset number for this simulation (runs within a deteset get unique seeds)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")

parser.add_option("-d", "--detector",default="IC86",
                  dest="DETECTOR", help="either IC86 or IC79")


parser.add_option("--generator",default="nugen",choices=('nugen', 'muongun', 'corsika-nu','corsika-icetop'),
                  dest="GENERATOR", help="either nugen or muongun")

# nugen
parser.add_option("-n", "--numevents", type="int", default=100,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("--flavor", type="str", default="NuMu",
                  dest="FLAVOR", help="The neutrino flavor (NuE, NuMu or NuTau)")
parser.add_option("--no-auto-extend-muon-volume", action="store_false", default=True,
                  dest="AUTOEXTEND", help="Auto-extend the muon generation volume with energy")

# muongun

parser.add_option("--low-energy-mode",  action="store_true", default=False,
                  dest="LOWENERGYMODE", help="simulate at lower energies (applies only to MuonGun)")


parser.add_option("--from-energy", type="float", default=1.*I3Units.TeV / I3Units.GeV,
                  dest="FROMENERGY", help="lower neutrino energy bound in GeV")
parser.add_option("--to-energy", type="float", default=10.*I3Units.PeV / I3Units.GeV,
                  dest="TOENERGY", help="upper neutrino energy bound in GeV")

# clsim
parser.add_option("-p", "--max-parallel-events", type="int", default=1,
                  dest="MAXPARALLELEVENTS", help="maximum number of events(==frames) that will be processed in parallel")

parser.add_option("--icemodel",default="SpiceMie",
                  dest="ICEMODEL", help="Either Spice1, SpiceMie or SpiceLea")

parser.add_option("--no-hybrid",  action="store_false", default=True,
                  dest="HYBRIDSIMULATION", help="do not perform a hybrid simulation (i.e. use clsim only)")
parser.add_option("--ignore-muon-light",  action="store_true", default=False,
                  dest="IGNOREMUONLIGHT", help="ignore light from (bare) muons in hybrid mode. (i.e. do not run clsim)")
parser.add_option("--use-gpu",  action="store_true", default=False,
                  dest="USEGPU", help="simulate using GPUs instead of CPU cores")
parser.add_option("--without-photon-prop",  action="store_false", default=True,
                  dest="CLSIM", help="do not propagate photons")

# detector sim
parser.add_option("--keep-mchits", action="store_true", default=False,
                  dest="KEEPMCHITS", help="Keep I3MCHits before writing the output file")
parser.add_option("--with-icetop", action="store_true", default=False,
                  dest="DO_ICETOP", help="Simulate IceTop response in addition to InIce")

# icetop
parser.add_option("--top-response", action="store", default='param',
                  dest="TOPRESPONSE", help="IceTop tank response to use. One of 'param' or 'g4' (g4 is incompatible with CLSim)")

# used in simulation and reconstrcution
parser.add_option("--unshadowed-fraction", type="float", default=0.9,
                  dest="UNSHADOWEDFRACTION", help="Fraction of a DOM that is not shadowed by the cable (used in simulation and reconstruction)")


# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if options.GENERATOR.startswith('corsika'):
    infiles = list(args)
    del args[:]
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)


if not os.path.isdir(cascadeTablePath) and options.HYBRIDSIMULATION:
    print " *** cascade table directory %s not found." % cascadeTablePath
    print " *** please change the script code to point \"cascadeTablePath\" to where your spline tables are"
    exit(-1)


import os
import sys
import shutil

print os.uname()

from icecube import icetray, dataclasses, dataio, phys_services

if options.DETECTOR not in ("IC79", "IC86"):
    raise RuntimeError("unknown detector %s. use IC79 or IC86" % options.DETECTOR)


eventGenerationTimes = {'IC86': dataclasses.I3Time(2011, 158100000000000000),   # 2011-07-02 23:40:00.000,000,000,0 UTC
                        'IC79': dataclasses.I3Time(2010, 158100000000000000)}   # 2010-07-02 23:40:00.000,000,000,0 UTC
eventGenerationTime = eventGenerationTimes[options.DETECTOR]

if options.GCDFILE=="auto":
    gcdFilesDefault = {'IC86': expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz"),
                       'IC79': expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC79.55380_L2a.i3.gz"),}
    gcdFile = gcdFilesDefault[options.DETECTOR]
else:
    gcdFile = options.GCDFILE

if options.GENERATOR.lower() == "nugen":
    scratchDir = "%s/%s_%s_%u_%08u_%010u" % (scratchDirBase, options.DETECTOR, options.FLAVOR, options.DATASETNUMBER, options.RUNNUMBER, options.SEED)
elif options.GENERATOR.lower() == "muongun":
    scratchDir = "%s/%s_MuonGun_%u_%08u_%010u" % (scratchDirBase, options.DETECTOR, options.DATASETNUMBER, options.RUNNUMBER, options.SEED)
elif options.GENERATOR == "corsika-nu":
    scratchDir = "%s/%s_corsika-nu_%u_%08u_%010u" % (scratchDirBase, options.DETECTOR, options.DATASETNUMBER, options.RUNNUMBER, options.SEED)
    if len(infiles) == 0:
        parser.error("You must supply at least 1 CORSIKA file")
elif options.GENERATOR == "corsika-icetop":
    scratchDir = "%s/%s_corsika-icetop_%u_%08u_%010u" % (scratchDirBase, options.DETECTOR, options.DATASETNUMBER, options.RUNNUMBER, options.SEED)
    if len(infiles) == 0:
        parser.error("You must supply at least 1 CORSIKA file")
    # some checks just because everything IceTop is not the default.
    if not options.DO_ICETOP:
        parser.error("If you use the --corsika-icetop, you should also use --with-icetop.")
    if options.TOPRESPONSE == 'g4' and options.CLSIM:
        parser.error("If you use --top-response='g4', you also need to use --no-photon-prop.")
else:
    raise RuntimeError("unknown generator %s" % options.GENERATOR)

if os.path.exists(scratchDir):
    print "WARN: ** scratch dir already exists:", scratchDir
    if not os.path.isdir(scratchDir):
        raise RuntimeError("scratch directory exists and is not a directory")
else:
    os.makedirs(scratchDir)


outfile_basename, outfile_ext = os.path.splitext(options.OUTFILE)
if outfile_ext==".bz2":
    delayed_zip = True
    outfile = scratchDir + "/" + os.path.basename(outfile_basename)
    final_outfile = options.OUTFILE
else:
    delayed_zip = False
    outfile = scratchDir + "/" + os.path.basename(options.OUTFILE)
    final_outfile = options.OUTFILE

print ""
print "dataset", options.DATASETNUMBER
print "run", options.RUNNUMBER
print "number of events", options.NUMEVENTS

print "generator is", options.GENERATOR
if options.GENERATOR.lower() == "nugen":
    print "flavor", options.FLAVOR

if delayed_zip:
    print "writing to %s," % outfile
    print "then zipping to %s" % final_outfile
else:
    print "writing to", outfile
print "detector is", options.DETECTOR
print "using GCD file", gcdFile
print "seed", options.SEED
print "scratch space at", scratchDir
print ""

if options.HYBRIDSIMULATION:
    print "loading spline tables.."
    cascade_tables = segments.LoadCascadeTables(IceModel = options.ICEMODEL, TablePath = cascadeTablePath)
    print "done."
else:
    cascade_tables = None

tray = I3Tray()

maxDatasetNumber=10000
maxRunNumber=100000
maxInternalRunNumber = maxDatasetNumber*maxRunNumber


if options.RUNNUMBER < 0:
    raise RuntimeError("negative run numbers are not supported")
elif options.RUNNUMBER >= maxRunNumber:
    raise RuntimeError("run numbers > %u are not supported" % maxRunNumber)

if options.DATASETNUMBER < 0:
    raise RuntimeError("negative dataset numbers are not supported")
elif options.DATASETNUMBER >= maxDatasetNumber:
    raise RuntimeError("dataset numbers > %u are not supported" % maxDatasetNumber)

internalRunNumber = options.DATASETNUMBER*maxRunNumber + options.RUNNUMBER

# make sure we use the same seed for all RNGs (different streams are okay!)

# set up a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = maxInternalRunNumber*2,
    streamnum = internalRunNumber + maxInternalRunNumber)
# re-use the same RNG for modules that need it on the context
tray.context['I3RandomService'] = randomService

# a special random service for muon propagation
randomServiceForPropagators = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = maxInternalRunNumber*2,
    streamnum = internalRunNumber)

if options.GENERATOR.lower() == "nugen":

    tray.AddModule("I3InfiniteSource","TheSource",
                   Prefix=gcdFile,
                   Stream=icetray.I3Frame.DAQ)    
    tray.AddSegment(segments.GenerateNeutrinos, "GenerateNeutrinos",
        RandomService = randomService,
        NumEvents = options.NUMEVENTS,
        Flavor = options.FLAVOR,
        AutoExtendMuonVolume = options.AUTOEXTEND,
        GammaIndex = 2.,
        
        FromEnergy = options.FROMENERGY*I3Units.GeV,
        ToEnergy   = options.TOENERGY*I3Units.GeV,

        # PointSourceDecl = -29.008*I3Units.deg, # Sgr A*
        # PointSourceRA   = 266.417*I3Units.deg,
        )
elif options.GENERATOR.lower() == "muongun":
    if options.LOWENERGYMODE:
        MuonGunFromEnergy     = 100.*I3Units.GeV
        MuonGunToEnergy       = 10.*I3Units.TeV
    else:
        MuonGunFromEnergy     = 10.*I3Units.TeV
        MuonGunToEnergy       = 10.*I3Units.PeV

    tray.AddSegment(segments.GenerateCosmicRayMuons, "GenerateCosmicRayMuons",
        mctree_name = "I3MCTree_preMuonProp",
        gcd_file = gcdFile,
        num_events = options.NUMEVENTS,

        energy_min     = MuonGunFromEnergy,
        energy_max       = MuonGunToEnergy,

        gamma_index = 1.1,
        )
elif options.GENERATOR.startswith("corsika"):
    if options.GENERATOR == 'corsika-icetop':
        tray.AddSegment(segments.GenerateIceTopShowers, "GenerateIceTopShowers",
                        NSamples = options.NUMEVENTS,
                        Files = infiles,
                        GCDFile = options.GCDFILE,
                        x=0, y=0, r=800*I3Units.m, 
                        TankResponse=options.TOPRESPONSE,
                        RunID=options.RUNNUMBER,
                        RaiseObservationLevel=3*I3Units.m
                    )
    else:
        tray.AddSegment(segments.GenerateAirShowers, "GenerateAtmosphericNeutrinos",
                        NEvents = options.NUMEVENTS,
                        Files = infiles,
                        GCDFile = options.GCDFILE,
                        CylinderRadius = 500, CylinderHeight = 1000,
                        SimulateIceTop = options.DO_ICETOP,
                    )
        if options.GENERATOR == "corsika-nu":
            tray.AddSegment(segments.SelectNeutrino, "SelectNeutrino",
                            AutoExtendMuonVolume = options.AUTOEXTEND,
                            EnergyBiasPower = 1,
                            CylinderRadius = 500, CylinderHeight = 1000,
                        )
    tray.Add('Rename', Keys=['I3MCTree', 'I3MCTree_preMuonProp'])
else:
    raise RuntimeError("Unkown generator %s" % options.GENERATOR)

tray.AddSegment(segments.PropagateMuons, "PropagateMuons",
    RandomService = randomServiceForPropagators)
    
if options.CLSIM:
    tray.AddSegment(segments.PropagatePhotons, "PropagatePhotons",
                    RandomService = randomService,
                    MaxParallelEvents = options.MAXPARALLELEVENTS,
                    KeepIndividualMaps = False,
                    IceModel = options.ICEMODEL,
                    UnshadowedFraction = options.UNSHADOWEDFRACTION,
                    HybridMode = options.HYBRIDSIMULATION,
                    IgnoreMuons = options.IGNOREMUONLIGHT,
                    # OutputPhotonSeriesName = "I3PhotonSeriesMap",
                    #UseAllCPUCores = True,
                    UseGPUs = options.USEGPU,
                    CascadeService = cascade_tables)

if options.DO_ICETOP:
    tray.Add("Rename", Keys=['I3MCPESeriesMap', 'InIceMCPEs'])
    tray.Add("I3CombineMCPE", 
            InputResponses = ["IceTopMCPEs", "InIceMCPEs"],
            OutputResponse = "I3MCPESeriesMap")
    tray.Add("Delete", Keys=['InIceMCPEs', 'IceTopMCPEs'])

tray.AddSegment(segments.DetectorSim, "DetectorSim",
    RandomService = 'I3RandomService',
    RunID = options.RUNNUMBER,
    GCDFile = gcdFile,
    InputPESeriesMapName = "I3MCPESeriesMap",
    KeepMCHits = True,
    SkipNoiseGenerator = False)

if options.INCLUDEGCDINOUTFILE:
    tray.AddModule("I3Writer","writer",
        Filename = outfile,
        )
else:
    tray.AddModule("I3Writer","writer",
        Filename = outfile,
        # DropOrphanStreams=[icetray.I3Frame.DAQ],
        Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Stream('S'), icetray.I3Frame.Stream('M')],
        )



tray.Execute()


del tray

if delayed_zip:
    if os.path.exists(outfile+".bz2"):
        print "removing existing bzip2 output file in scratch directory.."
        os.remove(outfile+".bz2")

    print "bzipping output file.."
    os.system("/usr/bin/bzip2 %s" % (outfile))

    print "copying output to final destination.."
    shutil.copyfile(outfile + ".bz2", final_outfile)

    print "cleaning up.."
    os.remove(outfile + ".bz2")
    os.rmdir(scratchDir)

else:
    print "copying output to final destination.."
    shutil.copyfile(outfile, final_outfile)

    print "cleaning up.."
    os.remove(outfile)
    os.rmdir(scratchDir)

print "done."

