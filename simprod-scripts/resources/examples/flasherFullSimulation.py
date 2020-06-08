#!/usr/bin/env python

from __future__ import print_function

import os
import string
from os.path import expandvars

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio, phys_services, sim_services, clsim
from icecube import DOMLauncher
from icecube.simprod import segments
from icecube.sim_services import bad_dom_list_static
from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)

parser.add_option("-o", "--outfile",
                  default="test_flashes.i3", 
                  dest="OUTFILE", 
                  help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("--seed",
                  type="int",
                  default=12344, 
                  dest="SEED", 
                  help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",
                  default=expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55380_corrected.i3.gz"),
                  dest="GCDFILE", 
                  help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("--runnumber", 
                  type="int", 
                  default=1, 
                  dest="RUNNUMBER", 
                  help="The run number for this simulation")
parser.add_option("--numevents", 
                  type="int", 
                  default=5, 
                  dest="NUMEVENTS", 
                  help="The number of simulated flashes")
parser.add_option("--nstreams", 
                  type="int", 
                  default=10000, 
                  dest="NSTREAMS", 
                  help="The nstreams for a random number generator (I3SPRNGRandomService)")
parser.add_option("-s", "--fstr", 
                  type="int", 
                  default=36, 
                  dest="FSTR", 
                  help="flashing string number")
parser.add_option("-d", "--fdom", 
                  type="int", 
                  default=30, 
                  dest="FDOM", 
                  help="flashing dom number")
parser.add_option("-b", "--flasher-brightness", 
                  type="int", 
                  default=25, 
                  dest="FLASHERBRIGHTNESS", 
                  help="Flasher Brightness")
parser.add_option("-w", "--flasher-width", 
                  type="int", 
                  default=30, 
                  dest="FLASHERWIDTH", 
                  help="Flasher Width")
parser.add_option("--flasher-mask", 
                  default=0b000001000000, 
                  dest="FLASHERMASK", 
                  help="Flasher Mask, 0-5 are tilted, 6-11 are horizontal [cDOM LEDs are all horizantal]")
parser.add_option("--flasher-time", 
                  type="float", 
                  default=0, 
                  dest="FLASHERTIME", 
                  help="Flasher time")
parser.add_option("--dom-oversize-factor", 
                  type="float", 
                  default=1, 
                  dest="DOMOVERSIZEFACTOR", 
                  help="DOM Oversize Factor")
parser.add_option("--dom-unshadowed-fraction", 
                  type="float", 
                  default=0.99, 
                  dest="UNSHADOWEDFRACTION", 
                  help="DOM Unshadowed Fraction")
parser.add_option("--ice-model",
                  default=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"),
                  dest="ICEMODEL", 
                  help="Ice Model directory")
parser.add_option("--max-parallel-events", 
                  type="int", 
                  default=10, 
                  dest="MAXPARALLELEVENTS", 
                  help="maximum number of events(==frames) that will be processed in parallel")
parser.add_option("--no-photon-data", 
                  action="store_true", 
                  default=False, 
                  dest="REMOVEPHOTONDATA", 
                  help="Remove I3Photons before writing the output file (in addition to I3MCPEs)")
parser.add_option("--no-mc-hits", 
                  action="store_false", 
                  default=True, dest="KEEPMCHITS", 
                  help="Remove MC Hits")
parser.add_option("--skip-noise-generator", 
                  action="store_true", 
                  default=False, 
                  dest="SKIPNOISEGENERATOR", 
                  help="Skip Noise Generator")
parser.add_option("--no-dynamic-bad-doms", 
                  action="store_false", 
                  default=True, 
                  dest="DYNAMICBADDOMS", 
                  help="Do not retrieve a dynamic bad dom list from the database")
parser.add_option("--use-gpus", 
                  action="store_true", 
                  default=False, 
                  dest="GPUS", 
                  help="Use GPUs")
parser.add_option("--select-gpu", 
                  type="int", 
                  default=0, 
                  dest="GPUNUMBER", 
                  help="Selects GPU")
parser.add_option("--no-unweighted-photons", 
                  action="store_false", 
                  default=True, 
                  dest="UNWEIGHTEDPHOTONS", 
                  help="Use UnWeightedPhotons=False in clsim I3CLSimMakeHits")
parser.add_option("--use-generateflashers", 
                  action="store_true", 
                  default=False, 
                  dest="USEGENERATEFLASHERS", 
                  help="Don't use GenerateFlashers, instead use clsim.FakeFlasherInfoGenerator")
parser.add_option("--skip-calibration", 
                  action="store_false", 
                  default=True, 
                  dest="DOCALIBRATION", 
                  help='''Do not perform calibration (wavedeform, ...). 
                  use this if your meta-project does not have the required projects''')
parser.add_option("--no-gcd-in-outfile",  
                  action="store_false", 
                  default=True, 
                  dest="INCLUDEGCDINOUTFILE", 
                  help="no GCD frames in the output file")
parser.add_option("--no-calibrated-waveforms",  
                  action="store_true", 
                  default=False, 
                  dest="NOCALIBRATEDWAVEFORMS", 
                  help="no Calibrated Waveforms in the output file")
parser.add_option("--use-clsim",  
                  action="store_true", 
                  default=False, 
                  dest="USECLSIM", 
                  help="use ppc instead of clsim")

(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

########################
outdir=""
outfile=None
if options.OUTFILE:
        outfile = options.OUTFILE
        # did the user specify a directory? then use that and auto-generate
        if os.path.isdir(outfile):
            outdir = outfile
            outfile = None
        else:
            outdir, outfile = os.path.split(outfile)

# add a trailing slash to the output directory name if not already there
if outdir and outdir!="":
    if outdir[-1] != "/":
        outdir += "/"

if not outfile:
        # automatically generate the output filename
        infileRootDir, infileRootFile = os.path.split(infileRoot)
        outfile = infileRootFile + "_calib"
        outfile = outfile + infileExt
print("output dir is %s" % outdir)
print("output file is %s" % outdir + outfile)

########################


mask = options.FLASHERMASK
nled=bin(mask).count("1")

tray = I3Tray()

# set up a random number generator
randomService = phys_services.I3SPRNGRandomService(seed = options.SEED, 
                                                   nstreams = options.NSTREAMS, 
                                                   streamnum = options.RUNNUMBER)

# re-use the same RNG for modules that need it on the context
tray.context['I3RandomService'] = randomService

tray.AddModule("I3InfiniteSource",
               Prefix=options.GCDFILE, 
               Stream=icetray.I3Frame.DAQ)

if options.USEGENERATEFLASHERS:
    print("segments.GenerateFlashers")
    tray.AddModule(segments.GenerateFlashers,
                   Streams=[icetray.I3Frame.DAQ],
                   FlashString=options.FSTR,
                   FlashDOM=options.FDOM,
                   FlashBrightness=options.FLASHERBRIGHTNESS,
                   FlashWidth=options.FLASHERWIDTH,
                   FlashMask=options.FLASHERMASK)
else:
    print("clsim.FakeFlasherInfoGenerator")
    tray.AddModule(clsim.FakeFlasherInfoGenerator, 
                   FlashingDOM = icetray.OMKey(options.FSTR,options.FDOM), 
                   FlasherTime = options.FLASHERTIME*I3Units.ns, 
                   FlasherMask = options.FLASHERMASK, 
                   FlasherBrightness = options.FLASHERBRIGHTNESS, 
                   FlasherWidth = options.FLASHERWIDTH)

photonSeriesName = None
if not options.REMOVEPHOTONDATA:
    photonSeriesName = "PropagatedPhotons"
print("photonSeriesName = %s" % photonSeriesName)

if options.GPUS:
    print("Using GPUs")
    setGPUs = True
    setCPUs = False
else:
    print("Using CPUs")
    setGPUs = False
    setCPUs = True


if options.USECLSIM:
    print("CLSIM mode")
    tray.AddSegment(clsim.I3CLSimMakeHits, 
                    PhotonSeriesName = photonSeriesName, 
                    ParallelEvents = options.MAXPARALLELEVENTS, 
                    RandomService = randomService, 
                    UseGPUs=setGPUs, 
                    UseCPUs=setCPUs,
                    UnWeightedPhotons=options.UNWEIGHTEDPHOTONS, 
                    DOMOversizeFactor=options.DOMOVERSIZEFACTOR, 
                    UnshadowedFraction=options.UNSHADOWEDFRACTION, 
                    IceModelLocation=options.ICEMODEL, 
                    FlasherInfoVectName="I3FlasherInfo")
else:
    print("PPC mode, use GPUs")
    os.putenv("PPCTABLESDIR", expandvars("$I3_BUILD/ppc/resources/ice/lea"))
    os.putenv("OGPU", "1") # makes sure only GPUs are used (with OpenCL version)
    icetray.load("libxppc")
    icetray.load("libppc")
    
    fstr=options.FSTR
    nph=2.75*2.5e9  # fitted p_y in SPICE Lea * photon bunch
    nph/=0.1315     # DOM acceptance
    nph/=0.9*0.85   # shadowing * disc. threshold loss
    nph/=6*0.977316 # number of LEDs * correction for bri=127 wid=124 (as used in SPICE)
    # do not modify the above lines unless you think they contain an error!
    
    nph*=(0.0006753+0.00005593*options.FLASHERBRIGHTNESS)*(options.FLASHERWIDTH+13.9-57.5/(1+options.FLASHERBRIGHTNESS/34.4))
    
    nph*=nled
    nph*=0.1315
    nph/=(options.DOMOVERSIZEFACTOR*options.DOMOVERSIZEFACTOR) # DOM oversize factor squared; must match the value in cfg.txt
    nph*=options.UNSHADOWEDFRACTION # shadowing losses. Set efficiency correction in cfg.txt to 1.
    
    if (int(mask)<64):
        fstr=-fstr
    
    tray.AddModule("i3ppc", 
                   gpu = options.GPUNUMBER, 
                   nph = nph, 
                   wid = options.FLASHERWIDTH*0.5*I3Units.ns,
                   fla = icetray.OMKey(fstr, options.FDOM)) 
    # set fstr=-fstr for tilted flashers, fstr=0 and fdom=1,2 for SC1 and 2


tray.AddSegment(segments.DetectorSim, 
                RandomService = "I3RandomService", 
                GCDFile = options.GCDFILE, 
                InputPESeriesMapName = "MCPESeriesMap", 
                KeepMCHits = options.KEEPMCHITS, 
                SkipNoiseGenerator = options.SKIPNOISEGENERATOR)

badDOMsSLC = bad_dom_list_static.IC86_static_bad_dom_list()     # SLC only
badDOMsHLC = bad_dom_list_static.IC86_static_bad_dom_list_HLC() # SLC + HLC

badDOMsSLC = set(badDOMsSLC)
badDOMsHLC = set(badDOMsHLC)

badDOMsSLC.add(icetray.OMKey(49,42))
badDOMsHLC.add(icetray.OMKey(49,42))

badDOMsSLC = list(badDOMsSLC)
badDOMsHLC = list(badDOMsHLC)

if options.DOCALIBRATION:
    tray.AddSegment(segments.Calibration, 
                    BadDOMsHLC = badDOMsHLC, 
                    BadDOMsSLC = badDOMsSLC)

    # now we generate P-frames
    tray.AddModule("I3NullSplitter", InputPulseSeries='OfflinePulses')                   
    tray.AddModule("I3LCPulseCleaning",  
                   Input='OfflinePulses', 
                   OutputHLC='UncleanedInIcePulsesHLC', 
                   OutputSLC='UncleanedInIcePulsesSLC')

    # cleanup
    if options.NOCALIBRATEDWAVEFORMS:
       tray.AddModule('Delete', Keys = ["OfflinePulses"])


def check(frame):
    assert('I3MCPulseSeriesMap' in frame)
    if options.DOCALIBRATION:
        assert('UncleanedInIcePulsesHLC' in frame)
        assert('UncleanedInIcePulsesSLC' in frame)
        if not options.NOCALIBRATEDWAVEFORMS:
            assert('OfflinePulses' in frame)

tray.Add(check, Streams=[icetray.I3Frame.Physics])

if options.INCLUDEGCDINOUTFILE:
    tray.AddModule("I3Writer",
                   Filename = outdir+outfile,)
else:
    tray.AddModule("I3Writer",
                   Filename = outfile, 
                   Streams=[icetray.I3Frame.DAQ, 
                            icetray.I3Frame.Physics, 
                            icetray.I3Frame.Stream('S'), 
                            icetray.I3Frame.Stream('M')])

print("executing number of events = %s" % options.NUMEVENTS)
tray.Execute(options.NUMEVENTS+3)

















