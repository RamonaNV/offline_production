#!/usr/bin/env python

import os
import sys

from I3Tray import *
from os.path import expandvars

import argparse
parser = argparse.ArgumentParser(description = "Propagates photons with ppc.")
parser.add_argument('gpu', type = int) 
parser.add_argument('nevents', type = int) 
parser.add_argument('seed', type = int) 
parser.add_argument('outfile', default = 'flashers.i3.bz2') 
args = parser.parse_args()

gcdfile = expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55380_corrected.i3.gz")
os.putenv("PPCTABLESDIR", expandvars("$I3_BUILD/ppc/resources/ice"))
os.putenv("OGPU", "1") # makes sure only GPUs are used (with OpenCL version)

from icecube import icetray, dataclasses, dataio, phys_services, sim_services, ppc

tray = I3Tray()

tray.context["I3RandomService"] = phys_services.I3GSLRandomService(args.seed)

tray.AddModule("I3InfiniteSource", 
               Prefix = gcdfile,
               Stream = icetray.I3Frame.DAQ)

num = args.nevents # how many events

fstr = 63  # flashing string
fdom = 20  # flashing DOM
bri = 127  # flasher brightness (same setting as in flasher xml driver file)
wid = 124  # flasher width (same setting as in flasher xml driver file)
ori = 0    # azimuthal direction (in degrees) of the first LED
mask = 0b000000000001 # bitmask controlling which LEDs are turned on, 0 for off, 1 for on

nled = bin(mask).count("1")

# the average flasher photon yield derived as part of SPIE Lea.
nph = 2.75*2.5e9  # fitted p_y in SPICE Lea * photon bunch
nph /= 0.1315     # DOM acceptance
nph /= 0.9*0.85   # shadowing * disc. threshold loss
nph /= 6*0.977316 # number of LEDs * correction for bri=127 wid=124 (as used in SPICE)
# do not modify the above lines unless you think they contain an error!

nph *= (0.0006753+0.00005593*bri)*(wid+13.9-57.5/(1+bri/34.4))

from icecube.simprod.segments.GenerateFlashers import GenerateFlashers
tray.AddModule(GenerateFlashers,"GenerateFlashers",
               Streams=[icetray.I3Frame.DAQ],
	       FlashString=fstr,
	       FlashDOM=fdom,
	       FlashBrightness=bri,
	       FlashWidth=wid,
	       FlashMask=mask)

# Set FLDR=x+(n-1)*360, where 0<=x<360 and n>0 to simulate n LEDs in a
# symmetrical n-fold pattern, with first LED centered in the direction x.
# Negative or unset FLDR simulates a symmetric in azimuth pattern of light.

fldr = ori+(nled-1)*360   # set fldr=-1 for azimuthally-symmetric emission

os.putenv("FLDR", str(fldr))   # direction of the first flasher LED
os.putenv("WFLA", "405")  # flasher wavelength; set to 337 for standard candles

nph *= nled   # number of LEDs
nph *= 0.1315 # DOM acceptance at 405 nm. For SC use eff(337)=0.0354
nph *= 0.9    # shadowing losses. Set efficiency correction in cfg.txt to 1.

if (int(mask)<64):
    fstr=-fstr

tray.AddModule("i3ppc",
               gpu = args.gpu,
               nph = nph,
               wid = wid*0.5*I3Units.ns,
               fla = OMKey(fstr, fdom)) # set fstr=-fstr for tilted flashers, fstr=0 and fdom=1,2 for SC1 and 2

tray.AddModule("I3Writer", 
               streams = [icetray.I3Frame.DAQ],
               filename = args.outfile)

if(args.nevents):
    tray.Execute(args.nevents + 3)
else:
    tray.Execute()

