#!/usr/bin/env python

from argparse import ArgumentParser
from os.path import expandvars
parser = ArgumentParser()
parser.add_argument("-g", "--gcd", default=expandvars('$I3_TESTDATA/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'))
parser.add_argument("outfile")
args = parser.parse_args()

from icecube import icetray, dataclasses, dataio, phys_services
from I3Tray import I3Tray
from os.path import expandvars

tray = I3Tray()

try:
    randomService = phys_services.I3SPRNGRandomService(1, 10000, 1)
except AttributeError:
    randomService = phys_services.I3GSLRandomService(1)
tray.context['I3RandomService'] = randomService

from icecube.icetray import I3Units
from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw
from icecube.MuonGun.segments import GenerateBundles
from icecube.simprod.segments import PropagateMuons

# Use Hoerandel as a template for generating muons
model = load_model('Hoerandel5_atmod12_SIBYLL')
# Generate only single muons, no bundles
model.flux.max_multiplicity = 1
# Center the sampling surface on the barycenter of IC79 strings
surface = Cylinder(1600*I3Units.m, 800*I3Units.m, dataclasses.I3Position(31.25, 19.64, 0))
surface = Cylinder(1600*I3Units.m, 800*I3Units.m)
# Draw energies from an E^-2 power law broken at 1 TeV, from 10 TeV to 10 PeV
spectrum = OffsetPowerLaw(2, 1*I3Units.TeV, 10*I3Units.TeV, 10*I3Units.PeV)
# Set up the generator. This gets stored in a special frame for later reference
generator = StaticSurfaceInjector(surface, model.flux, spectrum, model.radius)

tray.AddSegment(GenerateBundles, 'MuonGenerator', Generator=generator, NEvents=10, GCDFile=args.gcd)

icetray.logging.set_level_for_unit("I3PropagatorService", "TRACE")
tray.Add("Rename", Keys=["I3MCTree", "I3MCTree_preMuonProp"])
tray.Add(PropagateMuons, RandomService=randomService)

tray.AddModule('I3Writer', 'writer',
    Streams=list(map(icetray.I3Frame.Stream, "SQP")),
    filename=args.outfile)


tray.Execute()

