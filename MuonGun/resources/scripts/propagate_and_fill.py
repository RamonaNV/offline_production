#!/usr/bin/env python

# specialization for 5-component dCORSIKA

from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
from icecube import MuonGun
from I3Tray import I3Tray
from icecube.weighting.weighting import GenerationProbability, PowerLaw, EnergyWeight, FiveComponent
from icecube.weighting import fluxes
import numpy
import os, subprocess, operator
from optparse import OptionParser
parser = OptionParser(option_class=GenerationProbability.option())
parser.add_option("--flux", dest="flux", default="Hoerandel5", choices=("GaisserH3a", "GaisserH4a", "Hoerandel5", "CascadeOptimized", "5Comp"), help="CR flux model to weight to")
parser.add_option("--oversample", dest="oversample", default=1, type=int, metavar="N", help="Sample each air shower N times")
parser.add_option("--detcfg", dest="detcfg", type=float, default=1, help="height/diameter ratio of target volume. If unspecified, assume a uniform flux over all solid angle.")
parser.add_option("--mindepth", dest="mindepth", type=float, default=1.0, help="minimum vertical depth to propagate to [%default km]")
parser.add_option("--maxdepth", dest="maxdepth", type=float, default=2.8, help="maximum vertical depth to propagate to [%default km]")
parser.add_option("--steps", dest="steps", type=int, default=19, help="number of vertical depths [%default]")
parser.add_option("--n-standard", dest="n_standard", type=int, default=1, help="number of standard 5-comp files")
parser.add_option("--n-he", dest="n_he", type=int, default=1, help="number of high-threshold 5-comp files")

opts, args = parser.parse_args()
opts.n_standard *= opts.oversample
opts.n_he *= opts.oversample
if len(args) < 2:
	parser.error("You must supply at least one input and output file")
infiles, outfile = args[:-1], args[-1]

# If running in Condor, stage files in and out with block copies
if '_CONDOR_SCRATCH_DIR' in os.environ:
	cwd = os.environ['_CONDOR_SCRATCH_DIR']
	subprocess.call(['cp'] + infiles + [cwd])
	infiles = [os.path.join(cwd, os.path.basename(f)) for f in infiles]
	outdir = os.path.dirname(outfile)
	outfile = os.path.join(cwd, os.path.basename(outfile))

ptype = dataclasses.I3Particle.ParticleType
elements = [('H', ptype.PPlus), ('He', ptype.He4Nucleus), ('N', ptype.N14Nucleus), ('Al', ptype.Al27Nucleus), ('Fe', ptype.Fe56Nucleus)]

from utils import IsotropicWeight, VolumeCorrWeight, MultiFiller

# Set up (common) normalization term
standard = FiveComponent(nevents=opts.n_standard*2.5e6, emin=6e2, emax=1e11)
he = FiveComponent(nevents=opts.n_he*1e6, emin=5e4, emax=1e11, normalization=[3., 2., 1., 1., 1.], gamma=[-2.6]*5)
generator = standard + he
if opts.detcfg:
	spectrum = VolumeCorrWeight(opts.detcfg, generator)
else:
	spectrum = IsotropicWeight(generator)

def make_pseudoflux(kind):
	# Weight to a pseudo-flux that integrates to NEvents/(Area*SolidAngle)
	r = 800
	l = 1600
	area = numpy.pi**2*r*(r+l)
	nevents = 1./area
	if kind == "CascadeOptimized":
		flux_components = FiveComponent(gamma=[-2.6]*5, normalization=[3., 2., 1., 1., 1.], emin=3e4, emax=1e9, nevents=nevents)
	elif kind == "Standard":
		flux_components = FiveComponent(nevents=nevents, emin=6e2, emax=1e11)
	return flux_components

# Set up fluxes to weight to
target_fluxes = dict()
for k in "GaisserH3a", "GaisserH4a", "Hoerandel5":
	target_fluxes[k] = getattr(fluxes, k)()
for k in "CascadeOptimized", "Standard":
	target_fluxes[k+"5Comp"] = make_pseudoflux(k)

tray = I3Tray()

from os.path import expandvars
icetray.load('ucr-icetray')
ucr_opts = expandvars('$I3_BUILD/bin/ucr-icetray-ucr ')
for fn in infiles:
	ucr_opts += fn + ' '
ucr_opts += ("-DEPTH=1950 -LENGTH=1600 -RADIUS=800 -over=%d" % (opts.oversample))
print('ucr_opts: "%s"' % ucr_opts)
tray.AddModule('I3InfiniteSource', 'driver')
tray.AddModule('I3GeneratorUCR', 'reader', EventsToIssue=int(1e9), UCROpts=ucr_opts)

from icecube.MuonGun import MuonPropagator, Crust, Sphere

import random
seed = random.randint(0, (1<<31) - 1)
print('Seed: %d' % seed)
MuonPropagator.set_seed(seed)

crust = Crust(MuonPropagator("air", ecut=-1, vcut=5e-2, rho=0.673))
crust.add_layer(Sphere(1948, 6374134), MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=0.832))
crust.add_layer(Sphere(1748, 6373934), MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=1.005))

tray.AddModule('Muonitron', 'propagator',
    Depths=list(numpy.linspace(opts.mindepth, opts.maxdepth, opts.steps)*I3Units.km),
    Propagator=MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=1.005),
    Crust=crust,
)

tray.AddModule(MultiFiller, 'filler',
    Outfile=outfile, Fluxes=target_fluxes, GenerationSpectra=spectrum,
    MinDepth=opts.mindepth, MaxDepth=opts.maxdepth, DepthSteps=opts.steps)
	

tray.Execute()


from subprocess import call
call(['h5repack', '-f', 'GZIP=9', outfile, outfile+'.z'])
call(['mv', outfile+'.z', outfile])
if '_CONDOR_SCRATCH_DIR' in os.environ:
	subprocess.call(['cp', outfile, outdir+'/'+os.path.basename(outfile)])
	
