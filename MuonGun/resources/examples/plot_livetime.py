#!/usr/bin/env python

"""
Use muon flux weights to calculate an effective livetime for combined
CORSIKA samples as a function of energy.
"""

from argparse import ArgumentParser
from os.path import expandvars
parser = ArgumentParser()
parser.add_argument("outfile", help="save plot to file")
args = parser.parse_args()

import matplotlib
matplotlib.use('agg')
import pylab, numpy
from icecube import dataclasses, MuonGun

surface = MuonGun.Cylinder(1600, 800)
area = numpy.pi**2*surface.radius*(surface.radius+surface.length)

# 1 file of E^-2.6 5-component 3e4-1e9 GeV (3:2:1:1:1)
soft = 4e5*MuonGun.corsika_genprob('CascadeOptimized5Comp')
# 1 file of E^-2 5-component 6e2-1e11 GeV (10:5:3:2:1)
hard = 2.5e6*MuonGun.corsika_genprob('Standard5Comp')
# In order to compare to "unweighted" CORSIKA, turn the Hoerandel flux
# into a probability (since we happen to know the integral)
areanorm = 0.131475115*area
# 1 file of natural-spectrum ("unweighted") CORSIKA
unweighted = (2.5e7/areanorm)*MuonGun.corsika_genprob('Hoerandel5')

model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
model.flux.min_multiplicity = 1
model.flux.max_multiplicity = 1
# spectrum = MuonGun.OffsetPowerLaw(5.0, 5e2, 8e2, 10**4.5)
spectrum = MuonGun.OffsetPowerLaw(5, 8e2, 2e3, 1e5)

# spectrum = MuonGun.OffsetPowerLaw(1.1, 650, 800, 1e8)
gun = 1e5*MuonGun.EnergyDependentSurfaceInjector(surface, model.flux, spectrum, model.radius)
# gun = 1e5*MuonGun.StaticSurfaceInjector(surface, model.flux, spectrum, model.radius)

# gun.target_surface = lambda e: surface


def get_weight(weighter, energy, zenith=numpy.pi/8, scale=True):
	shape = energy.shape
	if scale:
		x = numpy.array([gun.target_surface(e).radius - 1 for e in energy])
	else:
		# x = numpy.ones(shape[0])*surface.radius - 1
		x = numpy.ones(shape[0])*surface.radius - 1
	
	# x = surface.radius*numpy.ones(shape) - 1
	y = numpy.zeros(shape)
	# z = z*numpy.ones(shape)
	if scale:
		z = numpy.array([gun.target_surface(e).center.z + gun.target_surface(e).length/2. for e in energy])
	else:
		z = numpy.ones(shape[0])*(surface.center.z + surface.length/2.)
	
	azimuth = numpy.zeros(shape)
	zenith = zenith*numpy.ones(shape)
	multiplicity = numpy.ones(shape, dtype=numpy.uint32)
	mmax = multiplicity.max()
	e = numpy.zeros(shape + (mmax,), dtype=numpy.float32)
	e[:,0] = energy
	r = numpy.zeros(shape + (mmax,), dtype=numpy.float32)
	try:
		return weighter(x, y, z, zenith, azimuth, multiplicity, e, r)
	except:
		# work around lack of vectorized pybindings
		@numpy.vectorize
		def weight(x,y,z,zenith,azimuth,multiplicity,e,r):
			axis = dataclasses.I3Particle()
			axis.pos = dataclasses.I3Position(x,y,z)
			axis.dir = dataclasses.I3Direction(zenith,azimuth)
			assert multiplicity == 1
			bundle = MuonGun.BundleConfiguration()
			bundle.append(MuonGun.BundleEntry(float(r),float(e)))
			return weighter(axis, bundle)
		return weight(x, y, z, zenith, azimuth, multiplicity, e[:,0], r[:,0])

e = numpy.logspace(1, 7, 101)

target = MuonGun.load_model('Hoerandel5_atmod12_SIBYLL')
# target = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')

generators = [
	('200k $E^{-2}$ 5-component CORSIKA', 2e5*hard),
	('300k $E^{-2.6}$ 5-component CORSIKA', (300e3)*soft),
	('200k unweighted CORSIKA', 2e5*unweighted),
	('1k MuonGun (100k muons each)', 1e3*gun),
	
	('total', 2e5*hard + 300e3*soft + 2e5*unweighted + 1e3*gun),
]


fig = pylab.figure(figsize=(6,4))
fig.subplots_adjust(bottom=0.15)

annum = 365*24*3600
for label, generator in generators:
	weighter = MuonGun.WeightCalculator(target, generator)
	pylab.plot(e, 1./(get_weight(weighter, e, scale=True)*annum), label=label)

pylab.loglog()
pylab.legend(loc='lower right', prop=dict(size='x-small'))
pylab.ylabel('Single-muon livetime [years]')
pylab.xlabel('Muon energy at sampling surface [GeV]')
pylab.grid()

pylab.savefig(args.outfile)
