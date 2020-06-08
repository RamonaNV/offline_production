#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser(usage="%prog [options] infile outfile")
parser.add_option("--single", dest="single", action="store_true", default=False, help="Fit single-muon energy distribution (dependent on depth and angle) rather than bundle energy distribution (additionally dependent on multiplicity and radius)")
opts, args = parser.parse_args()
try:
	infile, outfile = args
except ValueError:
	infile = "/data/uwa/jvansanten/projects/2012/muongun/corsika/SIBYLL/Hoerandel5/atmod_12.hdf5"
	outfile = "Hoerandel5_atmod12_SIBYLL.%s.fits" % (['bundle_energy', 'single_energy'] % opts.single)

import numpy, tables, dashi, os
from icecube.photospline import spglam as glam
from icecube.photospline import splinefitstable
from utils import load_espec, pad_knots

bias = 50
hes = load_espec(infile, opts.single, bias)

extents = [(e[1], e[-2]) for e in hes._h_binedges]

# construct a pseudo-chi^2 weight such that weight_i*(log(x_i + sigma_i) - log(x_i)) = 1
weights = 1./(numpy.log(numpy.exp(hes.bincontent - bias) + numpy.sqrt(hes._h_squaredweights[hes._h_visiblerange])) - (hes.bincontent - bias))
weights[numpy.logical_not(numpy.isfinite(weights))] = 0
# weights[0,:,:] = 0 # ingore the horizon bin

knots = [
	pad_knots(numpy.linspace(0, 1, 11), 2),  # cos(theta)
	pad_knots(numpy.linspace(1, 3, 11), 2),  # vertical depth [km]
	pad_knots(numpy.linspace(0, 20, 21), 2), # log(energy) [log(E/GeV)]
]
smooth = weights[weights > 0].mean() # make the smoothing term proportional to the weights
order = [2,2,2]
penalties = {2:[smooth/1e3, smooth/1e4, smooth/1e2]}    # Penalize curvature 

if opts.single:
	spline = glam.fit(hes.bincontent,weights,hes._h_bincenters,knots,order,smooth,penalties=penalties)
	spline.bias = bias
else:
	# The depth structure in the energy spectrum is mostly a function of
	# attenuation at the bundle rim, which we now have an extra dimension for.
	# Reduce the knot density in depth accordingly.
	knots[1] = pad_knots(numpy.linspace(1, 3, 7), 2)
	knots = knots[0:2] + [
		pad_knots(numpy.linspace(numpy.sqrt(2), numpy.sqrt(100), 5)**2), # bundle multiplicity
		pad_knots(numpy.linspace(0, numpy.sqrt(250), 11)**2), # distance to shower axis
		
	] + [knots[-1]]
	
	knots = [
		pad_knots(numpy.linspace(0, 1, 7), 2),  # cos(theta)
		pad_knots(numpy.linspace(1, 3, 7), 2),  # vertical depth [km]
		pad_knots(numpy.linspace(0, numpy.sqrt(100), 5)**2, 2), # bundle multiplicity
		pad_knots(numpy.linspace(0, 250**(1./2), 11)**2), # distance to shower axis
		pad_knots(numpy.linspace(0, 20, 11), 2), # log(energy) [log(E/GeV)]

	]
	
	order = order[0:2] + [2,2] + [order[-1]]
	penalties[2] = penalties[2][0:2] + [smooth/1e4, smooth/1e1] + [penalties[2][-1]]
	penalties = {2:[smooth/1e5, smooth/1e5, smooth*1e4, smooth/1e3, smooth/1e1]} # Penalize curvature
	# penalties[2][0] *= 1e1
	# penalties[2][1] *= 1e4
	# penalties[2][-1] *= 1e1

	spline = glam.fit(hes.bincontent,weights,hes._h_bincenters,knots,order,smooth,penalties=penalties)
	spline.bias = bias


if os.path.exists(outfile):
	os.unlink(outfile)
splinefitstable.write(spline, outfile)

