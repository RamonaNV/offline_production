#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] infile outfile")
parser.add_option("--single", dest="single", default=False, action="store_true")
opts, args = parser.parse_args()
if len(args) != 2:
	parser.error("You must supply exactly one input table and one output table name!")
infile, outfile = args

import numpy, tables, dashi, os
from icecube.photospline import spglam as glam
from icecube.photospline import splinefitstable
from utils import load_flux, pad_knots

bias=50
h = load_flux(infile, transform=True, bias=bias)
if opts.single:
	h = h[:,:,1]

extents = [(e[1], e[-2]) for e in h._h_binedges]

# construct a pseudo-chi^2 weight such that weight_i*(log(x_i + sigma_i) - log(x_i)) = 1
weights = 1./(numpy.log(numpy.exp(h.bincontent - bias) + numpy.sqrt(h._h_squaredweights[h._h_visiblerange])) - (h.bincontent - bias))
weights[numpy.logical_not(numpy.isfinite(weights))] = 0
if opts.single:
	pass
	# weights[:,:3] = 0 # ingore shallowest depths: edge effects from energy cut
else:
	weights[:,:,0] = 0 # ignore single-muon bin (fit separately)
	# weights[85:,:,:] = 0 # ignore the horizon bin
	weights[:,:3,:] = 0

knots = [
	pad_knots(numpy.linspace(0, 1, 21)**2, 2),  # cos(theta)
	pad_knots(numpy.linspace(1, 3, 11), 2),  # vertical depth [km]
	pad_knots(numpy.linspace(1, (100)**(1./2), 21)**2), # bundle multiplicity
]
smooth = weights[weights > 0].mean() # make the smoothing term proportional to the weights
order = [2,2,2]
penalties = {2:[smooth/1e12, smooth/1e8, smooth]} # Penalize curvature

if opts.single:
	knots = knots[:-1]
	order = order[:-1]
	penalties[2] = penalties[2][:-1]

spline = glam.fit(h.bincontent,weights,h._h_bincenters,knots,order,smooth,penalties=penalties)
spline.bias = bias

if os.path.exists(outfile):
	os.unlink(outfile)
splinefitstable.write(spline, outfile)

# def eval_slice(h, spline, idx):
# 	from icecube.photospline.glam.glam import grideval
# 	coords = []
# 	x = None
# 	for i, centers in zip(idx, h._h_bincenters):
# 		if isinstance(i, int):
# 			axis = centers[i-1]
# 		else:
# 			axis = centers[i]
# 		try:
# 			len(axis)
# 			axis = numpy.linspace(axis[0], axis[-1], 501)
# 			x = axis
# 		except TypeError:
# 			axis = [axis]
# 		coords.append(axis)
# 	return x, grideval(spline, coords).flatten()
# 
# def plot_slice(h, spline, idx, **kwargs):
# 	sub = h[idx]
# 	sub.scatter(**kwargs)
# 	pylab.plot(*eval_slice(h, spline, idx), **kwargs)
# 
# 
# 
# dashi.visual()
# import pylab
# pylab.figure()
# plot_slice(h, spline, (slice(None),4))
# plot_slice(h, spline, (slice(None),10))
# plot_slice(h, spline, (slice(None),15))
# plot_slice(h, spline, (slice(None),18))
# pylab.show()

