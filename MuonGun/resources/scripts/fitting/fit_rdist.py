#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser(usage="%prog infile outfile")
opts, args = parser.parse_args()
try:
	infile, outfile = args
except ValueError:
	infile = "/data/uwa/jvansanten/projects/2012/muongun/corsika/SIBYLL/Hoerandel5/atmod_12.hdf5"
	outfile = "Hoerandel5_atmod12_SIBYLL.radius.fits"
print(infile, outfile)

import numpy, tables, dashi, os
from icecube.photospline import spglam as glam
from icecube.photospline import splinefitstable
from utils import load_radial_distribution, pad_knots

bias = 50
h = load_radial_distribution(infile, bias)

extents = [(e[1], e[-2]) for e in h._h_binedges]

# construct a pseudo-chi^2 weight such that weight_i*(log(x_i + sigma_i) - log(x_i)) = 1
weights = 1./(numpy.log(numpy.exp(h.bincontent - bias) + numpy.sqrt(h._h_squaredweights[h._h_visiblerange])) - (h.bincontent - bias))
weights[numpy.logical_not(numpy.isfinite(weights))] = 0
weights[:,:,0,:] = 0 # ignore the single-track bin
# weights[0,:,6,:] = 0 # ignore high-multiplicity bin at the horizon

knots = [
	pad_knots(numpy.linspace(0, 1, 11), 2),  # cos(theta)
	pad_knots(numpy.linspace(1, 3, 11), 2),  # vertical depth [km]
	pad_knots(numpy.linspace(0, numpy.sqrt(100), 15)**2, 2), # bundle multiplicity
	pad_knots(numpy.linspace(0, numpy.sqrt(250), 25)**2, 2), # radius [m]	
]
smooth = weights[weights > 0].mean() # make the smoothing term proportional to the weights
order = [2,2,2,2]
penalties = {2:[smooth/1e5, smooth/1e5, smooth/1e4, smooth/10]} # Penalize curvature

spline = glam.fit(h.bincontent,weights,h._h_bincenters,knots,order,smooth,penalties=penalties)
spline.bias = bias

if os.path.exists(outfile):
	os.unlink(outfile)
splinefitstable.write(spline, outfile)
