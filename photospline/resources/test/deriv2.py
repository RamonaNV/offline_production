#!/usr/bin/env python

import numpy, os
from icecube.photospline import I3SplineTable
from icecube.photospline.splinetable import SplineTable
from icecube.photospline.splinefitstable import write

fname = "constant.fits"

# make a constant spline surface
spline = SplineTable()
spline.ndim = 2
spline.order = [2 for i in range(spline.ndim)]
spline.knots = [numpy.linspace(0, 1, 10) for i in range(spline.ndim)]
nsplines = tuple(knots.size-order-1 for knots,order in zip(spline.knots, spline.order))
spline.coefficients = numpy.ones(nsplines)

if os.path.exists(fname):
	os.unlink(fname)
write(spline, fname)

x = [0.5]*spline.ndim

spline = I3SplineTable(fname)

assert spline.eval(x) == spline.eval(x, [0,0]) == 1.

# all second derivatives must be zero
for i in range(len(x)):
	derivs = [0]*2
	derivs[i] = 2
	assert spline.eval(x, derivs) == 0.

os.unlink(fname)
