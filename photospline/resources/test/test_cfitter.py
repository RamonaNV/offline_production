#!/usr/bin/env python

import numpy
from icecube.photospline import spglam as glam
from icecube.photospline.utils import pad_knots

# reproducibility considered good
numpy.random.seed(42)

z, edges = numpy.histogramdd(numpy.random.multivariate_normal([0,0,0], numpy.eye(3), size=1000),
    bins=[numpy.linspace(-3, 3, 101)]*3)
z = z.cumsum(axis=2).astype(float)
w = numpy.ones(z.shape)

# evaluate CDF at right-hand bin edges
centers = [0.5*(edges[i][:-1] + edges[i][1:]) for i in range(2)] + [edges[-1][1:]]
knots = [pad_knots(numpy.linspace(-3, 3, 5))]*3
order = [2,2,3]
smooth = 1.

spline = glam.fit(z, w, centers, knots, order, smooth, penalties={2:[smooth]*3}, monodim=2)

y = glam.grideval(spline, centers)
residual = ((z-y)**2).sum()
numpy.testing.assert_almost_equal(residual, 50791.31, 2)
