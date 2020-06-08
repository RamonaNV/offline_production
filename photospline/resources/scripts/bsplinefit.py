import numpy
from icecube.photospline.glam.bspline import *
import Gnuplot

try:
	input = raw_input
except NameError:
	pass

numpts = 500
knots=list(range(-2,30))
order=2

#x1 = numpy.sort(numpy.random.normal(15,4,size=numpts))
x1 = numpy.sort(numpy.random.uniform(0,25,size=numpts))
z = numpy.random.poisson(numpy.cos(x1)*numpy.cos(x1) + (x1 -3.0)*(x1-3.0) + 10)

A = splinebasis(knots,order,x1)
result = numpy.linalg.lstsq(A, z)

gp = Gnuplot.Gnuplot()
xfine = numpy.sort(numpy.random.uniform(0,25,size=1000))
rawdat = Gnuplot.Data(x1, z)
spline = Gnuplot.Data(xfine, [sum([result[0][n]*bspline(knots, x, n, order) for n in range(0,len(knots)-2-1)]) for x in xfine], with_="lines")
gp.plot(rawdat,spline)
input("Press ENTER to continue")
