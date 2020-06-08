import numpy
try:
	from icecube.photospline import spglam as glam
except ImportError:
	print("Could not get SPGLAM, falling back")
	from icecube.photospline.glam import glam
import Gnuplot
from icecube.photospline.glam.bspline import *

try:
	input = raw_input
except NameError:
	pass

numpts = 500
knots=[numpy.linspace(-8,35,30)]
order=2
smooth=3.14159

#x1 = numpy.sort(numpy.random.normal(15,4,size=numpts))
x1 = numpy.sort(numpy.random.uniform(-4,25,size=numpts))
z = numpy.random.poisson(numpy.cos(x1)*numpy.cos(x1) + (x1 -3.0)*(x1-3.0) + 10)

result = glam.fit(z,1.+z,[x1],knots,2,smooth,[0])

gp = Gnuplot.Gnuplot()
xfine = numpy.sort(numpy.random.uniform(-5,35,size=1000))
rawdat = Gnuplot.Data(x1, z, title="Data")
fitdat = Gnuplot.Data(x1, glam.grideval(result, [x1]),title="Fit")
spline = Gnuplot.Data(xfine, [sum([result.coefficients[n] *
    bspline(knots[0], x, n, order) for n in range(0,len(knots[0])-2-1)])
    for x in xfine], with_="lines",title="Spline")
gp.plot(rawdat,fitdat,spline)
input("Press ENTER to continue")
