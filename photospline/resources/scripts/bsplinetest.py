import numpy
import Gnuplot

try:
	input = raw_input
except NameError:
	pass

knots = list(range(3,15))
x = numpy.random.uniform(1,10,size=200)
x = numpy.sort(x)

def bspline(knots, x, i, n):
	if n == 0:
		if (x > knots[i] and x < knots[i+1]):
			return 1
		else:
			return 0

	a = (x - knots[i])*bspline(knots, x, i, n-1)/(knots[i+n] - knots[i])
	b = (knots[i+n+1] - x)*bspline(knots, x, i+1, n-1)/(knots[i+n+1] - knots[i+1])
	return a+b

gp = Gnuplot.Gnuplot()
spline = [Gnuplot.Data(x, [bspline(knots, x_ele, n, 2) for x_ele in x], with_="lines") for n in range(0,len(knots)-2-1)]
gp.plot(spline[0],spline[1],spline[2])
input("Press ENTER to continue")
